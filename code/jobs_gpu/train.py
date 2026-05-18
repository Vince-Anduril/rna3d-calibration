"""Train RNA structure model on the Stanford 3D Folding dataset.

Two-stage:
    Stage 1 (~30 min): MLM pretraining on v2 train sequences (≈5k unique).
    Stage 2 (~rest of budget): structure supervision on sequences with PDB
                                ground-truth, joint coord + confidence loss.

Designed to saturate a single RTX 5090 (32 GB) for ~3 hours. Auto-checkpoints
every 10 minutes, auto-commits + pushes the run directory contents (excluding
the largest .pt weights, which we save but don't push if oversized).
"""
from __future__ import annotations
import argparse, json, os, time, math, random, subprocess, shutil
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader

import sys
sys.path.insert(0, str(Path(__file__).parent))
from model import RNAStructureModel, ModelConfig
from loss import mlm_loss, structure_loss, confidence_loss
from data import (
    StanfordRNADataset, SyntheticRNADataset, encode,
    PAD_ID, MASK_ID, CLS_ID, make_mlm_batch,
)


def now() -> str:
    return time.strftime("%H:%M:%S")


def log(msg: str, fp=None):
    line = f"[{now()}] {msg}"
    print(line, flush=True)
    if fp is not None:
        fp.write(line + "\n")
        fp.flush()


def git_commit_push(repo_root: Path, msg: str, paths: list[str]):
    """Add the listed paths, commit, and push. Best-effort; never raises."""
    try:
        subprocess.run(["git", "-C", str(repo_root), "add", *paths], check=False)
        r = subprocess.run(
            ["git", "-C", str(repo_root), "commit", "-m", msg],
            capture_output=True, text=True
        )
        if "nothing to commit" not in (r.stdout + r.stderr).lower():
            subprocess.run(["git", "-C", str(repo_root), "push", "origin", "main"], check=False, capture_output=True)
    except Exception as e:
        print(f"[git_commit_push] failed (non-fatal): {e}")


def stage1_pretrain(model, device, csv_path: Path, out_dir: Path, log_fp, deadline: float,
                    max_len: int = 256, batch_size: int = 64):
    """Masked-LM pretraining."""
    log("=== Stage 1: MLM pretraining ===", log_fp)
    df = pd.read_csv(csv_path)
    seqs = df["sequence"].tolist()
    log(f"Loaded {len(seqs)} sequences for MLM pretraining", log_fp)

    # Build a simple iterable
    class _MLMDS(torch.utils.data.Dataset):
        def __len__(self): return len(seqs) * 4  # 4 windows per seq per epoch
        def __getitem__(self, i):
            s = seqs[i % len(seqs)]
            if len(s) > max_len - 1:
                start = random.randint(0, len(s) - (max_len - 1))
                s = s[start : start + max_len - 1]
            ids = encode(s, max_len)
            ids = ids + [PAD_ID] * (max_len - len(ids))
            return torch.tensor(ids, dtype=torch.long)

    loader = DataLoader(_MLMDS(), batch_size=batch_size, shuffle=True, num_workers=4, pin_memory=True, drop_last=True)
    optim = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=0.01)
    scaler = torch.amp.GradScaler("cuda") if device.type == "cuda" else None

    model.train()
    step = 0
    losses = []
    start = time.time()
    stage1_deadline = min(deadline, start + 30 * 60)  # 30 min cap

    it = iter(loader)
    while time.time() < stage1_deadline:
        try:
            ids = next(it)
        except StopIteration:
            it = iter(loader); ids = next(it)
        ids = ids.to(device, non_blocking=True)
        masked, labels = make_mlm_batch(ids)
        optim.zero_grad(set_to_none=True)
        if scaler is not None:
            with torch.amp.autocast("cuda", dtype=torch.bfloat16):
                out = model(masked)
                loss = mlm_loss(out["mlm_logits"], labels)
            scaler.scale(loss).backward()
            scaler.unscale_(optim)
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(optim); scaler.update()
        else:
            out = model(masked)
            loss = mlm_loss(out["mlm_logits"], labels)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optim.step()
        losses.append(float(loss.item()))
        step += 1
        if step % 100 == 0:
            log(f"stage1 step={step} loss={np.mean(losses[-100:]):.4f} elapsed={(time.time()-start)/60:.1f}min", log_fp)

    torch.save({"step": step, "state_dict": model.state_dict(), "losses": losses},
               str(out_dir / "stage1_final.pt"))
    log(f"Stage 1 done: {step} steps", log_fp)
    return losses


def stage2_structure(model, device, dataset: StanfordRNADataset, out_dir: Path, log_fp,
                     deadline: float, batch_size: int = 8):
    """Structure supervision: coord + confidence joint loss on sequences with ground truth."""
    log("=== Stage 2: structure supervision ===", log_fp)
    # Filter to samples WITH ground-truth available
    log("Filtering dataset for samples with valid ground truth...", log_fp)
    valid_indices = []
    for i in range(len(dataset)):
        sample = dataset[i]
        if sample["valid"].any():
            valid_indices.append(i)
        if i % 100 == 0:
            log(f"  scanned {i+1}/{len(dataset)}, {len(valid_indices)} with ground truth so far", log_fp)
    log(f"Final: {len(valid_indices)} samples with valid coords (of {len(dataset)})", log_fp)

    class _SubDS(torch.utils.data.Dataset):
        def __len__(self): return len(valid_indices)
        def __getitem__(self, i): return dataset[valid_indices[i]]

    def _collate(batch):
        return {
            "token_ids": torch.stack([b["token_ids"] for b in batch]),
            "coords":    torch.stack([b["coords"]    for b in batch]),
            "valid":     torch.stack([b["valid"]     for b in batch]),
            "target_id": [b["target_id"] for b in batch],
            "split":     [b["split"]     for b in batch],
            "category":  [b["category"]  for b in batch],
        }

    # Separate train (split=="train") from val/test for supervision
    train_indices = [i for i in valid_indices if dataset.df.iloc[i]["split"] == "train"]
    log(f"Training pool (split=train): {len(train_indices)}", log_fp)
    class _TrainDS(torch.utils.data.Dataset):
        def __len__(self): return len(train_indices)
        def __getitem__(self, i): return dataset[train_indices[i]]

    loader = DataLoader(_TrainDS(), batch_size=batch_size, shuffle=True, num_workers=2,
                        pin_memory=True, drop_last=True, collate_fn=_collate)

    optim = torch.optim.AdamW(model.parameters(), lr=1e-4, weight_decay=0.01)
    scaler = torch.amp.GradScaler("cuda") if device.type == "cuda" else None
    model.train()
    step, losses = 0, []
    start = time.time()
    last_ckpt = start
    last_commit = start

    it = iter(loader)
    while time.time() < deadline:
        try: batch = next(it)
        except StopIteration:
            it = iter(loader); batch = next(it)
        ids = batch["token_ids"].to(device, non_blocking=True)
        coords = batch["coords"].to(device, non_blocking=True)
        valid = batch["valid"].to(device, non_blocking=True)

        optim.zero_grad(set_to_none=True)
        if scaler is not None:
            with torch.amp.autocast("cuda", dtype=torch.bfloat16):
                out = model(ids)
                struct = structure_loss(out["coords"].float(), coords, valid)
                conf = confidence_loss(out["confidence"].float(), out["coords"].float(), coords, valid)
                loss = struct["total"] + 0.5 * conf
            scaler.scale(loss).backward()
            scaler.unscale_(optim)
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(optim); scaler.update()
        else:
            out = model(ids)
            struct = structure_loss(out["coords"], coords, valid)
            conf = confidence_loss(out["confidence"], out["coords"], coords, valid)
            loss = struct["total"] + 0.5 * conf
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optim.step()
        losses.append({"total": float(loss.item()), "coord": float(struct["coord_mse"].item()), "pair": float(struct["pair_mse"].item()), "conf": float(conf.item())})
        step += 1
        if step % 25 == 0:
            recent = pd.DataFrame(losses[-25:]).mean()
            log(f"stage2 step={step} total={recent['total']:.3f} coord={recent['coord']:.2f} pair={recent['pair']:.2f} conf={recent['conf']:.3f} elapsed={(time.time()-start)/60:.1f}min", log_fp)

        if time.time() - last_ckpt > 600:  # 10 min
            torch.save({"step": step, "state_dict": model.state_dict()}, str(out_dir / f"stage2_step_{step:06d}.pt"))
            pd.DataFrame(losses).to_csv(out_dir / "stage2_losses.csv", index=False)
            last_ckpt = time.time()
            log(f"checkpoint saved (step {step})", log_fp)

        if time.time() - last_commit > 900:  # 15 min commit cadence
            try:
                # Don't push the .pt (too big); push the CSV log + summary
                git_commit_push(Path("/workspace/rna3d"),
                                f"chore(GPU_training): incremental log @ step {step}",
                                [".research/30_experiments/runs/GPU_training_stage2/"])
            except Exception as e:
                log(f"git push warn: {e}", log_fp)
            last_commit = time.time()

    torch.save({"step": step, "state_dict": model.state_dict(), "losses": losses}, str(out_dir / "stage2_final.pt"))
    pd.DataFrame(losses).to_csv(out_dir / "stage2_losses.csv", index=False)
    log(f"Stage 2 done: {step} steps in {(time.time()-start)/60:.1f}min", log_fp)
    return losses


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hours", type=float, default=2.8, help="total training budget in hours")
    ap.add_argument("--data-root", type=str, default="/workspace/rna3d/data")
    ap.add_argument("--out-dir", type=str, default="/workspace/rna3d/.research/30_experiments/runs/GPU_training_stage2")
    ap.add_argument("--repo-root", type=str, default="/workspace/rna3d")
    ap.add_argument("--max-len", type=int, default=256)
    ap.add_argument("--batch-mlm", type=int, default=64)
    ap.add_argument("--batch-struct", type=int, default=8)
    ap.add_argument("--cpu", action="store_true", help="force CPU (for Mac smoke test)")
    ap.add_argument("--synthetic", action="store_true", help="use synthetic data (smoke test)")
    ap.add_argument("--smoke", action="store_true", help="quick 30s smoke run (no commit)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_fp = open(out_dir / "train.log", "a")
    log(f"=== START args={vars(args)} ===", log_fp)

    device = torch.device("cuda" if torch.cuda.is_available() and not args.cpu else "cpu")
    log(f"Device: {device}", log_fp)
    if device.type == "cuda":
        log(f"GPU: {torch.cuda.get_device_name(0)}", log_fp)

    cfg = ModelConfig(max_len=args.max_len)
    model = RNAStructureModel(cfg).to(device)
    log(f"Model: {model.num_params/1e6:.2f}M params", log_fp)

    start = time.time()
    if args.smoke:
        deadline = start + 30  # 30 second smoke
    else:
        deadline = start + args.hours * 3600

    if args.synthetic:
        # Single combined stage with synthetic data for smoke test
        ds = SyntheticRNADataset(n=32, max_len=cfg.max_len)
        # Run a few structure-loss steps
        loader = DataLoader(ds, batch_size=4, shuffle=True)
        optim = torch.optim.AdamW(model.parameters(), lr=1e-4)
        model.train()
        steps = 0
        for batch in loader:
            ids = batch["token_ids"].to(device)
            coords = batch["coords"].to(device)
            valid = batch["valid"].to(device)
            optim.zero_grad()
            out = model(ids)
            struct = structure_loss(out["coords"], coords, valid)
            conf = confidence_loss(out["confidence"], out["coords"], coords, valid)
            loss = struct["total"] + 0.5 * conf
            loss.backward()
            optim.step()
            steps += 1
            log(f"smoke step={steps} loss={loss.item():.3f}", log_fp)
            if time.time() > deadline: break
        log("=== smoke DONE ===", log_fp)
        return

    # Real pod run
    stage1_pretrain(model, device, Path(args.data_root) / "kaggle_raw" / "train_sequences.v2.csv",
                    out_dir, log_fp, deadline=deadline, max_len=args.max_len, batch_size=args.batch_mlm)

    if time.time() >= deadline:
        log("Deadline reached after stage 1; skipping stage 2", log_fp)
    else:
        category_csv = Path(args.repo_root) / ".research/30_experiments/runs/01b_description_categorization_v2/results/category_assignment.csv"
        ds = StanfordRNADataset(data_root=Path(args.data_root), splits=("train", "val", "test"),
                                max_len=args.max_len, category_csv=category_csv)
        stage2_structure(model, device, ds, out_dir, log_fp, deadline, batch_size=args.batch_struct)

    # Final commit
    git_commit_push(Path(args.repo_root),
                    f"feat(GPU_training): final — completed in {(time.time()-start)/60:.0f} min",
                    [args.out_dir])
    log("=== DONE ===", log_fp)
    log_fp.close()


if __name__ == "__main__":
    main()
