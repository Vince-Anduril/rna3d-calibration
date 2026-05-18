"""01b_description_categorization_v2 — v2 regex refinements.

v2 changes vs v1 (Q-005 follow-ups absorbed via Q-007):
  - ribozyme rule additions: `\bRibonuclease P RNA\b` (case-insensitive).
  - viral_rna rule additions: `\bpicornaviral\b|\bbacteriophage\b|\bcoliphage\b`.

Mirrors the v1 layout exactly. Writes a v2-vs-v1 delta table.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).parent
RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)

DATA_DIR = Path("/workspace/rna3d/data/kaggle_raw")
SPLITS = {
    "train": DATA_DIR / "train_sequences.csv",
    "val": DATA_DIR / "validation_sequences.csv",
    "test": DATA_DIR / "test_sequences.csv",
}
NMR_LENGTH_THRESHOLD = 60
VIABILITY_THRESHOLD = 10

RULES: list[tuple[str, str]] = [
    ("ribosome_subunit",
     r"\bribosom(e|al|os)\b"
     r"|\b(small|large)\s+subunit\b"
     r"|\b(30|40|50|60|70|80)s\b"
     r"|\b(5|5\.8|16|18|23|28)s\s*r?rna\b"
     r"|\brrna\b"),
    ("riboswitch",
     r"\briboswitch(es)?\b|\baptamer\b"
     r"|\b(sam|tpp|fmn|cobalamin|glycine|lysine|preq[12]?|guanine|adenine)\b.*\b(riboswitch|aptamer)\b"
     r"|\b(riboswitch|aptamer)\b.*\b(sam|tpp|fmn|cobalamin|glycine|lysine|preq[12]?|guanine|adenine)\b"),
    ("ribozyme",
     r"\bribozyme(s)?\b|\bhammerhead\b|\bhairpin\s+ribozyme\b|\bhdv\b"
     r"|\bgroup\s+i+\b|\brnase\s*p\b|\bspliceosom"
     r"|\bribonuclease\s+p\s+rna\b|\bribonuclease\s+p\b"),  # v2 rescue
    ("tRNA", r"\btrna\b|\btransfer\s+rna\b"),
    ("viral_rna",
     r"\bhiv\b|\bhcv\b|\bsars\b|\bires\b|\bframeshift\w*\b|\btar\s+rna\b|\bviral\s+rna\b"
     r"|\bpicornaviral\b|\bbacteriophage\b|\bcoliphage\b"),  # v2 rescue
    ("pseudoknot", r"\bpseudoknot(s)?\b"),
    ("nmr_solution_motif", r"\b(nmr|solution)\b"),
    ("loop_motif", r"\b(loop|hairpin|kink[-\s]?turn|sarcin[-\s]?ricin|tetraloop)\b"),
    ("complex_with_protein",
     r"\bcomplex\b.*\b(protein|with)\b"
     r"|\b(bound\s+to|in\s+complex\s+with)\b.*\bprotein\b"
     r"|\brna[-\s]protein\s+complex\b"),
    ("synthetic_designed", r"\b(designed|engineered|synthetic|chimer\w*)\b"),
]

COMPILED = [(name, re.compile(pat, re.IGNORECASE | re.DOTALL)) for name, pat in RULES]


def categorize(description: str, seq_length: int) -> tuple[str, str, str]:
    if not isinstance(description, str):
        description = ""
    matched: list[tuple[str, re.Pattern[str]]] = []
    for name, pat in COMPILED:
        if name == "nmr_solution_motif":
            if pat.search(description) and seq_length < NMR_LENGTH_THRESHOLD:
                matched.append((name, pat))
        else:
            if pat.search(description):
                matched.append((name, pat))
    if not matched:
        return ("other", "", "low")
    chosen_name, chosen_pat = matched[0]
    snippet = chosen_pat.pattern
    if len(snippet) > 80:
        snippet = snippet[:77] + "..."
    confidence = "high" if len(matched) == 1 else "medium"
    return (chosen_name, snippet, confidence)


def main() -> None:
    rows = []
    for split_name, path in SPLITS.items():
        df = pd.read_csv(path)
        print(f"[load] {split_name}: {len(df)} rows from {path}")
        for _, r in df.iterrows():
            seq = r.get("sequence", "") or ""
            desc = r.get("description", "") or ""
            length = len(seq)
            category, snippet, confidence = categorize(desc, length)
            desc_trunc = (desc[:80] + "...") if isinstance(desc, str) and len(desc) > 80 else desc
            rows.append({
                "split": split_name,
                "target_id": r.get("target_id", ""),
                "sequence_length": length,
                "description": desc,
                "description_truncated_80chars": desc_trunc,
                "category": category,
                "matched_rule": snippet,
                "confidence": confidence,
            })

    assign = pd.DataFrame(rows)
    assign.to_csv(RESULTS / "category_assignment.csv", index=False)
    uncat = assign[assign["category"] == "other"]
    uncat.to_csv(RESULTS / "uncategorized.csv", index=False)

    category_order = [name for name, _ in RULES] + ["other"]
    pivot = (
        assign.groupby(["category", "split"])
        .size()
        .unstack(fill_value=0)
        .reindex(category_order, fill_value=0)
    )
    for s in ["train", "val", "test"]:
        if s not in pivot.columns:
            pivot[s] = 0
    pivot = pivot[["train", "val", "test"]]
    pivot["total"] = pivot.sum(axis=1)
    pivot["viable"] = pivot["total"] >= VIABILITY_THRESHOLD

    # v1 totals (from prior run) for delta-table
    V1_TOTALS = {
        "ribosome_subunit": 199,
        "riboswitch": 21,
        "ribozyme": 58,
        "tRNA": 66,
        "viral_rna": 77,
        "pseudoknot": 11,
        "nmr_solution_motif": 131,
        "loop_motif": 52,
        "complex_with_protein": 43,
        "synthetic_designed": 4,
        "other": 206,
    }

    fig, ax = plt.subplots(figsize=(12, 6))
    x = range(len(category_order))
    width = 0.27
    ax.bar([i - width for i in x], pivot["train"], width=width, label="train", color="#3b6db5")
    ax.bar([i for i in x], pivot["val"], width=width, label="val", color="#d49a31")
    ax.bar([i + width for i in x], pivot["test"], width=width, label="test", color="#5e9c5c")
    ax.set_xticks(list(x))
    ax.set_xticklabels(category_order, rotation=35, ha="right")
    ax.set_ylabel("Number of sequences")
    ax.set_xlabel("Category (priority order, first match wins)")
    ax.set_title("Description-based RNA categorization v2 (Stanford Kaggle splits)")
    ax.axhline(VIABILITY_THRESHOLD, color="gray", linestyle="--", linewidth=0.8,
               label=f"viability threshold (>={VIABILITY_THRESHOLD} total)")
    ax.legend()
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    plt.tight_layout()
    plt.savefig(RESULTS / "category_distribution.png", dpi=130)

    lines = []
    lines.append("# Run 01b — Description categorization v2 — SUMMARY")
    lines.append("")
    lines.append("Scheme v2: same as v1 + Q-005 rescues:")
    lines.append("- `ribozyme` rescues `Ribonuclease P RNA` (case-insensitive).")
    lines.append("- `viral_rna` rescues `picornaviral|bacteriophage|coliphage`.")
    lines.append("")
    lines.append(f"Splits: train ({(assign.split=='train').sum()}), val ({(assign.split=='val').sum()}), test ({(assign.split=='test').sum()}). Total: {len(assign)}.")
    lines.append("")
    lines.append("## Per-category counts (v2)")
    lines.append("")
    lines.append("| category | train | val | test | total | v1 total | delta | viable |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|:--:|")
    for cat in category_order:
        row = pivot.loc[cat]
        v1 = V1_TOTALS.get(cat, 0)
        delta = int(row["total"]) - v1
        delta_str = f"+{delta}" if delta > 0 else (f"{delta}" if delta < 0 else "0")
        flag = "OK" if bool(row["viable"]) else "NO"
        lines.append(f"| {cat} | {int(row['train'])} | {int(row['val'])} | {int(row['test'])} | {int(row['total'])} | {v1} | {delta_str} | {flag} |")
    lines.append("")
    viable_cats = [c for c in category_order if bool(pivot.loc[c, "viable"]) and c != "other"]
    lines.append(f"**Viable categories (v2):** {', '.join(viable_cats) if viable_cats else '(none)'}.")
    lines.append("")
    conf = assign.groupby("confidence").size().to_dict()
    lines.append("## Confidence distribution")
    lines.append("")
    for c in ["high", "medium", "low"]:
        lines.append(f"- {c}: {conf.get(c, 0)}")
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- Delta vs v1: rescues should *reduce* `other` and *increase* `ribozyme` + `viral_rna`. See delta column.")
    lines.append("- Otherwise identical to v1 methodology; no other rules changed.")
    (RESULTS / "SUMMARY.md").write_text("\n".join(lines))
    print("[done] v2 categorization complete.")


if __name__ == "__main__":
    main()
