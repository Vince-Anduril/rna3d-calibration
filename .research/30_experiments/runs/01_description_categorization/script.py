"""01_description_categorization — apply 11-class biological scheme to Stanford CSVs.

Reads train/validation/test sequence CSVs, applies 11 priority-ordered regex
rules (case-insensitive) on the `description` column, with a length gate on
`nmr_solution_motif` (requires len(sequence) < 60). Writes:
  - results/category_assignment.csv  (all rows)
  - results/category_distribution.png (grouped bar chart)
  - results/uncategorized.csv         (rows with category == 'other')

Confidence:
  - "high"   : exactly one of the 11 specific rules matched (other not counted)
  - "medium" : >=2 specific rules matched (priority resolved the conflict)
  - "low"    : no specific rule matched -> fell back to 'other'

Re-run:
  source /workspace/rna3d/venv/bin/activate
  cd /workspace/rna3d/.research/30_experiments/runs/01_description_categorization
  python script.py
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

# Priority-ordered rules (first match wins). 'other' is the fallback.
RULES: list[tuple[str, str]] = [
    ("ribosome_subunit",
     r"\bribosom(e|al|os)\b"
     r"|\b(small|large)\s+subunit\b"
     r"|\b(30|40|50|60|70|80)s\b"
     r"|\b(5|5\.8|16|18|23|28)s\s*r?rna\b"
     r"|\brrna\b"),
    ("riboswitch",
     r"\briboswitch(es)?\b"
     r"|\baptamer\b"
     r"|\b(sam|tpp|fmn|cobalamin|glycine|lysine|preq[12]?|guanine|adenine)\b.*\b(riboswitch|aptamer)\b"
     r"|\b(riboswitch|aptamer)\b.*\b(sam|tpp|fmn|cobalamin|glycine|lysine|preq[12]?|guanine|adenine)\b"),
    ("ribozyme",
     r"\bribozyme(s)?\b"
     r"|\bhammerhead\b"
     r"|\bhairpin\s+ribozyme\b"
     r"|\bhdv\b"
     r"|\bgroup\s+i+\b"
     r"|\brnase\s*p\b"
     r"|\bspliceosom"),
    ("tRNA",
     r"\btrna\b|\btransfer\s+rna\b"),
    ("viral_rna",
     r"\bhiv\b|\bhcv\b|\bsars\b|\bires\b|\bframeshift\w*\b|\btar\s+rna\b|\bviral\s+rna\b"),
    ("pseudoknot",
     r"\bpseudoknot(s)?\b"),
    ("nmr_solution_motif",
     r"\b(nmr|solution)\b"),  # ALSO requires len(sequence) < 60
    ("loop_motif",
     r"\b(loop|hairpin|kink[-\s]?turn|sarcin[-\s]?ricin|tetraloop)\b"),
    ("complex_with_protein",
     r"\bcomplex\b.*\b(protein|with)\b"
     r"|\b(bound\s+to|in\s+complex\s+with)\b.*\bprotein\b"
     r"|\brna[-\s]protein\s+complex\b"),
    ("synthetic_designed",
     r"\b(designed|engineered|synthetic|chimer\w*)\b"),
]

# Pre-compile
COMPILED = [(name, re.compile(pat, re.IGNORECASE | re.DOTALL)) for name, pat in RULES]


def categorize(description: str, seq_length: int) -> tuple[str, str, str]:
    """Return (category, matched_rule_regex_snippet, confidence)."""
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

    chosen_name, chosen_pat = matched[0]  # priority = list order
    snippet = chosen_pat.pattern
    if len(snippet) > 80:
        snippet = snippet[:77] + "..."
    confidence = "high" if len(matched) == 1 else "medium"
    return (chosen_name, snippet, confidence)


def main() -> None:
    rows: list[dict] = []
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
    out_csv = RESULTS / "category_assignment.csv"
    assign.to_csv(out_csv, index=False)
    print(f"[write] {out_csv}  ({len(assign)} rows)")

    uncat = assign[assign["category"] == "other"]
    out_uncat = RESULTS / "uncategorized.csv"
    uncat.to_csv(out_uncat, index=False)
    print(f"[write] {out_uncat}  ({len(uncat)} rows)")

    # Distribution table
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
    print("\n=== Category distribution ===")
    print(pivot.to_string())

    # Bar chart (grouped by split)
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
    ax.set_title("Description-based RNA categorization (Stanford Kaggle splits)")
    ax.axhline(VIABILITY_THRESHOLD, color="gray", linestyle="--", linewidth=0.8,
               label=f"viability threshold (>={VIABILITY_THRESHOLD} total)")
    ax.legend()
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    plt.tight_layout()
    out_png = RESULTS / "category_distribution.png"
    plt.savefig(out_png, dpi=130)
    print(f"[write] {out_png}")

    # SUMMARY.md
    lines: list[str] = []
    lines.append("# Run 01 — Description categorization — SUMMARY")
    lines.append("")
    lines.append(f"Scheme: 11 priority-ordered regex rules (v1, from Q-2026-05-18-002).")
    lines.append(f"Splits processed: train ({(assign.split=='train').sum()}), "
                 f"val ({(assign.split=='val').sum()}), test ({(assign.split=='test').sum()}).")
    lines.append(f"Total rows: {len(assign)}.")
    lines.append(f"Viability threshold: >= {VIABILITY_THRESHOLD} across train+val+test combined.")
    lines.append("")
    lines.append("## Per-category counts")
    lines.append("")
    lines.append("| category | train | val | test | total | viable |")
    lines.append("|---|---:|---:|---:|---:|:--:|")
    for cat in category_order:
        row = pivot.loc[cat]
        flag = "OK" if bool(row["viable"]) else "NO"
        lines.append(f"| {cat} | {int(row['train'])} | {int(row['val'])} | "
                     f"{int(row['test'])} | {int(row['total'])} | {flag} |")
    lines.append("")
    viable_cats = [c for c in category_order if bool(pivot.loc[c, "viable"]) and c != "other"]
    lines.append(f"**Viable categories (>= {VIABILITY_THRESHOLD} total, excluding `other`):** "
                 f"{', '.join(viable_cats) if viable_cats else '(none)'}.")
    lines.append("")
    # Confidence stats
    conf = assign.groupby("confidence").size().to_dict()
    lines.append("## Confidence distribution")
    lines.append("")
    for c in ["high", "medium", "low"]:
        lines.append(f"- {c}: {conf.get(c, 0)}")
    lines.append("")
    # Ambiguity: medium confidence rows = first-match resolved conflict
    medium = assign[assign["confidence"] == "medium"]
    lines.append(f"## Ambiguity / edge cases")
    lines.append("")
    lines.append(f"- {len(medium)} rows had >=2 rules match before priority resolution "
                 f"(medium confidence). Inspect by filtering `category_assignment.csv` "
                 f"on `confidence == 'medium'`.")
    lines.append(f"- {len(uncat)} rows fell back to `other` (low confidence). "
                 f"See `uncategorized.csv`.")
    if len(uncat):
        sample = uncat["description_truncated_80chars"].head(5).tolist()
        lines.append("- Sample uncategorized descriptions:")
        for s in sample:
            lines.append(f"    - {s!r}")
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- `nmr_solution_motif` requires both (`nmr` or `solution`) AND "
                 "`len(sequence) < 60`; without the length gate that category absorbs "
                 "most cryo-EM ribosome assemblies that mention `solution`.")
    lines.append("- Priority order (ribosome -> riboswitch -> ribozyme -> tRNA -> "
                 "viral -> pseudoknot -> nmr_solution_motif -> loop -> "
                 "complex_with_protein -> synthetic_designed -> other) is taken "
                 "verbatim from the scientist's Q-002 answer.")
    lines.append("- All regexes case-insensitive. DOTALL is enabled so multi-line "
                 "descriptions (validation/test rows wrap lines inside quotes) are "
                 "matched correctly.")
    lines.append("")
    summary_path = RESULTS / "SUMMARY.md"
    summary_path.write_text("\n".join(lines))
    print(f"[write] {summary_path}")
    print("\n=== Done ===")


if __name__ == "__main__":
    main()
