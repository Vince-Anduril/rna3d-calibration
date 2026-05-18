"""Pod smoke test: validate pod execution path + extract a tiny project-meaningful signal.

Prints torch+cuda info, loads train_sequences.csv, reports basic stats, and prints the
15 most common single words (lowercase, punctuation-stripped) in the description column
as a quick proxy for represented biological classes.
"""
import re
from collections import Counter

import pandas as pd
import torch


def main() -> None:
    print("=== Torch / CUDA ===")
    print(f"torch.__version__       = {torch.__version__}")
    print(f"torch.cuda.is_available = {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"torch.cuda.device_count = {torch.cuda.device_count()}")
        print(f"torch.cuda.get_device_name(0) = {torch.cuda.get_device_name(0)}")

    csv_path = "/workspace/rna3d/data/kaggle_raw/train_sequences.csv"
    print(f"\n=== Loading {csv_path} ===")
    df = pd.read_csv(csv_path)
    print(f"columns                 = {list(df.columns)}")
    print(f"row count               = {len(df)}")

    if "target_id" in df.columns:
        print(f"unique target_id count  = {df['target_id'].nunique()}")
    else:
        print("target_id column not present")

    if "sequence" in df.columns:
        seq_lens = df["sequence"].astype(str).str.len()
        print(f"average sequence length = {seq_lens.mean():.2f}")
        print(f"min / max seq length    = {seq_lens.min()} / {seq_lens.max()}")
    else:
        print("sequence column not present")

    print("\n=== 15 most common description words ===")
    if "description" in df.columns:
        word_re = re.compile(r"[a-z0-9]+")
        counter: Counter = Counter()
        for desc in df["description"].dropna().astype(str):
            counter.update(word_re.findall(desc.lower()))
        for word, count in counter.most_common(15):
            print(f"{count:>8d}  {word}")
    else:
        print("description column not present")


if __name__ == "__main__":
    main()
