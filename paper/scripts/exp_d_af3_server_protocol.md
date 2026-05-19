# Extension D — AlphaFold Server (AF3) submission protocol

> Manual user action required. AF3 is not redistributable; the AlphaFold
> Server is the official Google interface (free Google account, ~20 jobs/day).

## Sequences to submit (7 jobs total)

Copy each FASTA below into a separate AlphaFold Server job. Use the
**RNA** entity type and set the sequence as a single chain.

### 1. R1107_human (CPEB3 ribozyme — Homo sapiens)
```
>R1107_human
GGGGGCCACAGCAGAAGCGUUCACGUCGCAGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU
```

### 2. R1108_chimp (CPEB3 ribozyme — Pan troglodytes)
```
>R1108_chimp
GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU
```

### 3. null_pos5 (G5A peripheral control)
```
>null_pos5
GGGGACCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU
```

### 4. null_pos20 (U20A peripheral control)
```
>null_pos20
GGGGGCCACAGCAGAAGCGAUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU
```

### 5. null_pos41 (A41C peripheral control)
```
>null_pos41
GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCCGCCAUUGCACUCCGGCUGCGAAUUCUGCU
```

### 6. null_pos55 (U55A peripheral control)
```
>null_pos55
GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCAGCGAAUUCUGCU
```

### 7. null_pos64 (U64A peripheral control)
```
>null_pos64
GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUACUGCU
```

## Output to capture

For each job, AF3 Server gives you a downloadable ZIP containing:

- `model_0.cif` (the predicted structure)
- `confidences.json` (per-residue pLDDT + ipTM, etc.)
- `summary_confidences.json` (aggregate confidence)
- `terms_of_use.md`

## Where to put the files on the pod

```bash
# After downloading each job's ZIP, on the Mac:
mkdir -p /Users/leduigouvincent/rna3d-mirror/data/af3/{R1107_human,R1108_chimp,null_pos5,null_pos20,null_pos41,null_pos55,null_pos64}
# Unzip each into the matching subdirectory.
# Then SCP up to pod:
scp -r /Users/leduigouvincent/rna3d-mirror/data/af3/ root@<POD_IP>:-p <POD_PORT>:/workspace/rna3d/data/af3/
```

## Analysis script

`paper/scripts/exp_d_af3_compare.py` will:

1. Load `model_0.cif` for each target via biotite.
2. Compute per-residue Kabsch deltas R1107 vs R1108.
3. Extract top-5 divergent residues.
4. Compare to DRfold2 top-5: same anchors? Different?
5. Also load `confidences.json` for per-residue pLDDT and report at cascade residues.
6. Write `af3_results.md` with a side-by-side comparison.

## What would the AF3 result mean?

| AF3 outcome | Story |
|---|---|
| AF3 puts **9 AND 60** in top-5 (like DRfold2) | Strong cross-model anchor — F-001 robust |
| AF3 top-5 = totally different residues | DRfold2 has model-specific bias; can't generalize the claim |
| AF3 puts NEITHER 9 nor 60 in top-5 | **Publishability multiplier:** "the strongest open SOTA does not respond to a 4×-activity-changing single-nt variant" |
| AF3 high pLDDT at residue 30 in human (where structure is "wrong" per Skilandat) | Calibration failure — high-stakes claim |

Any of these is publishable. The third (AF3 misses the variant entirely) would be
particularly attention-grabbing.
