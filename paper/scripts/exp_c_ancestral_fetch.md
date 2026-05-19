# Extension C — Ancestral CPEB3 data fetch instructions

> Manual one-time setup. After this, `exp_c_ancestral.py` runs on the pod.

## Source

Bendixsen, Pollock, Peri & Hayden (2021). "Experimental Resurrection of
Ancestral Mammalian CPEB3 Ribozymes Reveals Deep Functional Conservation."
*Molecular Biology and Evolution* 38(7):2843-2861. PMID:33720319.

Open-access at [PMC8233481](https://pmc.ncbi.nlm.nih.gov/articles/PMC8233481/).

Their code + supplementary data is at
[gitlab.com/devinbendixsen/cpeb3_phylo_fls](https://gitlab.com/devinbendixsen/cpeb3_phylo_fls).

## What to fetch

You need a FASTA file with one record per reconstructed ancestor at the
phylogenetic tree nodes they reported, **labeled by node identity**, e.g.:

```
>ancestral_mammalian
GGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGC
>ancestral_eutherian
GGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGC
>ancestral_primate
...
```

Bendixsen 2021 reports the **majority of ancestral nodes** carry one
identical 67-nt sequence (same as 41 of the 97 extant mammals they
analyzed). You can probably get away with:

- 1 consensus ancestor (the dominant reconstructed mammalian sequence)
- 2-3 alternative reconstructions at deeper nodes (if reported separately)
- The 41 extant mammalian sequences sharing the consensus (for a robustness
  test)

## Step-by-step

```bash
# On your Mac (won't be needed on pod since /workspace is small):
mkdir -p /Users/leduigouvincent/rna3d-mirror/data/ancestral
cd /Users/leduigouvincent/rna3d-mirror/data/ancestral

# Option 1: clone the repo and look for the ancestral fasta or csv
git clone https://gitlab.com/devinbendixsen/cpeb3_phylo_fls.git
ls cpeb3_phylo_fls/  # explore the file tree

# Option 2: download supp file S1 from MBE article page
# Visit https://academic.oup.com/mbe/article/38/7/2843/6171150
# Click "Supplementary data" → download data_S1 (likely Excel or CSV with genotypes)

# Once you have the sequences extracted:
# build sequences.fasta with the node labels as record IDs.
# Length will be 67 nt; the pod-side script pads to 69 nt only if needed
# (it should not be needed — we keep ancestor length as-is).

# Build activity.csv from the same paper's supp data:
# ancestor_id,fraction_cleaved
# ancestral_mammalian,0.92
# ... (other ancestors and rates if reported)
```

## SCP to pod (when next pod is up)

```bash
POD_IP=...; POD_PORT=...
ssh -p $POD_PORT root@$POD_IP 'mkdir -p /workspace/rna3d/data/ancestral'
scp -P $POD_PORT data/ancestral/sequences.fasta root@$POD_IP:/workspace/rna3d/data/ancestral/
scp -P $POD_PORT data/ancestral/activity.csv root@$POD_IP:/workspace/rna3d/data/ancestral/
```

## Then run on pod

```bash
ssh -p $POD_PORT root@$POD_IP 'cd /workspace/rna3d && source venv/bin/activate && python paper/scripts/exp_c_ancestral.py'
```

## Expected output

- `.research/30_experiments/runs/cpeb3_focused/ancestral_results.csv`
- `.research/30_experiments/runs/cpeb3_focused/ancestral_results.md`

## What the result would mean

If structural divergence (RMSD vs human R1107) anti-correlates with measured
activity (more divergent = less active), then DRfold2's predicted structure
tracks experimentally measured function across millions of years of evolution
— a much stronger claim than the single-pair observation.
