# CPEB3 Probe — Research Roadmap (post-n=20 → 4 extensions)

> Status as of 2026-05-19. n=20 NULL controls completes overnight (in flight).
> This document lists the 4 next angles that would graduate the work from
> "sanity check" to genuine discovery, with cost/risk/payoff estimates.

## Where we are

| Done | What |
|---|---|
| ✅ | DRfold2 inference on R1107 + R1108 + 6 NULL controls |
| ✅ | Top-5 cascade analysis (residues 9 and 60 specific to real mutation) |
| ✅ | Exp #1: 3D distance check (cascade at 9/60 = direct contact, 22/24/51 = DRfold2 attractor) |
| ✅ | Exp #2 (n=20 NULL controls, p-value) — running on pod overnight |
| ✅ | Honest paper draft v0.4 + Canva-style slides (FR + EN) |

## The 4 extensions

Ordered by **leverage / cost** (best first):

### 🥇 Extension A — Confidence channel (DRfold2 .ret decode)
**Effort:** Low (pod-side decoding + Python analysis, ~30 min total)
**Payoff:** Could be high — turns the methodology note into a genuine calibration claim.

DRfold2 saves `.ret` files per config-model (~17 per config × 4 configs = ~70 per
target). They are pickled Python dicts. We need to:
1. Decode the dict schema (find confidence/pLDDT/predicted-error keys).
2. For R1107 and R1108, extract per-residue confidence at the cascade
   residues (9, 22, 24, 51, 60).
3. Test: **Is DRfold2 confident at residues 9 and 60?**
   - HIGH confidence + LARGE delta = "confidently wrong" → miscalibration story
   - LOW confidence + LARGE delta = "the model itself knows it's uncertain" → calibration is correct
4. Per-position confidence heatmap R1107 vs R1108.

Script: `paper/scripts/exp_a_confidence_decode.py` (drafted).
Needs pod for one ~15 min CPU run, no GPU.

### 🥈 Extension B — Full 69-position mutational scan
**Effort:** Medium (~3.5h pod time, ~€1).
**Payoff:** High — replaces the n=5/n=20 NULL framing with a complete map.

For each of the 69 positions, generate one mutation and run DRfold2 inference,
then compute the (delta @ 9, delta @ 60) tuple. Output: 69×2 matrix that
visualizes as a heatmap. The question becomes: **which positions, if any,
activate the P1 anchor (both 9 AND 60 in top-5) when mutated?**

Expected positions: 30 (real) and possibly its direct P1.1 neighbors (29, 31)
plus the 5' P1 stem (1-9). If ONLY the position-30 mutation activates both
anchors out of 69, we have a clean specificity statement.

Script: `paper/scripts/exp_b_full_69_scan.sh` (drafted).

### 🥉 Extension C — Ancestral resurrection
**Effort:** Medium-high (data fetch + 5-10 DRfold2 runs, ~1h).
**Payoff:** Could elevate to "structural genomics" angle if signal holds.

Use Bendixsen *et al.* 2021 (MBE, [PMC8233481](https://pmc.ncbi.nlm.nih.gov/articles/PMC8233481/))
ancestral CPEB3 reconstructions. They published a consensus mammalian
ancestral sequence (67 nt) and node-by-node activity rates. Run DRfold2 on
each, compute structural divergence vs the ancestral, and test whether
predicted divergence tracks measured activity decline.

**Data fetch (manual step):**
1. Visit https://gitlab.com/devinbendixsen/cpeb3_phylo_fls
2. Download `genotypes.csv` (or equivalent) from their supp data S1
3. Place sequences at `/workspace/rna3d/data/ancestral/sequences.fasta`

Script: `paper/scripts/exp_c_ancestral.py` (drafted; needs sequences).

### 🏅 Extension D — AlphaFold3 server bench
**Effort:** Low compute, but manual browser action required.
**Payoff:** Highest publishability — AF3 is the gold standard. If AF3 fails
to discriminate R1107 from R1108, the title is "AF3 does not respond to a
4×-activity-changing single-nt variant", which is a real reviewer-grabbing claim.

**Manual steps for the user:**
1. Sign in at https://alphafoldserver.com/ (free Google account, ~20 jobs/day)
2. Submit each sequence as a single-RNA job (R1107, R1108, the 5 NULL pos)
3. Download the JSON output for each
4. Place files at `data/af3/<target_id>.json`

A parser script (`paper/scripts/exp_d_af3_compare.py`, drafted) then
ingests them and runs the same cascade analysis as for DRfold2 and RhoFold+.

## Combined budget if all 4 are done

| Item | Cost |
|---|---|
| A (confidence) | ~€0.05 |
| B (69-scan) | ~€1.20 |
| C (ancestral) | ~€0.30 |
| D (AF3 server) | €0 + user time |
| **Total** | **~€1.55 + ~30 min user submission** |

## Risk-adjusted recommendation

1. **Do A first** (low cost, could change the story alone if calibration is off).
2. **Do B in parallel** as the overnight workhorse — strongest paper claim if signal holds.
3. **Then C** — converts the work from a model-behavior note to evolutionary genomics if it works.
4. **D is the publishability multiplier** — even null AF3 result is interesting.
