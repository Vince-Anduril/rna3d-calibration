# 30 — Experiment registry

> Owner: rna-code. One row per experimental run. Each run gets a directory
> under `30_experiments/runs/<run_id>/` containing code, config, results, figures.

| run_id | date | hypothesis | model(s) | status | notes |
|---|---|---|---|---|---|
| _none yet_ |  |  |  |  |  |
| 00_pod_smoke | 2026-05-18 | smoke | synthetic | done | Pod execution validated on RTX 5090 (torch 2.8 cu128); 844 train rows / 844 unique target_ids / mean seq len 162.43; top description words dominated by structure/rna/complex/solution/cryo/em/loop/ribosome/nmr. |
| 01_description_categorization | 2026-05-18 | H-001 (focus selection) | regex/pandas | done | 11-class scheme (v1 from Q-002) applied to train+val+test (868 rows). 10/11 categories pass >=10 viability; only `synthetic_designed` fails (n=4). 206 rows -> `other` (mostly RNA-protein complexes lacking the keyword `complex`+`protein` co-occurrence and bacteriophage/picornaviral entries). See runs/01_description_categorization/results/SUMMARY.md. |
