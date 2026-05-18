# 30 — Experiment registry

> Owner: rna-code. One row per experimental run. Each run gets a directory
> under `30_experiments/runs/<run_id>/` containing code, config, results, figures.

| run_id | date | hypothesis | model(s) | status | notes |
|---|---|---|---|---|---|
| _none yet_ |  |  |  |  |  |
| 00_pod_smoke | 2026-05-18 | smoke | synthetic | done | Pod execution validated on RTX 5090 (torch 2.8 cu128); 844 train rows / 844 unique target_ids / mean seq len 162.43; top description words dominated by structure/rna/complex/solution/cryo/em/loop/ribosome/nmr. |
