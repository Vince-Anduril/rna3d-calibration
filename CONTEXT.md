# Project Context Snapshot

> Maintained at the project root so any new session — Claude or human — can pick
> up the state in one read. Update at every major decision or pivot.

## Last update: 2026-05-18 (afternoon — code authoring complete, awaiting pod)

## Where we are

- **Bootstrap complete.** Pod (RTX 5090, 32 GB) had:
  - Dataset extracted (~60 GB across 868 sequences + 8,672 PDB structures + 2,534 MSAs)
  - `.research/` knowledge base with 16 seed files
  - 3 Sprint 0 Week 1 cycles ran (categorization run 01, then 01b v2 regex)
  - The full agent team produced FOCUS_SHORTLIST.md (3 candidates: ribozyme STRONG, riboswitch STRONG, viral_rna MEDIUM)
- **Pod is currently OFF** (user stopped it because the first autonomous run did not saturate GPU; RunPod also reports availability issues for restart).
- **Authored Mac-side and pushed:** complete GPU training pipeline at `code/jobs_gpu/` (model + loss + data + calibration + train + eval + pod_launch.sh) with 9 unit tests passing on Mac CPU.
- **Three Claude Code subagents** are implemented at `~/.claude/agents/rna-{research,scientist,code}.md` on the Mac (Claude Code reads them from there).
- **Focus pick is FROZEN: ribozyme** (autonomous decision per user grant — see `docs/superpowers/specs/FOCUS_DECISION.md`).

## Thesis direction

NOT a general "calibration of SOTA models" attack-style paper. The framing has
been refined:

- **Pick ONE under-studied RNA family or motif** from the Stanford Kaggle test
  set, chosen **opportunistically from data signals** (e.g., largest inter-model
  error variance, motif most frequent but worst-predicted, structures where SOTA
  diverges most). This selection happens in Sprint 0 Week 1.
- **Deep-dive that family**: structural, biological, prediction-error properties.
- **Find a novel pattern** worth communicating to the field.
- Calibration analysis is a *tool* inside this focus, not the title.

Working title pattern: *"Structural prediction challenges in [chosen family X]:
a focused analysis of the Stanford RNA 3D benchmark."*

## Hard constraints

- Humble framing, never adversarial toward Stanford or SOTA.
- Pod-only execution. No computation on Mac.
- Data stays on pod. Mac mirrors only lightweight artifacts (figures, FINDINGS,
  paper).
- First paper for the author (intermediate ML, no prior publication). Bias toward
  rigor and clarity over methodological prouesse.

## Lost in the migration (must re-create on this pod)

- The original `docs/superpowers/specs/2026-05-17-rna-calibration-research-design.md`
  (354 lines, 3 review iterations). Substantively superseded by this CONTEXT.md
  for the active thesis; a formal revised spec will be re-written when time allows.
- The original `docs/superpowers/plans/2026-05-17-sprint-0-bootstrap.md` (1757 lines,
  reviewed and approved). Bootstrap is being re-executed pragmatically; a formal
  revised plan can be re-derived if needed.
- The empty `.research/` knowledge base seeds — will be re-created in the next
  bootstrap step.

## Active to-do (top of mind)

**Blocker:** pod is OFF, RunPod reports availability issues. Cannot launch heavy GPU work until user restarts the pod.

**When pod restarts (user action needed):**
1. User gives new IP/port (changes on every pod restart).
2. SCP nothing — pod just runs `git pull origin main` and gets the latest pipeline at `code/jobs_gpu/`.
3. Launch via `nohup bash code/jobs_gpu/pod_launch.sh > pod_run.log 2>&1 &` — runs ~3h autonomously.
4. Pipeline auto-commits and auto-pushes incrementally (every 15 min of training).

**While pod is off (autonomous Mac work continues):**
- Read the rna-research papers/ deep-dives if I get them in writable form.
- Refine the paper outline as findings start coming back.
- Tighten the eval script if anything in the Stage 2 logs surprises us.

**Files NOT to touch:**
- `.research/00_thesis.md` — read-only for non-scientist agents (decision log handled by scientist agent only)
- The Kaggle data on pod (gitignored anyway)
- Anything in `.research/30_experiments/runs/01*/` — those are scientist+code's owned territory and may have results from prior runs

## Memory of decisions made (chronological highlights)

- 2026-05-17 — Brainstorm: angle C (error analysis) chosen, later refined to
  calibration. Three-agent architecture validated.
- 2026-05-17 — Pod-only execution rule established.
- 2026-05-18 — Old pod lost. New pod RTX 5090.
- 2026-05-18 — Thesis pivoted from "calibration analysis" to "focused finding on
  under-studied family, calibration as tool not headline". User explicit
  framing: "find a focused point/pattern, do not position as Stanford being
  wrong."
