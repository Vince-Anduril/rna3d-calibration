# Project Context Snapshot

> Maintained at the project root so any new session — Claude or human — can pick
> up the state in one read. Update at every major decision or pivot.

## Last update: 2026-05-18

## Where we are

- Project bootstrap is in the **re-execution** phase after pod loss. The previous
  pod (RTX 4090) became unavailable and `/workspace` was lost. All work is being
  re-created on the current pod (RTX 5090, 32 GB).
- Sprint 0 of the research project (see Section 4.2 of the original spec) is
  starting. No real inference has run yet.
- Three Claude Code subagents are *designed* but not yet *implemented* in
  `~/.claude/agents/` on the Mac.

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

1. Re-extract Kaggle dataset (download in progress).
2. Set up Python venv with pinned deps.
3. Bootstrap `.research/` knowledge base with the **revised thesis** in
   `00_thesis.md`.
4. Write the three subagent definitions in `~/.claude/agents/` on the Mac.
5. Smoke-test each subagent (boundary enforcement, mailbox protocol).
6. Set up private GitHub mirror (pod pushes).
7. Start Sprint 0 Week 1: rna-research surfaces 3-5 candidate under-studied
   families; rna-scientist+rna-code measure data-signals on each; user picks the
   focus.

## Memory of decisions made (chronological highlights)

- 2026-05-17 — Brainstorm: angle C (error analysis) chosen, later refined to
  calibration. Three-agent architecture validated.
- 2026-05-17 — Pod-only execution rule established.
- 2026-05-18 — Old pod lost. New pod RTX 5090.
- 2026-05-18 — Thesis pivoted from "calibration analysis" to "focused finding on
  under-studied family, calibration as tool not headline". User explicit
  framing: "find a focused point/pattern, do not position as Stanford being
  wrong."
