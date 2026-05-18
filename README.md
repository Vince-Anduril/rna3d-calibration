# RNA Structure Prediction — Focused Pattern Analysis

Research project targeting a first published paper from the Stanford RNA 3D Folding
Kaggle benchmark.

**Thesis direction (2026-05-18):** focused deep-dive on ONE specific, under-studied
RNA family/molecule, identified opportunistically from data signals in Sprint 0,
where we find a novel structural pattern. Calibration analysis is a tool inside
this focus, not the headline.

**Framing:** humble and additive. *NOT* "Stanford / SOTA is wrong." *Instead*:
"Here is a focused look at family X — a class that has received less attention —
and a pattern we believe is worth the community's notice."

**Status:** Sprint 0 (exploratory, ~2 weeks).

## Layout

- `.research/` — shared knowledge base (literature, hypotheses, experiments, findings, inter-agent mailbox)
- `code/` — Python library
- `notebooks/` — exploratory Jupyter notebooks
- `scripts/` — orchestration helpers (e.g., regenerate_inbox.sh)
- `docs/superpowers/specs/` — design specs
- `docs/superpowers/plans/` — implementation plans
- `data/` — Kaggle dataset and external annotations (mostly gitignored)
- `models/` — model weights (gitignored)
- `paper/` — manuscript drafts

## Execution model

- **Pod (this machine)** — canonical for everything: git repo, source, data, weights,
  experiments, knowledge base.
- **Mac (orchestrator)** — runs the Claude Code session, hosts the three subagent
  definitions in `~/.claude/agents/`, and mirrors important outputs (figures,
  FINDINGS.md, paper drafts). No computation on the Mac. The Kaggle dataset never
  touches the Mac.

## Hardware

- GPU: RTX 5090, 32 GB VRAM (since 2026-05-18; previous pod with RTX 4090 24 GB
  was lost when /workspace failed to persist across availability migration).

## Subagent team

Three Claude Code subagents collaborate via the file-based knowledge base:

- **rna-research** — research scout: literature review + cross-domain opportunity hunting.
- **rna-scientist** — owns falsifiable hypotheses, protocols, and findings interpretation.
- **rna-code** — implements protocols on the pod via SSH; reports raw results.

The orchestrator (main Claude Code session) routes work via `99_protocols/INBOX.md`
and `questions_for_<agent>.md` mailboxes.
