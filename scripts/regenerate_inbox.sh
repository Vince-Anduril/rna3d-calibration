#!/usr/bin/env bash
# Rebuild .research/99_protocols/INBOX.md as the index of all [OPEN] entries
# across the three mailbox files. Run before deciding which agent to invoke.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INBOX="${ROOT_DIR}/.research/99_protocols/INBOX.md"
MAILBOXES=(
  "${ROOT_DIR}/.research/99_protocols/questions_for_research.md"
  "${ROOT_DIR}/.research/99_protocols/questions_for_scientist.md"
  "${ROOT_DIR}/.research/99_protocols/questions_for_code.md"
)

{
  echo "# INBOX — Open inter-agent questions"
  echo ""
  echo "_Regenerated $(date -u +'%Y-%m-%dT%H:%M:%SZ') by scripts/regenerate_inbox.sh._"
  echo ""
  for mb in "${MAILBOXES[@]}"; do
    [ -f "$mb" ] || continue
    name="$(basename "$mb" .md)"
    matches=$(grep -nE '^##\s+\[OPEN\]' "$mb" || true)
    if [ -n "$matches" ]; then
      echo "## $name"
      echo ""
      while IFS= read -r line; do
        lineno="${line%%:*}"
        rest="${line#*:}"
        echo "- $rest  (\`$name.md:$lineno\`)"
      done <<< "$matches"
      echo ""
    fi
  done
} > "$INBOX"

echo "Wrote $INBOX"
wc -l "$INBOX"
