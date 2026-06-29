#!/usr/bin/env bash
#
# doclint.sh -- lint the Bertini 2 C++ Doxygen documentation.
#
# This is the first foothold of a broader linting initiative.  It is deliberately
# cheap: it runs Doxygen only (no compilation, no wheel, no graphviz) and finishes
# in well under a minute, so it can gate every PR without the cost of a full build.
#
# It enforces two independent things, built from the canonical root Doxyfile with
# a few overrides piped in on stdin:
#
#   Pass 1 -- ERROR GATE (zero tolerance)
#     With EXTRACT_ALL=YES, Doxygen only warns about *wrong* documentation:
#     @param names that don't match the signature, doc blocks attached to a
#     signature that no longer exists, duplicate @param, unresolved \ref/\cite,
#     etc.  Missing docs are NOT flagged here, so the undocumented backlog never
#     blocks a build.  Any warning in this pass fails the lint.
#
#   Pass 2 -- UNDOCUMENTED RATCHET (monotonic)
#     With EXTRACT_ALL=NO, Doxygen flags every undocumented entity.  We count them
#     and compare against tools/doc_undocumented_baseline.txt.  The count may only
#     decrease: if it goes up, the lint fails; when it goes down, run with
#     --update-baseline to lock in the gain.  Target is 0, at which point the gate
#     can be made absolute.
#
# Usage:
#   tools/doclint.sh                  # run both passes (CI mode)
#   tools/doclint.sh --update-baseline  # rewrite the baseline from the current count
#
# Override the Doxygen binary with $DOXYGEN if it isn't on PATH.

set -euo pipefail

DOXYGEN="${DOXYGEN:-doxygen}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASELINE_FILE="$REPO_ROOT/tools/doc_undocumented_baseline.txt"
cd "$REPO_ROOT"

if ! command -v "$DOXYGEN" >/dev/null 2>&1; then
  echo "error: '$DOXYGEN' not found on PATH (set \$DOXYGEN to override)." >&2
  exit 127
fi

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

# Settings shared by both passes: generate no output at all (we only want the
# warnings from parsing -- this is the fastest mode and, unlike generating a
# format, Doxygen emits each warning exactly once).  HAVE_DOT off so graphviz
# isn't needed.
common_overrides() {
  cat <<EOF
OUTPUT_DIRECTORY = $WORKDIR/out
GENERATE_HTML    = NO
GENERATE_LATEX   = NO
GENERATE_XML     = NO
HAVE_DOT         = NO
QUIET            = YES
EOF
}

# A warning line worth counting: a real Doxygen "<file>:<line>: warning:" line,
# excluding two things that are not code/doc mismatches --
#   * the harmless "No output formats selected" notice (we deliberately emit none)
#   * \cite numbering, which depends on bibtex being installed (a doc-*build*
#     concern handled in build_docs.yml, not a signature/parameter mismatch)
real_warnings() {
  grep ': warning:' "$1" 2>/dev/null \
    | grep -v 'No output formats selected' \
    | grep -v '\\cite command' || true
}

run_doxygen() {
  # $1 = path for the warning log; remaining stdin lines (after the Doxyfile and
  # common overrides) are extra config appended by the caller.
  local logfile="$1"; shift
  { cat Doxyfile; common_overrides; cat; echo "WARN_LOGFILE = $logfile"; } \
    | "$DOXYGEN" - >/dev/null 2>&1 || true
}

###############################################################################
# Pass 1 -- error gate
###############################################################################
PASS1_LOG="$WORKDIR/pass1.log"
run_doxygen "$PASS1_LOG" <<'EOF'
EXTRACT_ALL           = YES
WARN_IF_DOC_ERROR     = YES
WARN_IF_INCOMPLETE_DOC = YES
WARN_NO_PARAMDOC      = NO
EOF

# Count real warning lines (each Doxygen warning begins "<file>:<line>: warning:").
PASS1_COUNT=$(real_warnings "$PASS1_LOG" | grep -c ': warning:' || true)

echo "== doclint pass 1: documentation correctness =="
if [[ "$PASS1_COUNT" -gt 0 ]]; then
  echo "FAIL: $PASS1_COUNT documentation error(s) -- docs describe code that does not exist:" >&2
  echo >&2
  real_warnings "$PASS1_LOG" >&2
  PASS1_OK=0
else
  echo "OK: no documentation correctness errors."
  PASS1_OK=1
fi

###############################################################################
# Pass 2 -- undocumented ratchet
###############################################################################
PASS2_LOG="$WORKDIR/pass2.log"
run_doxygen "$PASS2_LOG" <<'EOF'
EXTRACT_ALL           = NO
WARN_IF_UNDOCUMENTED  = YES
WARN_NO_PARAMDOC      = NO
EOF

# Entity-level undocumented count -- stable across runs for a given Doxygen version.
UNDOC_COUNT=$(grep -c 'is not documented' "$PASS2_LOG" || true)

if [[ "${1:-}" == "--update-baseline" ]]; then
  echo "$UNDOC_COUNT" > "$BASELINE_FILE"
  echo "== baseline updated: $UNDOC_COUNT undocumented entit(ies) =="
  # Still honor the pass-1 result so we never bless a baseline atop broken docs.
  [[ "$PASS1_OK" -eq 1 ]] || exit 1
  exit 0
fi

BASELINE=$(tr -d '[:space:]' < "$BASELINE_FILE" 2>/dev/null || echo "")
echo "== doclint pass 2: undocumented ratchet =="
if [[ -z "$BASELINE" ]]; then
  echo "note: no baseline file; current undocumented count is $UNDOC_COUNT."
  echo "      run 'tools/doclint.sh --update-baseline' to seed it."
  PASS2_OK=1
elif [[ "$UNDOC_COUNT" -gt "$BASELINE" ]]; then
  echo "FAIL: undocumented entities rose from $BASELINE to $UNDOC_COUNT (+$((UNDOC_COUNT - BASELINE)))." >&2
  echo "      document the new code, or this is genuinely new public surface." >&2
  PASS2_OK=0
else
  if [[ "$UNDOC_COUNT" -lt "$BASELINE" ]]; then
    echo "OK: undocumented entities dropped from $BASELINE to $UNDOC_COUNT (-$((BASELINE - UNDOC_COUNT)))."
    echo "    run 'tools/doclint.sh --update-baseline' to lock in the gain."
  else
    echo "OK: undocumented entities at baseline ($UNDOC_COUNT)."
  fi
  PASS2_OK=1
fi

[[ "$PASS1_OK" -eq 1 && "$PASS2_OK" -eq 1 ]]
