#!/usr/bin/env bash
# render_runbook_direct.sh -- render a .qmd runbook on NEMO.
#
# Differs from the reproducible-runbook skill's render_runbook.sh in one respect,
# and it is not cosmetic: the skill documents `module load quarto/1.5.57-x64`,
# which FAILS on this host as of 2026-07-29 --
#   ERROR: Unable to locate a modulefile for 'quarto/1.5.57-x64'
# The binary is installed and reports 1.5.57; only the modulefile is missing. So
# quarto is invoked by absolute path here.
#
# Three things must line up, each failing differently:
#   1. quarto            -- by path, since the module does not resolve
#   2. QUARTO_PYTHON     -- must have bash_kernel installed (envs/runbook does)
#   3. pandoc            -- not on PATH; lives in envs/reanalysis_R/bin
#
# Run this as a SLURM job on ncpu, not on a login node: chunks that read BAMs or
# BED files hit the login-node memory cgroup, which kills the process SILENTLY,
# leaving a truncated HTML and a cheerful exit code.
#
# Usage: render_runbook_direct.sh <runbook.qmd>

set -uo pipefail

QMD="${1:?usage: render_runbook_direct.sh <runbook.qmd>}"
[ -f "$QMD" ] || { echo "no such file: $QMD" >&2; exit 2; }

QUARTO_BIN=/camp/apps/eb/software/quarto/1.5.57-x64/bin/quarto
export QUARTO_PYTHON=/nemo/lab/turnerj/working/guangxin/envs/runbook/bin/python
export PATH=/nemo/lab/turnerj/working/guangxin/envs/reanalysis_R/bin:$PATH

echo "quarto  : $("$QUARTO_BIN" --version)"
echo "python  : $QUARTO_PYTHON"
echo "pandoc  : $(command -v pandoc)"
"$QUARTO_PYTHON" -c "import bash_kernel; print('bash_kernel: OK')"

"$QUARTO_BIN" render "$QMD" --to html
echo "QUARTO_EXIT=$?"

# --------------------------------------------------------------------------
# THE EXIT-0 TRAP. Quarto embeds a failing chunk's error text in the HTML and
# still exits 0, so the exit status proves nothing. Verify the ARTEFACT.
# --------------------------------------------------------------------------
HTML="${QMD%.qmd}.html"
if [ ! -f "$HTML" ]; then
    echo "VERIFY: FAIL -- no HTML produced" >&2
    exit 1
fi

# Scan CHUNK OUTPUT only -- and specifically NOT the echoed chunk source.
#
# Two false-positive classes bit this on the way here, both worth naming:
#   1. Grepping the whole HTML matches quarto's inlined ClipboardJS, whose
#      minified source contains "Error" and "Invalid" (embed-resources puts the
#      whole library in the file).
#   2. Grepping every <pre> matches the DISPLAYED SOURCE of the chunks. This
#      runbook's own scripts contain the literal strings 'MISMATCH' and
#      'NO FILES' as the failure branches they print -- so the source echo trips
#      a grep looking for exactly those.
#
# Quarto marks source blocks with class="sourceCode" and output blocks with
# class="cell-output*", so keep only the latter.
CHUNKTXT=$(mktemp)
python3 - "$HTML" > "$CHUNKTXT" <<'PY'
import html, re, sys
doc = open(sys.argv[1], encoding='utf-8', errors='replace').read()
for m in re.finditer(r'<pre([^>]*)>(.*?)</pre>', doc, flags=re.S):
    attrs, body = m.group(1), m.group(2)
    if 'sourceCode' in attrs:          # echoed chunk source, not its output
        continue
    print(html.unescape(re.sub(r'<[^>]+>', '', body)))
PY

nerr=$(grep -ciE 'No such file or directory|command not found|Traceback \(most recent|MISSING /|NO FILES|MISMATCH|^FAIL' "$CHUNKTXT" || true)
nsent=$(grep -c 'PHASE_COMPLETE' "$CHUNKTXT" || true)
npass=$(grep -c 'PASS -- every number above is reproduced' "$CHUNKTXT" || true)
nmatch=$(grep -c 'MATCH' "$CHUNKTXT" || true)

echo ""
echo "VERIFY on $HTML  (chunk output only, $(wc -l < "$CHUNKTXT") lines)"
echo "  error signatures      : $nerr   (expect 0)"
echo "  PHASE_COMPLETE        : $nsent  (expect >=1; absent means it died early)"
echo "  predicate self-check  : $npass  (expect 1: reproduces step3.log)"
echo "  21.39% self-check     : $nmatch (expect >=1: MATCH)"
echo "  size                  : $(stat -c %s "$HTML") bytes"

rc_verify=0
if [ "$nsent" -lt 1 ]; then
    echo "VERIFY: FAIL -- sentinel absent, document did not run to the end" >&2
    rc_verify=1
fi
if [ "$nerr" -gt 0 ]; then
    echo "VERIFY: FAIL -- error signatures in chunk output:" >&2
    grep -inE 'No such file or directory|command not found|Traceback \(most recent|MISSING /|NO FILES|MISMATCH|^FAIL' \
        "$CHUNKTXT" | head -8 >&2
    rc_verify=1
fi
if [ "$npass" -lt 1 ]; then
    echo "VERIFY: FAIL -- the predicate self-check against step3.log did not pass" >&2
    rc_verify=1
fi

if [ "$rc_verify" -eq 0 ]; then
    echo "VERIFY: PASS -- document ran to the end, self-checks passed, no error signatures"
fi
rm -f "$CHUNKTXT"
exit $rc_verify
