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

# Pin LC_NUMERIC so the rendered numbers are reproducible. Added 2026-07-30.
# The document's chunks use printf "%'d", whose thousands separator comes from
# the locale: under a login shell it gave "1,154,910", under a batch job
# (LANG=C.UTF-8) the same chunk printed "1154910". Same value, different text,
# so a diff of two renders was full of noise that looked like drift.
# Only LC_NUMERIC is set -- NOT LC_ALL, because LC_COLLATE would change `sort`
# order inside the chunks and that WOULD alter results.
export LC_NUMERIC=en_US.utf8

echo "quarto  : $("$QUARTO_BIN" --version)"
echo "python  : $QUARTO_PYTHON"
echo "pandoc  : $(command -v pandoc)"
"$QUARTO_PYTHON" -c "import bash_kernel; print('bash_kernel: OK')"

# --------------------------------------------------------------------------
# WARM THE SHELL-STARTUP CACHE BEFORE RENDERING. Added 2026-07-30.
#
# bash_kernel starts bash through pexpect with pexpect's own bashrc.sh as
# --rcfile, and that file SOURCES ~/.bashrc. Our ~/.bashrc runs `conda init`,
# which reads a large tree under /camp/apps/eb/software/Anaconda3/2024.10-1/.
# On a compute node with a cold page cache that takes ~45 s -- and pexpect's
# prompt timeout is 30 s. The kernel dies before replying to kernel_info and
# quarto exits 1 without writing anything.
#
# Measured on cn087, job 51055067, three consecutive real kernel boots:
#     1st (cold) 46.8 s -> FAILED
#     2nd (warm) 22.9 s -> ok
#     3rd (warm)  8.9 s -> ok
#
# It is purely a cold-cache effect, so touching the same path once first is
# enough. This is why the document renders on a login node (cache always warm
# there) but not on a fresh compute node -- and rendering on a login node is
# exactly what the runbook forbids.
PEXPECT_RC=$("$QUARTO_PYTHON" -c \
    "import pexpect,os;print(os.path.join(os.path.dirname(pexpect.__file__),'bashrc.sh'))")
echo "warmup  : $PEXPECT_RC"
for attempt in 1 2 3; do
    t0=$(date +%s)
    bash --rcfile "$PEXPECT_RC" -i -c 'true' >/dev/null 2>&1 || true
    el=$(( $(date +%s) - t0 ))
    echo "warmup  : attempt $attempt took ${el}s (need < 30s for pexpect)"
    [ "$el" -lt 20 ] && break
done
if [ "$el" -ge 30 ]; then
    echo "VERIFY: FAIL -- shell startup still ${el}s after $attempt warmups;" \
         "pexpect will time out at 30s and the kernel will die." >&2
    exit 1
fi

HTML="${QMD%.qmd}.html"

# STALE-ARTEFACT GUARD, added 2026-07-30 after this script reported
# "VERIFY: PASS" on a render that never happened.
#
# What went wrong: quarto died with `Kernel died before replying to kernel_info`
# and exited 1, so no new HTML was written -- but the PREVIOUS run's HTML was
# still sitting there. Every check below then ran against that file, found the
# sentinel and the self-checks (because the OLD render really had passed), and
# printed PASS. Job 51052670.
#
# So record what was on disk before, and refuse to verify a file the render did
# not touch. "The output file exists" is not evidence the stage ran.
HTML_BEFORE=""
[ -f "$HTML" ] && HTML_BEFORE=$(stat -c %Y "$HTML")

"$QUARTO_BIN" render "$QMD" --to html
QUARTO_EXIT=$?
echo "QUARTO_EXIT=$QUARTO_EXIT"

if [ ! -f "$HTML" ]; then
    echo "VERIFY: FAIL -- no HTML produced" >&2
    exit 1
fi

HTML_AFTER=$(stat -c %Y "$HTML")
if [ -n "$HTML_BEFORE" ] && [ "$HTML_BEFORE" = "$HTML_AFTER" ]; then
    echo "VERIFY: FAIL -- $HTML was NOT rewritten (mtime unchanged:" \
         "$(date -d @"$HTML_AFTER" '+%F %T')). This is a STALE artefact from an" \
         "earlier run; nothing below it would mean anything." >&2
    exit 1
fi

# A non-zero quarto exit is a real failure even though a zero one proves nothing.
if [ "$QUARTO_EXIT" -ne 0 ]; then
    echo "VERIFY: FAIL -- quarto exited $QUARTO_EXIT" >&2
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

# The completion sentinel is whatever the DOCUMENT emits, not a name hardcoded
# here. Hardcoding it made a full render "fail" once because that runbook echoes
# RUNBOOK_COMPLETE while this script was looking for PHASE_COMPLETE -- a verifier
# that reports a false failure is as bad as one that reports a false pass.
# Convention: the last chunk echoes a bare <SOMETHING>_COMPLETE token.
SENTINEL=$(grep -oE '\b[A-Z][A-Z0-9_]*_COMPLETE\b' "$QMD" | sort -u | head -1)
SENTINEL=${SENTINEL:-PHASE_COMPLETE}

nerr=$(grep -ciE 'No such file or directory|command not found|Traceback \(most recent|MISSING /|NO FILES|MISMATCH|^FAIL' "$CHUNKTXT" || true)
nsent=$(grep -c "$SENTINEL" "$CHUNKTXT" || true)

echo ""
echo "VERIFY on $HTML  (chunk output only, $(wc -l < "$CHUNKTXT") lines)"
echo "  error signatures      : $nerr   (expect 0)"
echo "  sentinel '$SENTINEL' : $nsent  (expect >=1; absent means it died early)"
echo "  size                  : $(stat -c %s "$HTML") bytes"

# Self-checks are asserted only when the document declares them, so this script
# works on any runbook rather than only the one it was written against.
npass=1; nmatch=1
if grep -q 'PASS: same predicate' "$QMD"; then
    npass=$(grep -c 'PASS: same predicate' "$CHUNKTXT" || true)
    echo "  predicate vs step3.log: $npass  (expect >=1)"
fi
if grep -q "reproduce step3.log whole-library" "$QMD"; then
    nmatch=$(grep -c '\-> *MATCH' "$CHUNKTXT" || true)
    echo "  21.39% recomputation  : $nmatch (expect >=1: MATCH)"
fi

rc_verify=0
if [ "$nsent" -lt 1 ]; then
    echo "VERIFY: FAIL -- sentinel '$SENTINEL' absent; document did not run to the end" >&2
    rc_verify=1
fi
if [ "$nerr" -gt 0 ]; then
    echo "VERIFY: FAIL -- error signatures in chunk output:" >&2
    grep -inE 'No such file or directory|command not found|Traceback \(most recent|MISSING /|NO FILES|MISMATCH|^FAIL' \
        "$CHUNKTXT" | head -8 >&2
    rc_verify=1
fi
if [ "$npass" -lt 1 ]; then
    echo "VERIFY: FAIL -- the document declares a step3.log predicate check but it did not pass" >&2
    rc_verify=1
fi
if [ "$nmatch" -lt 1 ]; then
    echo "VERIFY: FAIL -- the document declares a 21.39% recomputation but it did not MATCH" >&2
    rc_verify=1
fi

if [ "$rc_verify" -eq 0 ]; then
    echo "VERIFY: PASS -- document ran to the end, self-checks passed, no error signatures"
fi
rm -f "$CHUNKTXT"
exit $rc_verify
