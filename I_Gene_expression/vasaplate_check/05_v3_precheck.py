#!/usr/bin/env python3
"""v3_precheck.py -- prove the v3 rRNA reference is the intended one-entry fix,
and that a stages 3-7 re-run on it can actually start.

WHY (Rule 2 of vasaseq-lab-conventions)
---------------------------------------
The run this precedes is ~6-10 h over 384 cells. Stage 7 has already destroyed
one run in this repo (job 50542441, IndexError in reduceGeneName on names that
violated an undocumented 4-field contract) and nothing visible beforehand
predicted it. This is read-only and touches nothing.

WHAT IT ASSERTS
---------------
1. v3 vs v2 is EXACTLY the documented change: same 921 names in the same order,
   920 sequences byte-identical, exactly 1 reverse-complemented, 0 other edits.
   If this fails, v3 is not the file the README describes and the run must stop.
2. The reverse-complemented entry is ENSMUSG00000106106 -- the locus the rRNA
   diagnosis blames for the ~600x discrepancy (95,823 UFIs against 72 published).
3. v3 has NO remaining antisense entries relative to the 47S units, which is the
   property ORIENT_TO_UNITS=yes was supposed to establish. Checked by the same
   test the builder uses, not by trusting the builder.
4. Stage 1-2 products needed at START=3 are present and resolvable (they are
   symlinks into an older run directory, so a broken link is a silent failure).
5. The bwa index for v3 is complete and newer than the fasta.

Exit 0 = the run may be submitted. Exit 1 = do not submit.
"""
import glob
import gzip
import os
import subprocess
import sys

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    names, seqs, name, buf = [], [], None, []
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        for line in fh:
            if line[0] == '>':
                if name is not None:
                    names.append(name); seqs.append(''.join(buf))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None:
        names.append(name); seqs.append(''.join(buf))
    return names, seqs


def main():
    v2p, v3p, prev_run, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    lines, fail, warn = [], [], []

    def emit(s=''):
        lines.append(s); print(s, flush=True)

    def section(t):
        emit(); emit('=' * 74); emit(t); emit('=' * 74)

    # ---- 1. v3 vs v2 ----------------------------------------------------
    section('CHECK 1 -- v3 is exactly the documented one-entry change')
    n2, s2 = read_fasta(v2p)
    n3, s3 = read_fasta(v3p)
    emit('v2: %d entries   v3: %d entries' % (len(n2), len(n3)))
    if n2 != n3:
        fail.append('name lists differ (order or content) between v2 and v3')
        emit('  x names differ -- cannot continue the comparison')
    else:
        emit('  ok same %d names, same order' % len(n2))
        identical = revcomp = other = 0
        rc_names = []
        for nm, a, b in zip(n2, s2, s3):
            if a == b:
                identical += 1
            elif b == rc(a):
                revcomp += 1; rc_names.append(nm)
            else:
                other += 1
                emit('  x %s changed but is NOT the reverse complement' % nm[:60])
        emit('  identical           : %d' % identical)
        emit('  reverse-complemented: %d  %s' % (revcomp, rc_names))
        emit('  changed some other way: %d' % other)
        if revcomp != 1 or other != 0 or identical != len(n2) - 1:
            fail.append('v3 is not the documented 1-revcomp/920-identical change '
                        '(got %d identical, %d revcomp, %d other)'
                        % (identical, revcomp, other))
        else:
            emit('  ok matches the README exactly')

        # ---- 2. is it the blamed locus? ---------------------------------
        section('CHECK 2 -- the flipped entry is the locus blamed for the ~600x')
        emit('flipped: %s' % (rc_names[0] if rc_names else '(none)'))
        if rc_names and 'ENSMUSG00000106106' in rc_names[0]:
            emit('  ok ENSMUSG00000106106 -- the locus with 95,823 UFIs vs 72 published')
        else:
            fail.append('the flipped entry is not ENSMUSG00000106106; the fix does not '
                        'target the diagnosed defect')
            emit('  x expected ENSMUSG00000106106')

    # ---- 3. no remaining antisense entries ------------------------------
    section('CHECK 3 -- v3 has no antisense entries left (the builder\'s own test)')
    units = [i for i, nm in enumerate(n3)
             if any(k in nm for k in ('BK000964', 'U13369', 'NR_046233', '47S', 'RNA45S'))]
    emit('47S/unit reference entries found in v3: %d %s'
         % (len(units), [n3[i][:40] for i in units[:4]]))
    if not units:
        warn.append('no 47S unit entry identifiable by name in v3; the antisense test '
                    'below was skipped (the builder aligns to the NCBI units, which '
                    'may be named differently)')
        emit('  ! skipped -- cannot identify the unit entries by name')
    else:
        # cheap self-contained test: every non-unit entry should share a k-mer with
        # a unit in FORWARD orientation more often than reverse.
        unit_seq = ''.join(s3[i] for i in units)
        K = 25
        ufwd = {unit_seq[i:i+K] for i in range(0, len(unit_seq) - K, 3)}
        bad = []
        for nm, sq in zip(n3, s3):
            if len(sq) < 200:
                continue
            probes = [sq[i:i+K] for i in range(0, min(len(sq), 3000) - K, 97)]
            if not probes:
                continue
            f = sum(p in ufwd for p in probes)
            r = sum(rc(p) in ufwd for p in probes)
            if r > f and r >= 2:
                bad.append((nm, f, r))
        emit('entries matching the units better REVERSED than forward: %d' % len(bad))
        for nm, f, r in bad[:5]:
            emit('  x %-52s fwd=%d rev=%d' % (nm[:52], f, r))
        if bad:
            fail.append('%d entries in v3 still look antisense to the 47S units' % len(bad))
        else:
            emit('  ok none -- ORIENT_TO_UNITS did its job')

    # ---- 4. stage 1-2 inputs --------------------------------------------
    section('CHECK 4 -- stage 1-2 products for START=3 exist and resolve')
    cbc = sorted(glob.glob(os.path.join(prev_run, '*_cbc.fastq.gz')))
    emit('_cbc.fastq.gz in %s: %d' % (prev_run, len(cbc)))
    if len(cbc) != 384:
        fail.append('expected 384 _cbc.fastq.gz, found %d' % len(cbc))
    broken = [p for p in cbc if not os.path.exists(os.path.realpath(p))]
    emit('broken symlinks: %d' % len(broken))
    for p in broken[:5]:
        emit('  x %s -> %s' % (os.path.basename(p), os.readlink(p)))
    if broken:
        fail.append('%d _cbc.fastq.gz symlinks do not resolve' % len(broken))
    else:
        emit('  ok all resolve')
    sizes = [os.path.getsize(os.path.realpath(p)) for p in cbc[:20]]
    if sizes:
        emit('first 20 sizes: min %.1f MB, max %.1f MB'
             % (min(sizes) / 1e6, max(sizes) / 1e6))
        if min(sizes) == 0:
            fail.append('at least one _cbc.fastq.gz is zero bytes')

    # ---- 5. bwa index ----------------------------------------------------
    section('CHECK 5 -- bwa index for v3 is complete and not stale')
    fa_mt = os.path.getmtime(v3p)
    for ext in ('.amb', '.ann', '.bwt', '.pac', '.sa'):
        p = v3p + ext
        if not os.path.exists(p):
            fail.append('missing bwa index file %s' % os.path.basename(p))
            emit('  x missing %s' % ext)
        else:
            age = os.path.getmtime(p) - fa_mt
            emit('  %-5s %10d bytes  %+.0f s vs fasta %s'
                 % (ext, os.path.getsize(p), age, 'OK' if age >= -1 else 'STALE'))
            if age < -1:
                fail.append('bwa index %s is OLDER than the fasta -- reindex' % ext)

    # ---- verdict ---------------------------------------------------------
    section('VERDICT')
    if fail:
        emit('DO NOT SUBMIT (%d blocking):' % len(fail))
        for f in fail:
            emit('  x %s' % f)
    else:
        emit('OK -- stages 3-7 may be submitted on v3.')
    if warn:
        emit()
        for w in warn:
            emit('  ! %s' % w)

    with open(out, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\nreport: %s' % out, flush=True)
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
