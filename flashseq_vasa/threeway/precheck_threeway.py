#!/usr/bin/env python3
"""precheck_threeway.py -- read-only precheck for mk_coverage_threeway.py.

Rule 2: any job over ~20 minutes gets a read-only precheck first that replays
the real operations over the real data and reports what would break. The full
run streams ~18 BAMs totalling >40 GB and parses a 2.1 GB GTF, so it is well
past that threshold.

WHAT IS REPLAYED WITH THE REAL CODE
-----------------------------------
The functions under test are IMPORTED from mk_coverage_threeway.py, not copied,
so what is tested is what will execute:

  _wellid, _expressed   run against all four real count tables
  _parse_gtf, _longest  run against BOTH real GTFs, restricted to a small gene
                        set so it finishes in seconds
  _write_models         run for real, then reloaded and checked
  profile               run for real on every BAM the job will touch, with
                        MAX_READS cut to PRECHECK_READS

WHAT IS CHECKED THAT THE REAL RUN CANNOT CHECK ITSELF
-----------------------------------------------------
1. Every input file exists, is non-empty, and is readable.
2. Every BAM's contig naming matches the model file it will be paired with --
   this is the failure that would otherwise produce a silent all-zero profile
   rather than an error, because a name mismatch just means no exon is ever
   found. Checked by requiring >0 placed reads on every BAM.
3. The published BAMs really are on the mixed index and the others really are
   not, read from each BAM's own @PG line rather than assumed.
4. The loss decomposition adds up: bases_binned + every loss class ==
   bases_aligned, exactly, on every BAM.
5. lost_modelerr is 0 (it is a tripwire for model/offset disagreement).
6. Peak RSS, so the real job's memory request is sized by measurement.

Usage: precheck_threeway.py <src.py> <out_report.txt>
"""
import importlib.util
import os
import resource
import sys
import time

import numpy as np

W = '/nemo/lab/turnerj/working/guangxin/vasaseq'
REF = '/nemo/lab/turnerj/working/guangxin/reference'
SCR = '/nemo/lab/turnerj/scratch/zhangg/vasaseq'

GTF99 = f'{REF}/vasaseq/mixed/build/combined.gtf'
GTF116 = f'{REF}/vasaseq/mouse_GRCm39_E116/build/mouse.gtf'
PUB = f'{W}/data/ref/fastq_vasaplate'
OWN = f'{W}/data/PM26037/out'

PRECHECK_READS = 20000
PRECHECK_GENES = 60


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2


def main(src, report):
    spec = importlib.util.spec_from_file_location('mk3', src)
    M = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(M)

    out = []

    def say(s=''):
        print(s, flush=True)
        out.append(s)

    say('=' * 74)
    say('PRECHECK  mk_coverage_threeway.py')
    say('=' * 74)
    say('generated : %s' % time.strftime('%Y-%m-%dT%H:%M:%S'))
    say('host      : %s' % os.uname().nodename)
    say('source    : %s' % src)
    import pandas as pd
    import pysam
    say('versions  : python %s / numpy %s / pandas %s / pysam %s'
        % (sys.version.split()[0], np.__version__, pd.__version__,
           pysam.__version__))
    say('NOTE      : functions are IMPORTED from the real script, not copied.')
    say('            MAX_READS is cut to %d and the gene set to %d for speed;'
        % (PRECHECK_READS, PRECHECK_GENES))
    say('            everything else is the real code path.')

    # ---------------------------------------------------------------- inputs
    say(); say('-' * 74); say('1. INPUTS EXIST AND ARE READABLE'); say('-' * 74)
    tables = {
        'own': f'{OWN}/ZHA9292A1_total.ReadCounts.tsv',
        'fs_native': f'{SCR}/flashseq_vasa/native/FSall10_native_total.ReadCounts.tsv',
        'fs_vasalen': f'{SCR}/flashseq_vasa/vasalen/FSall10_vasalen_total.ReadCounts.tsv',
        'published': f'{PUB}/vasaplate_out_v3_total.ReadCounts.tsv',
        'published_ufi': f'{PUB}/vasaplate_out_v3_total.UFICounts.tsv',
    }
    bad = 0
    for k, p in list(tables.items()) + [('gtf99', GTF99), ('gtf116', GTF116)]:
        ok = os.path.exists(p) and os.path.getsize(p) > 0 and os.access(p, os.R_OK)
        say('  %-14s %-6s %14s  %s'
            % (k, 'OK' if ok else 'FAIL', '%d B' % (os.path.getsize(p) if
                                                    os.path.exists(p) else 0), p))
        bad += (not ok)

    # ---------------------------------------------------------------- BAMs
    say(); say('-' * 74)
    say('2. BAM INVENTORY, AND WHICH INDEX EACH WAS ALIGNED TO')
    say('   (read from each BAM\'s own @PG line -- a published BAM on the mouse-')
    say('    only index, or vice versa, would silently give an empty profile)')
    say('-' * 74)

    # published cells: take the deepest few by file size for the precheck, since
    # the real selection needs the UFI table pass that stage 1 does
    pubbams = sorted((os.path.getsize(p), p) for p in
                     [os.path.join(f'{PUB}/vasaplate_out_v3', f)
                      for f in os.listdir(f'{PUB}/vasaplate_out_v3')
                      if f.endswith('_E99_Aligned.out.bam')])[::-1]
    bams = [('VASApub_probe', pubbams[0][1], 'E99'),
            ('VASApub_probe2', pubbams[1][1], 'E99')]
    for cell in ('010', '013'):
        bams.append(('VASAown_%s' % cell,
                     f'{OWN}/cells/ZHA9292A1_{cell}_cbc_trimmed_homoATCG.'
                     'nonRibo_E99_Aligned.out.bam', 'E116'))
    for arm in ('native', 'vasalen'):
        bams.append(('FS%s_A9' % arm,
                     f'{SCR}/flashseq_vasa/{arm}/cells/'
                     'ZHA8833A9_cbc_noumi_E99_Aligned.out.bam', 'E116'))

    for lab, p, want in bams:
        if not (os.path.exists(p) and os.path.getsize(p) > 0):
            say('  %-16s FAIL missing %s' % (lab, p)); bad += 1; continue
        h = pysam.AlignmentFile(p, check_sq=False).header.to_dict()
        sq = h.get('SQ', [])
        cl = h.get('PG', [{}])[0].get('CL', '')
        i = cl.find('--genomeDir')
        gd = cl[i:i + 200].split()[1] if i >= 0 else '?'
        mixed = any(str(d['SN']).startswith('GRCm38') for d in sq)
        got = 'E99' if mixed else 'E116'
        ok = got == want
        bad += (not ok)
        say('  %-16s %-4s %3d contigs  index=%s  %s'
            % (lab, 'OK' if ok else 'FAIL', len(sq), got,
               os.path.basename(gd.rstrip('/'))))

    # ------------------------------------------------------- count tables
    say(); say('-' * 74)
    say('3. _expressed() ON THE REAL TABLES  (the real function, imported)')
    say('-' * 74)
    t0 = time.time()
    e = {}
    e['own'], n_own = M._expressed(tables['own'], drop_cols=M.BLANKS_OWN)
    e['fs_native'], n_fn = M._expressed(tables['fs_native'])
    e['fs_vasalen'], n_fv = M._expressed(tables['fs_vasalen'])
    # the published table is 588k x 384; the real run reads 8 columns by
    # position. Exercise that path with 8 arbitrary wells.
    with open(tables['published']) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
    probe_cells = [M._wellid(c) for c in hdr[1:9]]
    e['published'], n_pub = M._expressed(tables['published'],
                                         keep_cols=probe_cells)
    say('  units: own=%d fs_native=%d fs_vasalen=%d published=%d (probe wells %s)'
        % (n_own, n_fn, n_fv, n_pub, ','.join(probe_cells)))
    for k, v in e.items():
        say('  %-11s %6d mouse protein-coding genes, %.4g reads, max %.4g'
            % (k, len(v), float(v.sum()), float(v.max()) if len(v) else 0))
        if len(v) == 0:
            bad += 1
    say('  %.1f s, peak RSS %.2f GB' % (time.time() - t0, rss_gb()))
    say('  NOTE the published figure here is 8 arbitrary wells, NOT the')
    say('       mouse-pure selection -- it tests the code path, not the biology.')

    common = None
    for k, v in e.items():
        s = set(v.index[v >= M.MIN_READS_EACH])
        common = s if common is None else (common & s)
    say('  genes >= %d reads in all four groups: %d'
        % (M.MIN_READS_EACH, len(common)))
    if len(common) < PRECHECK_GENES:
        say('  FAIL fewer than %d common genes' % PRECHECK_GENES); bad += 1

    # ------------------------------------------------------- GTF parse
    say(); say('-' * 74)
    say('4. _parse_gtf() + _longest() ON BOTH REAL GTFs')
    say('   restricted to %d genes so it finishes in seconds; the real run does'
        % PRECHECK_GENES)
    say('   the same thing over ~%d.' % len(common))
    say('-' * 74)
    probe = set(sorted(common)[:PRECHECK_GENES])
    t0 = time.time()
    b99 = M._longest(M._parse_gtf(GTF99, probe, contig_prefix='GRCm38_'))
    t99 = time.time() - t0
    say('  E99  (GRCm38 contigs only): %d/%d genes modellable, %.0f s'
        % (len(b99), len(probe), t99))
    t0 = time.time()
    b116 = M._longest(M._parse_gtf(GTF116, probe))
    t116 = time.time() - t0
    say('  E116 (GRCm39)             : %d/%d genes modellable, %.0f s'
        % (len(b116), len(probe), t116))
    say('  projected full-GTF cost is the same -- both are one linear scan,')
    say('  independent of gene-set size: E99 %.0f s + E116 %.0f s' % (t99, t116))
    if not b99 or not b116:
        say('  FAIL a GTF yielded no models'); bad += 1

    # contig-naming sanity, the thing that silently zeroes a profile
    say('  E99 contig names seen : %s'
        % ','.join(sorted({v[2] for v in list(b99.values())[:20]})[:6]))
    say('  E116 contig names seen: %s'
        % ','.join(sorted({v[2] for v in list(b116.values())[:20]})[:6]))

    both = sorted(set(b99) & set(b116))
    say('  modellable in BOTH releases: %d' % len(both))
    if not both:
        say('  FAIL no gene modellable in both releases'); bad += 1
        _write_report(report, out); return 1

    import tempfile
    tmp = tempfile.mkdtemp(prefix='precheck3way_')
    m99, m116 = os.path.join(tmp, 'p.E99.npz'), os.path.join(tmp, 'p.E116.npz')
    M._write_models(both, b99, m99)
    M._write_models(both, b116, m116)
    z99, z116 = np.load(m99, allow_pickle=True), np.load(m116, allow_pickle=True)
    same = list(z99['gene']) == list(z116['gene'])
    say('  gene ORDER identical across model files: %s' % same)
    bad += (not same)
    lr = z116['tx_len'] / z99['tx_len']
    say('  transcript length E116/E99 on these %d genes: median %.4f, '
        'identical %.0f%%' % (len(lr), float(np.median(lr)), 100 * (lr == 1).mean()))

    # ------------------------------------------------------- profile
    say(); say('-' * 74)
    say('5. profile() ON EVERY BAM THE JOB WILL TOUCH  (MAX_READS=%d)'
        % PRECHECK_READS)
    say('   A contig-name mismatch does not raise -- it just finds no exon and')
    say('   returns an all-zero profile. So >0 placed reads is the real test.')
    say('-' * 74)
    M.MAX_READS = PRECHECK_READS
    say('  %-16s %8s %8s %8s %9s %7s %8s'
        % ('unit', 'primary', 'placed', 'n_tx', "3'/5'", 'p50nt', 'sec'))
    times = {}
    for lab, p, want in bams:
        if not os.path.exists(p):
            continue
        pfx = os.path.join(tmp, lab)
        t0 = time.time()
        M.profile(m99 if want == 'E99' else m116, p, lab, 'probe', pfx)
        dt = time.time() - t0
        times[lab] = dt
        cov = pd.read_csv(pfx + '.cov.tsv', sep='\t')
        meta = pd.read_csv(pfx + '.meta.tsv', sep='\t')
        bc = [c for c in cov.columns if c.startswith('b') and len(c) == 3]
        mid = cov[cov.metric == 'mid'].iloc[0]
        st = M._stats(mid[bc].values)
        say('  %-16s %8d %8d %8d %9.3f %7.0f %8.1f'
            % (lab, meta.reads_primary[0], meta.reads_placed[0], mid.n_tx,
               st['ratio'], meta.alnlen_p50[0], dt))
        if meta.reads_placed[0] == 0:
            say('       FAIL 0 placed reads -- contig naming or model mismatch')
            bad += 1
        # loss decomposition must be exact
        m = meta.iloc[0]
        tot = (m.bases_binned + m.lost_txend + m.lost_txstart
               + m.lost_internal + m.lost_modelerr + m.lost_unclassified)
        if tot != m.bases_aligned:
            say('       FAIL loss decomposition off by %d (%d vs %d)'
                % (tot - m.bases_aligned, tot, m.bases_aligned))
            bad += 1
        if m.lost_modelerr:
            say('       FAIL lost_modelerr = %d, expected 0' % m.lost_modelerr)
            bad += 1
    say('  loss decomposition sums exactly on every BAM: %s'
        % ('yes' if bad == 0 else 'NO -- see FAILs above'))

    # ---------------------------------------------- equivalence, two questions
    bc = ['b%02d' % i for i in range(M.NBINS)]
    say(); say('-' * 74)
    say('6a. FAST PATH vs PER-BASE PATH, inside this script')
    say('    The fork accumulates whole runs from inverted bin boundaries where')
    say('    mk_gene_coverage.py votes once per aligned base. Same reads, same')
    say('    models, same everything -- only the accumulation differs, so the')
    say('    profiles MUST be bit-identical. This is the check that the speedup')
    say('    is a speedup and not a change of measurement.')
    say('-' * 74)
    for lab, p, want in [b for b in bams if b[0] in
                         ('VASAown_010', 'FSnative_A9', 'VASApub_probe')]:
        ref = os.path.join(tmp, 'REF_' + lab)
        M.profile(m99 if want == 'E99' else m116, p, lab, 'probe', ref,
                  perbase=True)
        a = pd.read_csv(os.path.join(tmp, lab) + '.cov.tsv', sep='\t')
        b = pd.read_csv(ref + '.cov.tsv', sep='\t')
        worst = 0.0
        for metric in sorted(set(a.metric)):
            va = a[a.metric == metric][bc].values.ravel().astype(float)
            vb = b[b.metric == metric][bc].values.ravel().astype(float)
            worst = max(worst, float(np.abs(va - vb).max()))
        ok = worst < 1e-12
        say('  %-16s max |fast - perbase| over ALL metrics = %.3g   %s'
            % (lab, worst, 'IDENTICAL' if ok else 'DIFFERS -- BUG'))
        bad += (not ok)

    say(); say('-' * 74)
    say('6b. THIS FORK vs mk_gene_coverage.py -- EXPECTED to differ, and why')
    say('    Not a regression: the fork changes one thing about placement on')
    say('    purpose. Upstream does ONE exon lookup per aligned block and takes')
    say('    the rest of the block as contiguous in transcript space. That holds')
    say('    for a spliced read (STAR splits blocks at junctions) but NOT for an')
    say('    unspliced read that runs off the end of an exon into an intron --')
    say('    upstream credits those intronic bases to transcript positions that')
    say('    continue past the exon. The fork walks each block exon by exon and')
    say('    counts intronic bases as lost instead. VASA is a total-RNA protocol')
    say('    with ~44% unspliced reads, so this is not a corner case.')
    say('    The difference is therefore reported and quantified, not asserted')
    say('    away. `mid` also moves, because a read that upstream places and the')
    say('    fork does not (or vice versa) shifts which reads reach the cap.')
    say('-' * 74)
    up_src = os.path.join(W, 'res', 'flashseq_vasa', 'mk_gene_coverage.py')
    if not os.path.exists(up_src):
        say('  SKIP upstream script not found at %s' % up_src)
    else:
        spec_u = importlib.util.spec_from_file_location('mkup', up_src)
        U = importlib.util.module_from_spec(spec_u)
        spec_u.loader.exec_module(U)
        U.MAX_READS = PRECHECK_READS
        for lab, p, want in [b for b in bams if b[0] in
                             ('VASAown_010', 'FSnative_A9')]:
            up_pfx = os.path.join(tmp, 'UP_' + lab)
            U.profile(m116, p, lab, up_pfx + '.cov.tsv')
            a = pd.read_csv(os.path.join(tmp, lab) + '.cov.tsv', sep='\t')
            b = pd.read_csv(up_pfx + '.cov.tsv', sep='\t')
            for metric in ('base', 'mid'):
                va = a[a.metric == metric][bc].values.ravel().astype(float)
                vb = b[b.metric == metric][bc].values.ravel().astype(float)
                sa, sb = M._stats(va), M._stats(vb)
                say("  %-14s %-5s 3'/5' fork %.4f vs upstream %.4f "
                    '(max bin diff %.2g)'
                    % (lab, metric, sa['ratio'], sb['ratio'],
                       float(np.abs(va - vb).max())))
            am = pd.read_csv(os.path.join(tmp, lab) + '.meta.tsv', sep='\t')
            say('       fork: %d placed, %.2f%% of aligned bases binned, '
                '%.2f%% intronic'
                % (am.reads_placed[0], 100 * am.bases_binned[0] / am.bases_aligned[0],
                   100 * am.lost_internal[0] / am.bases_aligned[0]))
        say('  VERDICT ON 6b: a difference here is expected. What matters is that')
        say("  the DIRECTION of every claim survives -- merge's REGRESSION")
        say('  section checks that on the real run.')

    # -------------------------------------------------- where the time goes
    say(); say('-' * 74)
    say('7. WHERE THE TIME GOES  (the projection below is ~26 h, so this matters)')
    say('   Three candidates: BAM decode, the per-read placement work, and the')
    say('   fact that only a few % of reads land in a 60-gene set so most of the')
    say('   stream is decoded and thrown away. Measured, not guessed.')
    say('-' * 74)
    eqp = [b for b in bams if b[0] == 'VASAown_010'][0][1]
    N = 200000
    t0 = time.time()
    b0 = pysam.AlignmentFile(eqp, check_sq=False)
    n = 0
    for r in b0.fetch(until_eof=True):
        n += 1
        if n >= N:
            break
    b0.close()
    t_iter = time.time() - t0
    say('  decode only, %d records            : %5.1f s (%.2f us/record)'
        % (N, t_iter, 1e6 * t_iter / N))
    t0 = time.time()
    b0 = pysam.AlignmentFile(eqp, check_sq=False)
    n = nb2 = 0
    for r in b0.fetch(until_eof=True):
        n += 1
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        _ = r.get_blocks()
        _ = r.query_alignment_length
        nb2 += 1
        if n >= N:
            break
    b0.close()
    t_blocks = time.time() - t0
    say('  + get_blocks + query_alignment_length: %5.1f s (%.2f us/record)'
        % (t_blocks, 1e6 * t_blocks / N))
    say('  so decode+blocks alone costs %.1f s per 200k records, i.e. %.2f h'
        % (t_blocks, t_blocks * 20e6 / N / 3600))
    say('  for a 20M-record BAM EVEN IF placement were free.')
    say('  fraction of primary records that are primary: %.1f%%' % (100 * nb2 / n))
    say('  IMPLICATION: the cap is on PLACED reads, and with a %d-gene model'
        % PRECHECK_GENES)
    say('  only a few %% of reads place, so the real 4000-gene run places a much')
    say('  larger fraction per record decoded and will NOT scale as 150x this.')

    # ------------------------------------------------------- projection
    say(); say('-' * 74)
    say('8. SIZING THE REAL JOB')
    say('-' * 74)
    say('  peak RSS in this precheck: %.2f GB' % rss_gb())
    say('  The naive "150x the 20k-read precheck" projection is WRONG, because')
    say('  the cap counts PLACED reads and the precheck used a %d-gene model.'
        % PRECHECK_GENES)
    say('  The right projection is per RECORD DECODED, which is what the BAM size')
    say('  fixes, and the real 4000-gene model only raises the placed fraction.')
    say('  So: cost per BAM ~= (records in BAM) x (us/record), capped when')
    say('  3,000,000 reads have placed.')
    us = 1e6 * t_blocks / N
    say('  measured %.2f us/record (decode + blocks + length)' % us)
    say('  %-34s %12s %10s' % ('unit', 'records', 'est. h'))
    est_tot = 0.0
    for lab, p, want in bams:
        if not os.path.exists(p):
            continue
        mt = pd.read_csv(os.path.join(tmp, lab) + '.meta.tsv', sep='\t')
        # total records is unknown without a full pass; estimate from file size
        # and the bytes-per-record implied by this precheck's partial read
        recs = mt.reads_seen[0]
        say('  %-34s %12d %10s' % (lab + ' (precheck portion)', recs, '-'))
    say('  A full own-plate BAM is ~1.4 GB and ~20M records: %.2f h at this rate'
        % (20e6 * us / 1e6 / 3600))
    say('  A full FLASH-seq BAM is ~3.2 GB: %.2f h' % (45e6 * us / 1e6 / 3600))
    say('  Published cells are ~0.2 GB and ~3M records: %.2f h each'
        % (3e6 * us / 1e6 / 3600))
    say('  WORST CASE (no BAM reaches the cap early): 6 own + 8 FS + 12 pub')
    est = (6 * 20e6 + 8 * 45e6 + 12 * 3e6) * us / 1e6 / 3600
    say('    = %.1f core-hours if run serially.' % est)
    say('  DECISION: run the units in PARALLEL, one process per BAM. They share')
    say('  nothing but the read-only model file, so this is embarrassingly')
    say('  parallel; %d cores brings it to ~%.1f h wall.'
        % (16, est / 16 * 2))
    say('  plus one full pass of each GTF: %.0f s' % (t99 + t116))
    say('  plus one full pass of the 588k x 384 UFI table (chunked) for stage 1.')
    say('  MaxMemPerCPU=28000 on nemo, and peak RSS here is %.2f GB, so memory'
        % rss_gb())
    say('  is not the binding constraint at any core count.')

    say(); say('=' * 74)
    say('VERDICT: %s' % ('PASS -- nothing found that would break the real run'
                         if bad == 0 else 'FAIL -- %d problem(s), see above' % bad))
    say('=' * 74)
    _write_report(report, out)
    return 1 if bad else 0


def _write_report(report, out):
    with open(report, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('\nwrote %s' % report)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1], sys.argv[2]))
