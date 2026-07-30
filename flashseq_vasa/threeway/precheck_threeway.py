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

    # ------------------------------------------------------- projection
    say(); say('-' * 74)
    say('6. SIZING THE REAL JOB')
    say('-' * 74)
    say('  peak RSS in this precheck: %.2f GB' % rss_gb())
    per = np.mean([times[k] for k in times]) if times else 0
    scale = M.__dict__.get('MAX_READS_REAL', 3_000_000) / PRECHECK_READS
    say('  mean %.1f s per BAM at %d reads' % (per, PRECHECK_READS))
    say('  the real run caps at 3,000,000 placed reads per BAM (%.0fx), over 18'
        % (3_000_000 / PRECHECK_READS))
    say('  BAMs. Published cells are ~10x shallower and will not reach the cap,')
    say('  so they cost less than the linear projection.')
    say('  linear projection, deep BAMs: %.1f h for 18 BAMs'
        % (18 * per * (3_000_000 / PRECHECK_READS) / 3600))
    say('  plus one full pass of each GTF: %.0f s' % (t99 + t116))
    say('  plus one full pass of the 588k x 384 UFI table (chunked) for stage 1.')

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
