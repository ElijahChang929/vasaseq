#!/usr/bin/env python3
"""
merge_pickles.py -- combine per-library step-6 dicts into the single frame
step 7 expects, doing ONLY what upstream step 6's parent tail does.

WHY THIS EXISTS
---------------
Upstream countTables_2pickle_cellsSpliced.py runs once over the whole cells/
folder with an 8-wide pool. That shape does not fit ten FLASH-seq libraries: a
full library yields roughly 870 MB of *_genes.bed.gz, and the same script was
measured on the VASA side at 45.6 GB of RSS per GB of BED with a scaling exponent
of 1.02 (69 MB -> 28m36s/3.24 GB; 254 MB -> 1h47m27s/12.22 GB). Eight libraries
resident at once is ~250 GB in a single process, and any failure loses the whole
run. So step 6 is run once per library and the results are merged here.

WHY THE SPLIT IS EXACT AND NOT AN APPROXIMATION
-----------------------------------------------
Upstream's parent, after the pool returns, is exactly:

    for (cell, cnt) in pool.imap_unordered(get_cellDict, cells):
        if len(cnt) > 0:
            gcnt[cell] = cnt
    pickle.dump(gcnt, open(output + 'dict.pickle', 'wb'))
    cntdf = pd.DataFrame(gcnt)
    pickle.dump(cntdf, open(output + '.pickle', 'wb'))
    os.system('gzip ' + output + '.pickle')

get_cellDict is per-library independent: it globs cell + '*_genes.bed.gz' and
shares no state with any other cell. So the union of per-library dicts IS the
gcnt a combined run would have built, and pd.DataFrame over it is the same frame.
This script reproduces those five lines and nothing else -- no re-counting, no
re-aggregation, no filtering.

Three properties make the equivalence hold, and all three are ASSERTED rather
than assumed:

  1. Column keys are disjoint across inputs. Each per-library run produced one
     key, 'cells/<LIB>'. A collision would mean two runs claimed the same
     library, and silently merging them would double-count reads.
  2. Every input dict is non-empty. Upstream skips an empty cnt (`if len(cnt) > 0`),
     so an empty input here means that library produced no assignments at all --
     a finding, not something to quietly drop.
  3. Row keys are unioned, never merged. pd.DataFrame over a dict of dicts aligns
     on the union of inner keys and fills missing entries with NaN; step 7's
     counters all test `type(x) == dict` and return 0 otherwise, so NaN is the
     correct "not detected" representation. No per-gene value is ever combined
     across libraries, so there is no summation semantics to get wrong.

Usage: merge_pickles.py <out_prefix> <lib1dict.pickle> <lib2dict.pickle> ...

Writes <out_prefix>dict.pickle and <out_prefix>.pickle.gz, the two files step 7
reads (it takes the .pickle.gz).
"""
import gzip
import os
import pickle
import sys

import pandas as pd


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    outpfx = sys.argv[1]
    inputs = sys.argv[2:]

    gcnt = {}
    print('merging %d per-library dicts' % len(inputs), flush=True)
    for p in inputs:
        if not os.path.exists(p):
            sys.exit('FATAL: missing %s' % p)
        d = pickle.load(open(p, 'rb'))
        if not isinstance(d, dict):
            sys.exit('FATAL: %s does not hold a dict but a %s' % (p, type(d)))
        # Property 2: upstream drops an empty cnt, so an empty dict here is a
        # library that produced no assignments -- report it rather than hide it.
        if len(d) == 0:
            sys.exit('FATAL: %s is empty -- that library produced no gene '
                     'assignments at all' % p)
        for k, v in d.items():
            # Property 1: disjoint column keys.
            if k in gcnt:
                sys.exit('FATAL: column key %r appears in more than one input '
                         '(second occurrence in %s) -- merging would double count' % (k, p))
            gcnt[k] = v
            print('  %-40s %8d gene entries  (%s)'
                  % (k, len(v), os.path.basename(p)), flush=True)

    print('\ncolumns: %s' % sorted(gcnt), flush=True)

    pickle.dump(gcnt, open(outpfx + 'dict.pickle', 'wb'))

    # Property 3: the union of row keys, exactly as upstream's DataFrame call does.
    cntdf = pd.DataFrame(gcnt)
    print('merged frame shape: %s' % (cntdf.shape,), flush=True)
    n_union = len(set().union(*[set(v) for v in gcnt.values()]))
    if cntdf.shape[0] != n_union:
        sys.exit('FATAL: frame has %d rows but the union of gene keys is %d'
                 % (cntdf.shape[0], n_union))
    if cntdf.shape[1] != len(gcnt):
        sys.exit('FATAL: frame has %d columns but %d were merged'
                 % (cntdf.shape[1], len(gcnt)))
    del gcnt

    # Upstream writes <output>.pickle then shells out to gzip. Written directly
    # as .pickle.gz here: identical bytes to step 7's reader
    # (pickle.load(gzip.open(...))), and it cannot hit the bare-gzip
    # refuses-to-clobber trap that cost a whole VASA run.
    with gzip.open(outpfx + '.pickle.gz', 'wb') as fh:
        pickle.dump(cntdf, fh)

    # Read it back through step 7's own call, so a corrupt or half-written file
    # is caught here rather than an hour into step 7.
    back = pickle.load(gzip.open(outpfx + '.pickle.gz', 'rb'))
    if back.shape != cntdf.shape:
        sys.exit('FATAL: round-trip shape mismatch %s vs %s' % (back.shape, cntdf.shape))
    print('round-trip via pickle.load(gzip.open(...)) OK: %s' % (back.shape,), flush=True)
    print('wrote %sdict.pickle and %s.pickle.gz' % (outpfx, outpfx), flush=True)


if __name__ == '__main__':
    main()
