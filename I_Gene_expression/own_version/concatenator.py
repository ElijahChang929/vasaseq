#!/usr/bin/env python3
# Reads R1.fastq and R2.fastq files;
# selects reads with proper cell barcode;
# produces a new _cbc.fastq.gz file.
#
# BIG PICTURE: This script "concatenates" the barcode information (cell barcode
# + UMI, read from the barcode read) onto the NAME of the biological read, keeps
# only reads whose cell barcode is in the whitelist, and writes them out as new
# fastq(s). Downstream steps then recover cell + UMI directly from the read name.
#
# ###########################################################################
# OWN_VERSION FORK -- differs from I_Gene_expression/a_Mapping/concatenator.py
# ###########################################################################
# This dataset carries a fixed-length uninformative prefix at the 5' end of
# BOTH mates. That prefix is removed from the sequence AND its quality string
# before any other processing, so everything downstream (barcode/UMI parsing
# on the barcode read, and the emitted biological read) sees only real
# sequence.
#
# Controlled by --skip5 (both mates, default 21) with per-mate overrides
# --skip5r1 / --skip5r2. Pass --skip5 0 to get the upstream behaviour back.
#
# CONSEQUENCE FOR READ GEOMETRY: the barcode read must be long enough to still
# contain the barcode block after the prefix is removed, i.e.
#     len(barcode read) >= skip5 + lenumi + lencbc          (default 21+6+8 = 35)
# This is checked once, up front, against the first record of the first file
# pair, so a wrong --skip5 fails immediately instead of silently producing
# nonsense barcodes for an entire run.

import sys, os          # sys: exit()/argv ; os: filesystem checks + calling gzip/mkdir
import itertools as it  # used to enumerate barcode mismatch combinations/permutations
import argparse as argp # command-line argument parsing
import numpy as np      # only used here for np.mod() to track position within a 4-line fastq record
import gzip             # to read the gzipped input fastq files
import pandas as pd     # builds the barcode lookup table (DataFrame)
from pandas.io.parsers import read_csv  # reads the tab-separated barcode whitelist file
from collections import Counter         # counts how many whitelist barcodes each variant maps to
import glob             # finds the input fastq files by wildcard pattern

#### function to identify cells from barcodes, allowing some edit distances ####
def find_compatible_barcodes(barcode, HDmax = 0):
    """Given a barcode sequence and a maximum Hammin distance, it returns a list of compatible barcode sequences"""
    # Purpose: given one whitelist barcode, return every sequence that should be
    # treated as "the same" barcode when up to HDmax mismatches are tolerated.
    # If HDmax==0 we still allow a single 'N' substitution (sequencer no-calls);
    # otherwise we allow real base substitutions (C/T/G/A/N) up to HDmax positions.
    nt = ['N'] if HDmax == 0 else ['N','C','T','G','A']  # which characters may be substituted in
    HDmax = 1 if HDmax == 0 else HDmax                    # even at HDmax==0, allow 1 substitution slot (for N)

    compatible_barcodes = set([barcode])                 # start with the exact barcode itself
    for hd in range(1, HDmax+1):                          # for each number of simultaneous substitutions 1..HDmax
        comb = [''.join(l) for l in it.product(nt, repeat = hd)]  # every combination of replacement characters
        for c in comb:                                   # c = the actual replacement letters to insert
            for p in it.permutations(range(len(barcode)), hd):   # p = which positions in the barcode to replace
                s0 = barcode                             # copy of the original barcode to mutate
                for x, l in zip(p, c):                   # apply each (position x -> letter l) substitution
                    s0 = s0[:x] + l + s0[x+1:]           # rebuild the string with position x replaced by l
                compatible_barcodes.add(s0)              # record this mismatch variant
    return list(compatible_barcodes)                     # OUTCOME: list of all sequences equivalent to this barcode

#### check input variables ####
# Define and parse the command-line options. Each add_argument declares one flag,
# its help text, and (where relevant) a default value.
parser = argp.ArgumentParser(description = 'Concatenates bcread to bioread qname.')
parser.add_argument('--fqf', help = 'Fastq files names, without _Rx.fastq.gz')                              # PREFIX of the input fastqs
parser.add_argument('--bcread', '-bcr', help = 'read where to find the barcode (umi+cell)', choices = ['R1', 'R2'], default = 'R1')  # which mate holds the barcode
parser.add_argument('--bioread', '-bior', help = 'read where to find biological information', choices = ['R1', 'R2'], default = 'R2') # which mate holds the cDNA
parser.add_argument('--demux', '-dx', help = 'print different fastq file for each barcode', action = 'store_true')  # if set: one fastq per cell
parser.add_argument('--lencbc', '-lcbc', help = 'cell barcode length (integer)', type = int, default = 8)   # cell-barcode length in bases
parser.add_argument('--lenumi', '-lumi', help = 'umi length (integer)', type = int, default = 6)            # UMI length in bases
parser.add_argument('--umifirst', help = 'logical variable: umi before cel barcode', action = 'store_true') # if set: order is [UMI|CBC], else [CBC|UMI]
parser.add_argument('--cbcfile', '-cbcf', help = 'cell specific barcode file. Please, provide full name')   # path to barcode whitelist (tsv)
parser.add_argument('--cbchd', help = 'collapse cell barcodes with the given hamming distance', type = int, default = 0)  # mismatch tolerance
parser.add_argument('--outdir', help = 'output directory for cbc.fastq.gz and log files', type = str, default = './')     # where outputs go
# OWN_VERSION: length of the uninformative 5' prefix to strip from each mate
# before any parsing. --skip5 sets both; --skip5r1/--skip5r2 override one mate.
parser.add_argument('--skip5', help = 'bases to strip from the 5-prime end of BOTH reads before any processing (default 21)', type = int, default = 21)
parser.add_argument('--skip5r1', help = 'override --skip5 for R1 only', type = int, default = None)
parser.add_argument('--skip5r2', help = 'override --skip5 for R2 only', type = int, default = None)
args = parser.parse_args()   # OUTCOME: args holds all parsed values

# Copy each parsed argument into a short local variable for readability below.
fqr = args.fqf          # input fastq prefix
bcread = args.bcread    # 'R1' or 'R2' : where the barcode block is
bioread = args.bioread  # 'R1' or 'R2' : where the biological sequence is
lcbc = args.lencbc      # cell barcode length
lumi = args.lenumi      # UMI length
umifirst = args.umifirst  # bool: UMI comes before cell barcode
cbcfile = args.cbcfile  # whitelist file path
hd = args.cbchd         # Hamming distance tolerance
outdir = args.outdir    # output directory
demux = args.demux      # bool: split output per cell
# OWN_VERSION: resolve the per-mate 5' skip lengths (per-mate override wins).
skip5r1 = args.skip5 if args.skip5r1 is None else args.skip5r1  # bases to drop from R1
skip5r2 = args.skip5 if args.skip5r2 is None else args.skip5r2  # bases to drop from R2
if skip5r1 < 0 or skip5r2 < 0:
    sys.exit('--skip5/--skip5r1/--skip5r2 must be >= 0')
print('OWN_VERSION 5-prime skip: R1=%d nt, R2=%d nt' % (skip5r1, skip5r2))

#### Find input fastq files ####
# Expand the prefix into actual file lists via wildcards (handles multiple lanes).
fq1s = sorted(glob.glob(fqr + '*_R1*.fastq.gz'))  # all R1 files, sorted so R1/R2 line up
fq2s = sorted(glob.glob(fqr + '*_R2*.fastq.gz'))  # all R2 files, sorted to match
print(fq1s, fq2s)                                  # echo which files were found (log/debug)

if len(fq1s) != len(fq2s):                         # sanity check: R1 and R2 counts must match
    sys.exit("Please, different number of input and output fastq files")  # abort if mismatched

if len(fq1s) == len(fq2s) == 0:                    # nothing matched the prefix at all
    sys.exit('fastq files not found')              # abort

#### OWN_VERSION: up-front read-geometry check ####
# Peek at the first record of the first file pair and confirm that, after the
# 5' prefix is removed, the barcode read is still long enough to hold the
# barcode block. Doing this here means a wrong --skip5 (or a file pair that
# does not actually carry the prefix) fails in seconds, rather than silently
# parsing barcodes out of the wrong offset for hundreds of millions of reads.
def _peek_first_seq(path):
    """Return the sequence line of the first fastq record, as a str."""
    with gzip.open(path) as fh:
        fh.readline()                              # header
        return str(fh.readline().rstrip(), 'utf-8')  # sequence

_s1_probe = _peek_first_seq(fq1s[0])               # first R1 read
_s2_probe = _peek_first_seq(fq2s[0])               # first R2 read
_need = (args.lenumi + args.lencbc)                # bases the barcode block needs
_skip_bc  = skip5r1 if bcread == 'R1' else skip5r2   # skip applied to the barcode mate
_probe_bc = _s1_probe if bcread == 'R1' else _s2_probe
print('first-record lengths: R1=%d nt, R2=%d nt' % (len(_s1_probe), len(_s2_probe)))
if len(_probe_bc) < _skip_bc + _need:
    sys.exit(
        'barcode read (%s) is %d nt, too short to skip %d nt and still read a '
        '%d nt barcode block (umi %d + cbc %d); it needs >= %d nt. '
        'Check --skip5/--skip5r1/--skip5r2, --lenumi and --lencbc.'
        % (bcread, len(_probe_bc), _skip_bc, _need, args.lenumi, args.lencbc,
           _skip_bc + _need))
_skip_bio = skip5r1 if bioread == 'R1' else skip5r2  # skip applied to the biological mate
_probe_bio = _s1_probe if bioread == 'R1' else _s2_probe
if len(_probe_bio) <= _skip_bio:
    sys.exit(
        'biological read (%s) is %d nt, so skipping %d nt would leave nothing. '
        'Check --skip5/--skip5r1/--skip5r2.'
        % (bioread, len(_probe_bio), _skip_bio))

#### Read barcodes and expand set according to input hamming distance ####
if not os.path.isfile(cbcfile):                    # whitelist file must exist
    sys.exit("Barcode file not found")             # abort if missing

bc_df = read_csv(cbcfile, sep = '\t', names = ['bc','cellID'], index_col = 0)  # load whitelist: index=barcode seq, col=cellID
print(bc_df.head())                                # show first rows (log/debug)
bc_df['compatible_bcs'] = bc_df.apply(lambda x: find_compatible_barcodes(x.name, hd), axis = 1)  # per barcode, list its mismatch variants
cnt_allbcs = Counter([x for idx in bc_df.index for x in bc_df.loc[idx, 'compatible_bcs']])        # count how many whitelist entries claim each variant
# Build the final lookup table: keep a variant ONLY if it is unambiguous
# (maps to exactly one whitelist barcode, cnt==1). Rows: variant seq -> cellID + original barcode.
allbc_df = pd.DataFrame({x: {'cellID': bc_df.loc[idx,'cellID'], 'original': idx} for idx in bc_df.index for x in bc_df.loc[idx, 'compatible_bcs'] if cnt_allbcs[x]==1}).T
# OUTCOME: allbc_df is the fast lookup used per-read: given an observed cell
# barcode, return its cellID and the canonical whitelist barcode.

### Create output directory if it does not exist ####
if not os.path.isdir(outdir):                      # make the output folder if needed
    os.system('mkdir '+outdir)                     # (shell mkdir)

#### Read fastq files and assign cell barcode and UMI ####
# Open the output handle(s). Non-demux: one combined file. Demux: a dict of open
# file handles keyed by whitelist barcode, one file per cell (named by cellID).
if not demux:
    fout = open(outdir + '/' + fqr + '_cbc.fastq', 'w')  # single output fastq
else:
    fout = {idx:  open(outdir + '/' + fqr + '_' + str(bc_df.loc[idx, 'cellID']).zfill(3) + '_cbc.fastq', 'w') for idx in bc_df.index}  # one file per cell

ns = 0; nt = 0   # ns = reads with a valid barcode (kept); nt = total reads seen
nshort = 0       # OWN_VERSION: reads discarded for being shorter than the 5' prefix
for fq1, fq2 in zip(fq1s, fq2s):                   # iterate lane by lane, R1 paired with R2
    print(fq1, fq2)                                # log current file pair
    with gzip.open(fq1) as f1, gzip.open(fq2) as f2:   # open both gzipped fastqs
        for idx, (l1, l2) in enumerate(zip(f1, f2)):   # walk both files line-by-line in lockstep
            try:
                # Decode bytes->str and strip: keep only the first whitespace token
                # (e.g. drops the " 1:N:0:..." part of the fastq header line).
                l1, l2 = str(l1.rstrip().rsplit()[0], 'utf-8'), str(l2.rstrip().rsplit()[0], 'utf-8')
            except:
                print('check line '+str(idx))      # malformed/empty line -> note it
                l1 = ''; l2 = ''                   # blank out so logic below skips it
            l = np.mod(idx,4)                      # position within the 4-line fastq record (0=name,1=seq,2=+,3=qual)
            if l == 0:                             # -- header line --
                n1, n2 = l1, l2                    # store read names
                if not n1 == n2:                   # R1 and R2 names must match (same read)
                    print (n1, n2)                 # log the offending pair
                    sys.exit('fastq files not syncrhonized (@name)')  # abort: files out of sync
            if l == 1:                             # -- sequence line --
                s1, s2 = l1, l2                    # store R1/R2 sequences
            if l == 2:                             # -- '+' separator line --
                p1, p2 = l1[0], l2[0]              # first char of each (should be '+')
                if not p1 == p2:                   # both must match
                    print(l1, l2, p1, p2)          # log if not
                    sys.exit('fastq files not synchronized (+)')  # abort
            if l == 3 and len(l1) > 0:             # -- quality line: now we have a full record --
                q1, q2 = l1, l2                    # store R1/R2 quality strings
                nt += 1                            # count this complete read
                if len(q1) != len(s1) or len(q2) != len(s2):  # quality length must equal seq length
                    print('phred and read length do not match for '+n1)  # warn (does not abort)

                # --- OWN_VERSION: drop the uninformative 5' prefix from BOTH
                #     mates before anything else looks at them. Sequence and
                #     quality are sliced together so they stay the same length.
                #     After this point the rest of the script is unchanged: the
                #     barcode block is taken from the NEW start of the barcode
                #     read, and the emitted biological read is prefix-free.
                if skip5r1:
                    s1 = s1[skip5r1:]; q1 = q1[skip5r1:]
                if skip5r2:
                    s2 = s2[skip5r2:]; q2 = q2[skip5r2:]
                # A read that is shorter than its prefix leaves nothing usable;
                # count it and skip rather than emitting an empty/mangled record.
                if len(s1) < (lumi+lcbc if bcread == 'R1' else 1) or \
                   len(s2) < (lumi+lcbc if bcread == 'R2' else 1):
                    nshort += 1
                    continue

                # --- pull the barcode block off whichever read holds it, and
                #     trim it away so the remaining read is purely biological ---
                if bcread == 'R1':
                    bcseq = s1[:lumi+lcbc]         # first lumi+lcbc bases = barcode block (seq)
                    bcphred = q1[:lumi+lcbc]       # matching quality for the barcode block
                    s1 = s1[lumi+lcbc:]            # remainder of R1 = biological seq
                    q1 = q1[lumi+lcbc:]            # remainder of R1 quality
                elif bcread == 'R2':
                    bcseq = s2[:lumi+lcbc]         # same, if the barcode is on R2
                    bcphred = q2[:lumi+lcbc]
                    s2 = s2[lumi+lcbc:]
                    q2 = q2[lumi+lcbc:]
                # --- split the barcode block into cell barcode vs UMI, honoring order ---
                if not umifirst:                   # layout [cell-barcode | UMI]
                    cellbcseq = bcseq[:lcbc]       # first lcbc bases = cell barcode
                    umiseq = bcseq[lcbc:]          # rest = UMI
                    cellbcphred = bcphred[:lcbc]   # matching qualities
                    umiphred = bcphred[lcbc:]
                else:                              # layout [UMI | cell-barcode]  (VASA-plate)
                    cellbcseq = bcseq[lumi:]       # after the UMI = cell barcode
                    umiseq = bcseq[:lumi]          # first lumi bases = UMI
                    cellbcphred = bcphred[lumi:]
                    umiphred = bcphred[:lumi]

                try:
                    cellID, originalBC = allbc_df.loc[cellbcseq]  # look up observed barcode; KeyError if not in whitelist
                    ns += 1                        # matched -> count as a kept read
                    # Shift barcode/UMI quality chars up by 32 (ASCII) so they can't be
                    # mistaken for the biological read's phred string later on.
                    cellbcphred = ''.join([chr(ord(c)+32) for c in cellbcphred])
                    umiphred = ''.join([chr(ord(c)+32) for c in umiphred])

                    # Build the new read name: original name + ';'-joined tags.
                    # SS=observed CBC, CB=canonical whitelist CBC, QT=CBC qual,
                    # RX=UMI seq, RQ=UMI qual, SM=cellID (zero-padded to 3 digits).
                    name = ';'.join([n1] + [':'.join(x) for x in zip(['SS','CB','QT','RX','RQ','SM'], [cellbcseq, originalBC, cellbcphred, umiseq, umiphred, str(cellID).zfill(3)])])
                    s, q = (s1, q1) if bioread == 'R1' else (s2, q2)  # choose the biological seq+qual to emit
                    if not demux:
                        fout.write( '\n'.join([name, s, '+', q, '']))            # write one fastq record to the single file
                    else:
                        fout[cellbcseq].write('\n'.join([name, s, '+', q, '']))  # write to that cell's own file
                except:
                    continue                       # barcode not in whitelist -> drop this read, move on

#nt = (idx+1)/4
# Close all output file handle(s).
if not demux:
    fout.close()
else:
    for cbc in fout:
        fout[cbc].close()

#### LOG ####
# Write a summary log: the parameters used plus total reads and the fraction
# carrying a valid cell barcode (ns/nt) -- a key demultiplexing QC number.
fout = open(outdir + '/' + fqr + '.log', 'w')
fout.write('=> to generate cbc file <=\n')
fout.write(', '.join(['fastq file:', str(fqr),'\n']))
fout.write(', '.join(['full barcode in:', str(bcread),'\n']))
fout.write(', '.join(['biological read in:', str(bioread), '\n']))
fout.write(', '.join(['cell specific barcode length:', str(lcbc), '\n']))
fout.write(', '.join(['umi length:', str(lumi), '\n']))
fout.write(', '.join(['umi goes first:', str(umifirst),'\n']))
# OWN_VERSION: record the 5' prefix that was removed, so the log fully
# describes the read geometry this run assumed.
fout.write(', '.join(['5-prime skipped in R1:', str(skip5r1), '\n']))
fout.write(', '.join(['5-prime skipped in R2:', str(skip5r2), '\n']))
fout.write(', '.join(['total sequenced reads:', str(nt), '\n']))
fout.write(', '.join(['reads too short after 5-prime skip:', str(nshort), '\n']))
fout.write(', '.join(['reads with proper barcodes:', str(ns), str(1.0*ns/nt), '\n']))  # count + fraction kept
fout.close()

#### zip fastq file ####
# Compress the plain-text output(s) to .gz (the form downstream steps expect).
if not demux:
    os.system('gzip '+ outdir + '/' + fqr + '_cbc.fastq')      # single file -> _cbc.fastq.gz
else:
    os.system('gzip '+ outdir + '/' + fqr + '*_cbc.fastq')     # every per-cell file -> *_cbc.fastq.gz
# OUTCOME: outdir now holds <prefix>[_<cellID>]_cbc.fastq.gz plus <prefix>.log
