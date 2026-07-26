#!/bin/bash
###############################################################################
# summarise_trim.sh -- one table per cell from the bench_trim.sh outputs.
#
# Every variant starts from the same 300,000 reads, so the column that settles
# the argument is UNIQ = the ABSOLUTE number of uniquely mapped reads. A variant
# that discards more reads can always show a better mapping *rate*; it cannot
# fake a bigger UNIQ.
###############################################################################
# the scripts are versioned in the repo, the results are not -- go to the data
cd "${TRIMTEST_DIR:-/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest}" || exit 1
IN=300000

printf "%-5s %-4s %9s %7s %9s %6s %9s %6s %7s %7s\n" \
       cell var "kept" "kept%" "UNIQ" "uniq%" "MULTI" "mm%" "tooShrt" "avgLen"
for c in 011 007 001 016; do
  for v in ${VARS:-v1 v2 v3 v4 v5 v6 v7 v8}; do
    L="${v}/${c}_Log.final.out"
    [ -s "$L" ] || { printf "%-5s %-4s %9s\n" "$c" "$v" "-- no STAR log --"; continue; }
    kept=$(grep -m1 "Number of input reads" "$L" | awk '{print $NF}')
    uniq=$(grep -m1 "Uniquely mapped reads number" "$L" | awk '{print $NF}')
    mult=$(grep -m1 "Number of reads mapped to multiple loci" "$L" | awk '{print $NF}')
    shrt=$(grep -m1 "% of reads unmapped: too short" "$L" | awk '{print $NF}')
    alen=$(grep -m1 "Average mapped length" "$L" | awk '{print $NF}')
    awk -v c="$c" -v v="$v" -v i="$IN" -v k="$kept" -v u="$uniq" -v m="$mult" \
        -v s="$shrt" -v a="$alen" 'BEGIN{
      printf "%-5s %-4s %9d %6.1f%% %9d %5.1f%% %9d %5.1f%% %7s %7s\n",
             c,v,k,100*k/i,u,100*u/i,m,100*m/i,s,a }'
  done
  echo
done

cat <<'EOF'
variants
  v0  current step2: TrimGalore -> cutadapt 1.18, -m 15, -a XX{5} on A/C/G/T
  v1  identical settings, cutadapt 5.1              (isolates the version)
  v2  -m 20, homopolymer run required raised 6 -> 20
  v3  -m 20, models the real 3' construct: polyA + 12 nt (rc CBC+UMI) + rc(RA5),
      plus --poly-a and --nextseq-trim=20 for the 2-colour poly-G artefact
  v4  control: TrimGalore only, NO homopolymer pass

round 2 -- built up from v4 one modifier at a time
  v5  v4 + trim the real 3' construct: polyA + 12 nt (rc CBC+UMI) + rc(RA5)
  v6  v5 + --poly-a
  v7  v6 + --nextseq-trim=20   (2-colour poly-G, NovaSeq X)
  v8  v6 + -a polyG=G{20}
EOF
