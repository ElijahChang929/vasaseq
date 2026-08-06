#!/bin/bash
set -eu
OWN=/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/own_version
W=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/out75
PY=/nemo/lab/turnerj/working/guangxin/envs/vasa/bin/python3
TMP=$(mktemp -d); trap 'rm -rf $TMP' EXIT
export TRIM_MINLEN=15
CELLS=$(ls $W/cells/*_cbc.fastq.gz | sed 's/.*_\([0-9]*\)_cbc.fastq.gz/\1/')

mkdir -p $TMP/cls $TMP/att
echo "$CELLS" | xargs -P 16 -I{} $PY $OWN/diagnostics/classify_reads.py --cell {} $W/cells $TMP/cls/{}.tsv
head -1 $TMP/cls/001.tsv > $W/read_classes.tsv
for f in $TMP/cls/*.tsv; do tail -n +2 $f; done | sort -k1,1 -k2,2 >> $W/read_classes.tsv
echo "read_classes.tsv done"

echo "$CELLS" | xargs -P 16 -I{} $PY $OWN/diagnostics/attribute_pass2_loss.py --mode vasa --cell {} $W/cells $TMP/att/{}.tsv
head -1 $TMP/att/001.tsv > $W/pass2_adapter_attribution.tsv
for f in $TMP/att/*.tsv; do tail -n +2 $f; done | sort -k1,1 -k2,2 >> $W/pass2_adapter_attribution.tsv
echo "pass2_adapter_attribution.tsv done"
