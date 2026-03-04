#!/usr/bin/env bash
# Usage: run_meme.sh <fasta> <outdir> <threads> <logfile> [label]
set -euo pipefail

FASTA="$1"
OUTDIR="$2"
THREADS="$3"
LOGFILE="$4"
LABEL="${5:-FASTA}"

SEQ_COUNT=$(grep -c '^>' "$FASTA" || true)

if [ "$SEQ_COUNT" -eq 0 ]; then
    mkdir -p "$OUTDIR"
    echo "# MEME skipped: no sequences in $LABEL" > "$OUTDIR/meme.txt"
    echo "MEME skipped: no sequences in $LABEL" >> "$LOGFILE"
    exit 0
fi

meme "$FASTA" \
  -nmotifs 1 -minw 4 -maxw 12 \
  -dna -mod oops -nostatus \
  -p "$THREADS" -oc "$OUTDIR" 2>>"$LOGFILE"
