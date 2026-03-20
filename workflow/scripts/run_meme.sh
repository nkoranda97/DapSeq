#!/usr/bin/env bash
# Usage: run_meme.sh <fasta> <outdir> <threads> <logfile> [label] [nmotifs] [minw] [maxw] [mod]
set -euo pipefail

FASTA="$1"
OUTDIR="$2"
THREADS="$3"
LOGFILE="$4"
LABEL="${5:-FASTA}"
NMOTIFS="${6:-1}"
MINW="${7:-4}"
MAXW="${8:-12}"
MOD="${9:-oops}"

SEQ_COUNT=$(grep -c '^>' "$FASTA" || true)

if [ "$SEQ_COUNT" -eq 0 ]; then
    mkdir -p "$OUTDIR"
    echo "# MEME skipped: no sequences in $LABEL" > "$OUTDIR/meme.txt"
    echo "MEME skipped: no sequences in $LABEL" >> "$LOGFILE"
    exit 0
fi

meme "$FASTA" \
  -nmotifs "$NMOTIFS" -minw "$MINW" -maxw "$MAXW" \
  -dna -mod "$MOD" -nostatus \
  -p "$THREADS" -oc "$OUTDIR" 2>>"$LOGFILE"
