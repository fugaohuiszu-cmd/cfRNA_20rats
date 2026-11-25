#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-.}"
RAW_RUN1="$WORKDIR/rawdata"
RAW_RUN2="$WORKDIR/backup"
TRIM_DIR="$WORKDIR/result/trimming"
OUT_DIR="$WORKDIR/cut_merged"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MAPCSV="$SCRIPT_DIR/pairs.csv"
THREADS=16
CUTADAPT_CORES=2

mkdir -p "$TRIM_DIR"/{index,cutadapter,splitbar} "$OUT_DIR"

split_run() {
    local raw_dir="$1" column="$2" pool
    while IFS= read -r pool; do
        pigz -dc "$raw_dir/$pool"_R2.fastq.gz |
            awk 'NR % 4 == 2 {print substr($0, 1, 8)}' > "$TRIM_DIR/index/$pool"_R2_index.txt
        cutadapt --cores "$CUTADAPT_CORES" -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            -o "$TRIM_DIR/cutadapter/$pool"_cut_R1_temp.fastq.gz "$raw_dir/$pool"_R1.fastq.gz
        cutadapt --cores "$CUTADAPT_CORES" -a AAAAAAAAAA \
            -o "$TRIM_DIR/cutadapter/$pool"_cut_R1.fastq.gz "$TRIM_DIR/cutadapter/$pool"_cut_R1_temp.fastq.gz
        python3 "$SCRIPT_DIR/barcode_extract_new.py" \
            "$TRIM_DIR/cutadapter/$pool"_cut_R1.fastq.gz \
            "$TRIM_DIR/index/$pool"_R2_index.txt "$TRIM_DIR/splitbar"
    done < <(awk -F, -v column="$column" 'NR > 1 {sub(/_[ACGT]+$/, "", $column); print $column}' "$MAPCSV" | sort -u)
}

split_run "$RAW_RUN1" 1
split_run "$RAW_RUN2" 3

tail -n +2 "$MAPCSV" | tr -d '\r' |
    while IFS=, read -r sid1 sample1 sid2 sample2; do
        [[ "$sample1" == "$sample2" ]]
        barcode1="$(printf '%s\n' "$sid1" | awk -F_ '{print $NF}')"
        barcode2="$(printf '%s\n' "$sid2" | awk -F_ '{print $NF}')"
        [[ "$barcode1" == "$barcode2" ]]
        pigz -dc "$TRIM_DIR/splitbar/$sid1"_cut_R1.fastq.gz \
            "$TRIM_DIR/splitbar/$sid2"_cut_R1.fastq.gz |
            pigz -p "$THREADS" > "$OUT_DIR/$sample1"_"$barcode1"_cut_R1.fastq.gz
    done

awk -F, 'NR > 1 {n = split($1, a, "_"); print $2 "_" a[n]}' "$MAPCSV" \
    > "$WORKDIR/result/sample_barcode_list.txt"
