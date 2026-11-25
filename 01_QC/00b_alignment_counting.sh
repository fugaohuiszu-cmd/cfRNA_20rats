#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-.}"
RESULT_DIR="$WORKDIR/result"
TRIM_DIR="$RESULT_DIR/trimming"
MAP_DIR="$RESULT_DIR/mapping"
CALL_DIR="$RESULT_DIR/calling"
FINAL_DIR="$WORKDIR/01_input"
REFERENCE_DIR="$WORKDIR/Reference"
CUT_MERGED="$WORKDIR/cut_merged"
SAMPLE_LIST="$RESULT_DIR/sample_barcode_list.txt"
JOBS=20
BT2_THREADS=6
FC_THREADS=64

mkdir -p "$TRIM_DIR"/{decompress,trimmomatic_15nt,trimmomatic_19nt,trimmomatic_23nt}
mkdir -p "$MAP_DIR"/{rRNA,miRNA,piRNA,GRCr8_15nt,GRCr8_23nt}
mkdir -p "$CALL_DIR"/{miRNA,piRNA,mRNA,lncRNA,snRNA,snoRNA} "$FINAL_DIR"

process_sample() {
    set -euo pipefail
    local id="$1"
    pigz -dc "$CUT_MERGED/$id"_cut_R1.fastq.gz > "$TRIM_DIR/decompress/$id"_cut_R1.fastq
    trimmomatic SE "$TRIM_DIR/decompress/$id"_cut_R1.fastq \
        "$TRIM_DIR/trimmomatic_15nt/$id"_cut_trim_R1.fastq \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:15
    trimmomatic SE "$TRIM_DIR/trimmomatic_15nt/$id"_cut_trim_R1.fastq \
        "$TRIM_DIR/trimmomatic_19nt/$id"_cut_trim_R1.fastq \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:19
    trimmomatic SE "$TRIM_DIR/trimmomatic_19nt/$id"_cut_trim_R1.fastq \
        "$TRIM_DIR/trimmomatic_23nt/$id"_cut_trim_R1.fastq \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:23

    bowtie2 -p "$BT2_THREADS" -x "$REFERENCE_DIR/rRNA/Rat_rRNA" \
        -U "$TRIM_DIR/trimmomatic_15nt/$id"_cut_trim_R1.fastq \
        -S "$MAP_DIR/rRNA/$id"_cut_trim_map2rRNA.sam > "$MAP_DIR/rRNA/$id".log 2>&1
    bowtie2 -k 1 -p "$BT2_THREADS" -x "$REFERENCE_DIR/miRNA/rno_miRNAs" \
        -U "$TRIM_DIR/trimmomatic_19nt/$id"_cut_trim_R1.fastq \
        -S "$MAP_DIR/miRNA/$id"_map2miRNA.sam > "$MAP_DIR/miRNA/$id".log 2>&1
    awk '$6 ~ /^[0-9]*M$/ {print $3}' "$MAP_DIR/miRNA/$id"_map2miRNA.sam |
        sort | uniq -c | sort -nr > "$CALL_DIR/miRNA/$id"_miRNA.txt
    bowtie2 -k 1 -p "$BT2_THREADS" -x "$REFERENCE_DIR/piRNA/piR_rno" \
        -U "$TRIM_DIR/trimmomatic_23nt/$id"_cut_trim_R1.fastq \
        -S "$MAP_DIR/piRNA/$id"_map2piRNA.sam > "$MAP_DIR/piRNA/$id".log 2>&1
    awk '$6 ~ /^[0-9]*M$/ {print $3}' "$MAP_DIR/piRNA/$id"_map2piRNA.sam |
        sort | uniq -c | sort -nr > "$CALL_DIR/piRNA/$id"_piRNA.txt
    bowtie2 -p "$BT2_THREADS" -x "$REFERENCE_DIR/GRCr8/rat_GRCr8" \
        -U "$TRIM_DIR/trimmomatic_23nt/$id"_cut_trim_R1.fastq \
        -S "$MAP_DIR/GRCr8_23nt/$id"_map2GRCr8.sam > "$MAP_DIR/GRCr8_23nt/$id".log 2>&1
    bowtie2 -p "$BT2_THREADS" -x "$REFERENCE_DIR/GRCr8/rat_GRCr8" \
        -U "$TRIM_DIR/trimmomatic_15nt/$id"_cut_trim_R1.fastq \
        -S "$MAP_DIR/GRCr8_15nt/$id"_map2GRCr8.sam > "$MAP_DIR/GRCr8_15nt/$id".log 2>&1
}

export TRIM_DIR MAP_DIR CALL_DIR REFERENCE_DIR CUT_MERGED BT2_THREADS
export -f process_sample
parallel --halt soon,fail=1 -j "$JOBS" process_sample :::: "$SAMPLE_LIST"

set --
while IFS= read -r id; do
    set -- "$@" "$MAP_DIR/GRCr8_23nt/$id"_map2GRCr8.sam
done < "$SAMPLE_LIST"
featureCounts -T "$FC_THREADS" -O -s 1 -t transcript -g gene_id \
    -a "$REFERENCE_DIR/Gencode/Rattus_norvegicus.GRCr8.114.protein_coding.gtf" \
    -o "$CALL_DIR/mRNA/mRNA_S1_transcript.txt" "$@"
featureCounts -T "$FC_THREADS" -O -s 1 -t transcript -g gene_id \
    -a "$REFERENCE_DIR/Gencode/Rattus_norvegicus.GRCr8.114.lncRNA.gtf" \
    -o "$CALL_DIR/lncRNA/lncRNA_S1_transcript.txt" "$@"

set --
while IFS= read -r id; do
    set -- "$@" "$MAP_DIR/GRCr8_15nt/$id"_map2GRCr8.sam
done < "$SAMPLE_LIST"
featureCounts -T "$FC_THREADS" -O -s 1 -t transcript -g gene_id \
    -a "$REFERENCE_DIR/Gencode/Rattus_norvegicus.GRCr8.114.snRNA.gtf" \
    -o "$CALL_DIR/snRNA/snRNA_S1_transcript.txt" "$@"
featureCounts -T "$FC_THREADS" -O -s 1 -t transcript -g gene_id \
    -a "$REFERENCE_DIR/Gencode/Rattus_norvegicus.GRCr8.114.snoRNA.gtf" \
    -o "$CALL_DIR/snoRNA/snoRNA_S1_transcript.txt" "$@"

: > "$FINAL_DIR/Total_reads.txt"
: > "$FINAL_DIR/Trimmed_reads.txt"
: > "$FINAL_DIR/Unmapped_reads.txt"
: > "$FINAL_DIR/Mapped_rRNA_reads.txt"
while IFS= read -r id; do
    awk 'END {print NR / 4 "\t" FILENAME}' "$TRIM_DIR/decompress/$id"_cut_R1.fastq \
        >> "$FINAL_DIR/Total_reads.txt"
    awk 'END {print NR / 4 "\t" FILENAME}' "$TRIM_DIR/trimmomatic_15nt/$id"_cut_trim_R1.fastq \
        >> "$FINAL_DIR/Trimmed_reads.txt"
    unmapped="$(samtools view -c -f 4 "$MAP_DIR/GRCr8_15nt/$id"_map2GRCr8.sam)"
    rrna="$(samtools view -c -F 4 "$MAP_DIR/rRNA/$id"_cut_trim_map2rRNA.sam)"
    printf '%s\t%s\n' "$unmapped" "$id"_unmapped.sam >> "$FINAL_DIR/Unmapped_reads.txt"
    printf '%s\t%s\n' "$rrna" "$id"_mapped_rRNA.sam >> "$FINAL_DIR/Mapped_rRNA_reads.txt"
done < "$SAMPLE_LIST"
