import gzip
import sys
from contextlib import ExitStack
from pathlib import Path

read1_path, index_path, output_dir = sys.argv[1:]
pool = Path(read1_path).name.split("_")[0]
barcodes = [
    "AGCGTAGC", "TCGCCTTA", "CTAGTACG", "TTCTGCCT",
    "GCTCAGGA", "AGGAGTCC", "CATGCCTA", "GTAGAGAG",
    "CCTCTCTG", "CAGCCTCG", "GTAAGGAG", "TCCTCTAC",
    "ACTGCATA", "AAGGAGTA", "CTAAGCCT", "TCTCTCCG",
    "GCGTAAGA", "TAGATCGC", "CTCTCTAT", "AGAGTAGA",
    "TCATGAGC", "CTGCGCAT", "CCTAGAGT", "GAGCCTTA",
    "AGCAGATC", "GCGCCTCT", "TGCCTCTT", "GCAGCGTA",
    "GACGTCGA", "AAGGCTAT", "TTATTCGT", "AGGAAGTC",
    "AGGCTCCG", "TAGCTGCA", "TCGACTAG", "ACTCATAA",
    "TACTTACT", "GTAGCTCC", "GTCTTAGG", "GGCTCTGA",
    "AGGCGAAG", "TAGCGAGT", "CGCTCAGT", "CCTATCCT",
    "ATAGAGGC", "TACTACGC", "ACTGATCG", "CGTCTAAT",
    "TTCGGATG", "TAATCTTA", "CCTGAGAT", "GAGCGCTA",
    "TATCCTCT", "TTATGCGA", "TATAGCCT", "CGCGGCTA",
]

with ExitStack() as stack:
    reads = stack.enter_context(gzip.open(read1_path, "rt"))
    indexes = stack.enter_context(open(index_path))
    outputs = {
        barcode: stack.enter_context(
            gzip.open(Path(output_dir) / f"{pool}_{barcode}_cut_R1.fastq.gz", "wt")
        )
        for barcode in barcodes
    }
    for header in reads:
        sequence = reads.readline().rstrip("\n")
        separator = reads.readline()
        quality = reads.readline().rstrip("\n")
        index = indexes.readline().strip()
        if not separator or not index or len(sequence) != len(quality):
            raise ValueError("R1 records and R2 indexes do not match.")
        output = outputs.get(index[:8])
        if output is not None:
            read_id = header.split(" ")[0].split("/")[0].rstrip("\n")
            output.write(f"{read_id}\n{sequence}\n+\n{quality}\n")
    if indexes.readline():
        raise ValueError("R2 contains more indexes than R1 records.")
