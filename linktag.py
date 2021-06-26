#!/usr/bin/env python3

import argparse

from utils.utils import *

logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument('linktag_dir', type=str, help="Folder with linktag sequences")
parser.add_argument('output_dir', type=str, help="Folder to write output to")
parser.add_argument('--allowed_mismatches', '-mm', type=int, default=3,
                    help="Allowed mismatches in RC_AP2A + AP1A sequence. Default=3")
args = parser.parse_args()

linktag_dir = os.path.abspath(args.linktag_dir)
output_dir = os.path.abspath(args.output_dir)

logging.info(f'LT Folder: {linktag_dir!r}')
logging.info(f'Output Folder: {output_dir!r}')

# create output dir if it does not exist
os.makedirs(output_dir, exist_ok=True)

FilePair = namedtuple('FilePair', ['R1', 'R2'])

RC_AP2A = "GGTTGATACCGGTGT"
AP1A = "TCGGAATTTCCCAGC"
search_for = RC_AP2A + AP1A


for file_pair in get_paired_files_from_dir(linktag_dir):
    found = 0
    total = 0
    sample_id = os.path.splitext(os.path.basename(file_pair.R1))[0].split('_R1_')[0]
    logging.debug(f"Extract 15N sequences for {sample_id}")
    output_file = os.path.join(output_dir, f"{sample_id}.15Npairs.tsv")
    with open(output_file, 'w') as out_fh:
        out_fh.write("#ID\tSEQUENCE_A\tSEQUENCE_B\n")

        fastq = process_paired(file_pair.R1, file_pair.R2)
        for identifier, seq1, qual1, seq2, qual2 in zip(fastq["ids"], fastq["r1"]["seqs"], fastq["r1"]["quals"],
                                                        fastq["r2"]["seqs"], fastq["r2"]["quals"]):
            identifier = identifier.split()[0]
            total += 1

            start_pos = 130
            length = 30
            for offset in range(16):
                search_range = seq1[offset + start_pos:offset + start_pos + length]
                if primer_matches(search_range, [search_for], args.allowed_mismatches):
                    found += 1
                    a_15n = seq1[offset + start_pos - 15:offset + start_pos]
                    b_15n = seq1[offset + start_pos + length:offset + start_pos + length + 15]
                    out_fh.write(f'{identifier}\t{reverse_complement(a_15n)}\t{b_15n}\n')
                    continue

    logging.info(
        f"Found structure in {100 * found / total}% (total: {total}, "
        f"found: {found}) of sequences for sample {sample_id}")
