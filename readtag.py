#!/usr/bin/env python3

import argparse
from collections import Counter

import regex

from utils.utils import *

logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument('readtag_dir', type=str, help="Folder with readtag sequences")
parser.add_argument('output_dir', type=str, help="Folder to write output to")
args = parser.parse_args()

readtag_dir = os.path.abspath(args.readtag_dir)
output_dir = os.path.abspath(args.output_dir)

logging.info(f'RT Folder: {readtag_dir!r}')
logging.info(f'Output Folder: {output_dir!r}')

# create output dir if it does not exist
os.makedirs(output_dir, exist_ok=True)

# AP1A-15N-LoopFw
fw_pattern = regex.compile("((TCGGAATTTCCCAGC)([A,T,G,C,N]{13,17})(AGAGTTTGATCCTGGCTCAG)){e<=4}", flags=regex.BESTMATCH)
# AP2A-15N-LoopRv
rw_pattern = regex.compile("((ACACCGGTATCAACC)([A,T,G,C,N]{13,17})(GCTACCTTGTTACGACTT)){e<=4}", flags=regex.BESTMATCH)
ap1a_alt = Counter()
loopfw_alt = Counter()
ap2a_alt = Counter()
looprv_alt = Counter()

for file_pair in get_paired_files_from_dir(readtag_dir):
    contains_n = 0
    foundInR1 = 0
    foundInR2 = 0
    counts = Counter()
    sample_id = os.path.splitext(os.path.basename(file_pair.R1))[0].split('_R1_')[0]
    total = 0
    logging.debug(f"Extract 15N sequences for {sample_id}")
    output_file = os.path.join(args.output_dir, f"{sample_id}.15N.tsv")

    with open(output_file, 'w') as out_fh:
        out_fh.write("#id\t15N sequence\tfound in\n")

        fastq = process_paired(file_pair.R1, file_pair.R2)
        for identifier, seq1, qual1, seq2, qual2 in zip(fastq["ids"], fastq["r1"]["seqs"], fastq["r1"]["quals"],
                                                        fastq["r2"]["seqs"], fastq["r2"]["quals"]):
            total += 1
            identifier = identifier.split()[0]

            # see if we have AP1A-15N-LoopFw in R1
            match = regex.match(fw_pattern, seq1[:60])
            if match:
                counts[len(match.groups()[1])] += 1
                out_fh.write(f"{identifier}\t{match.groups()[2]}\tR1\n")
                ap1a_alt[match.groups()[1]] += 1
                loopfw_alt[match.groups()[3]] += 1

                foundInR1 += 1
                print(seq1[:60])
                print(match.groups())
            else:
                # check if we have AP2A-15N-LoopRv in R2
                match = regex.match(rw_pattern, seq2[:60])
                if match:
                    counts[len(match.groups()[1])] += 1
                    ap2a_alt[match.groups()[1]] += 1
                    looprv_alt[match.groups()[3]] += 1

                    out_fh.write(f"{identifier}\t{match.groups()[2]}\tR2\n")
                    foundInR2 += 1

    logging.info(
        f"Length distribution of 15N sequence: {counts}. If there are many matches at the edge cases consider allowing a broader range.")
    logging.info(f"Alternatives matches for AP1A: {ap1a_alt}")
    logging.info(f"Alternatives matches for LoopFw: {loopfw_alt}")
    logging.info(f"Alternatives matches for AP2A: {ap2a_alt}")
    logging.info(f"Alternatives matches for LoopRv: {looprv_alt}")

    logging.info(
        f"R1: Found structure in {100 * foundInR1 / total}% (total: {total}, found: {foundInR1}) of sequences for sample {sample_id}")
    logging.info(
        f"R2: Found structure in {100 * foundInR2 / total}% (total: {total}, found: {foundInR2}) of sequences for sample {sample_id}")
