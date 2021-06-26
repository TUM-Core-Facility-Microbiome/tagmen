#!/usr/bin/env python3

import argparse
import logging
import os
from collections import namedtuple, Counter

import regex

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

FilePair = namedtuple('FilePair', ['R1', 'R2'])


def get_paired_files_from_dir(folder: str):
    """
    Identify paired files in a given folder. Files must follow Illumina naming scheme.
    :param folder:
    :return:
    """
    for file in os.listdir(folder):
        abspath = os.path.join(folder, file)
        if os.path.isfile(abspath):
            if os.path.splitext(abspath)[1] in ['.fastq']:
                if len(abspath.split("_R1_")) == 2:  # R1 file
                    # for a R1 file test if there is a matching R2 file in the same dir
                    paired = abspath.replace('_R1_', '_R2_')
                    if os.path.isfile(paired):
                        yield FilePair(R1=abspath, R2=paired)
                    else:
                        logging.warning(f"Could not identify paired file for {abspath}; "
                                        f"expected it to be named {paired}")


def string_match(str1, str2, allowed_mismatches):
    """
    Function to test if two strings are equal within the bounds of the given number of allowed mismatches.
    Implementation is similar to hamming distance. If the allowed number of mismatches is exceeded, False is returned.
    Input strings have to be of the same length.

    :param str1: (String) Input 1
    :param str2: (String) Input 2
    :param allowed_mismatches: (int) Number of allowed mismatches
    :return: (bool) Truth values if the two strings match within the bounds of x allowed mismatches
    """
    if not len(str1) == len(str2):
        logging.fatal("Could not get string similarity as they differ in length.")
        exit(1)

    mismatch_counter = 0

    i = 0
    for char in str1:
        if mismatch_counter > allowed_mismatches:
            return False, mismatch_counter

        if char is not str2[i]:
            mismatch_counter += 1
        i += 1

    if not mismatch_counter > allowed_mismatches:
        return True, mismatch_counter
    else:
        return False, mismatch_counter


def primer_matches(query, primer_alternatives, allowed_mismatches):
    match_found = False

    # for each possible primer sequence match against the region and return a boolean for match and hamming distance
    for seq in primer_alternatives:
        match, _ = string_match(query, seq, allowed_mismatches)
        if match:
            return True

    return match_found


def read_fastq(fq_file):
    """
    Function for extracting IDs, sequences and qualities from a FASTQ file
    :param fq_file: (string) absolute path to a FASTQ file
    :return: (list, list, list) 3 lists holding IDs, sequences and qualities are returned.
    """
    ids = []
    seqs = []
    quals = []

    with open(fq_file, "r") as file:
        i = 0
        for line in file:
            if i == 0:
                # This should be the ID in a FASTQ file. Exit if not.
                # print("ID:  {}".format(line))
                # ID lines in FASTQ files are tagged with an @ char
                if not line[0] == "@":
                    exit(1)
                ids.append(line[:-1])  # substring as we have to remove newline character

            if i == 1:
                # sequence line in FASTQ
                # print("SEQ: {}".format(line))
                seqs.append(line[:-1])  # substring as we have to remove newline character

            if i == 2:
                # the delimiter line in FASTQ, only yields a plus sign
                if not line[:-1] == "+":
                    exit(1)

            if i == 3:
                # Phred qualsity score line in FASTQ
                quals.append(line[:-1])  # substring as we have to remove newline character

            if not i == 3:
                i += 1
            else:
                i = 0

    return ids, seqs, quals


def process_paired(r1, r2):
    """
    This function reads FASTQ entries for paired files. At the end it is ensured that from both files the same number
    of IDs was extracted
    :param r1: (string) absolute path to R1 file
    :param r2: (string) absolute path to R2 file
    :return: (dict) a dict holding ids and for sequences and qualities for R1 and R2
    """
    r1_ids, r1_seqs, r1_quals = read_fastq(r1)
    r2_ids, r2_seqs, r2_quals = read_fastq(r2)

    if len(r1_ids) != len(r2_ids):
        exit(1)

    ret = {
        "ids": r1_ids,
        "r1": {
            "seqs": r1_seqs,
            "quals": r1_quals
        },
        "r2": {
            "seqs": r2_seqs,
            "quals": r2_quals
        }
    }
    return ret


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
