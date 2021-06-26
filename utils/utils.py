from collections import namedtuple
import logging
import os
from typing import Generator, Tuple, List, Dict

FilePair = namedtuple('FilePair', ['R1', 'R2'])


def get_paired_files_from_dir(folder: str) -> Generator[FilePair, None, None]:
    """
    Identify paired files in a given folder. Files must follow Illumina naming scheme.
    :param folder: path of the folder
    :return: generator with paired files of the given folder as FilePair
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


def string_match(str1: str, str2: str, allowed_mismatches: int) -> Tuple[bool, int]:
    """
    Function to test if two strings are equal within the bounds of the given number of allowed mismatches.
    Implementation is similar to hamming distance. If the allowed number of mismatches is exceeded, False is returned.
    Input strings have to be of the same length.

    :param str1: (String) Input 1
    :param str2: (String) Input 2
    :param allowed_mismatches: (int) Number of allowed mismatches
    :return: (bool, int) Truth value if the two strings match within the bounds of x allowed mismatches and number
    of encountered mismatches. Note this is only an approximation of the minimal distance if the strings do not match.
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


def primer_matches(query: str, primer_alternatives: List[str], allowed_mismatches: int) -> bool:
    """
    Check if a query string is matched by any of the alternative primers within the allowed mismatches.
    """
    match_found = False

    for seq in primer_alternatives:
        match, _ = string_match(query, seq, allowed_mismatches)
        if match:
            return True

    return match_found


def read_fastq(fq_file: str) -> Tuple[List[str], List[str], List[str]]:
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


def process_paired(r1: str, r2: str) -> Dict:
    """
    This function reads FASTQ entries for paired files. At the end it is ensured that from both files the same number
    of IDs was extracted
    :param r1: (string) absolute path to R1 file
    :param r2: (string) absolute path to R2 file
    :return: (dict) a dict holding ids and sequences and qualities for R1 and R2
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


def complement(seq: str) -> str:
    """
    Generate the complement sequence
    """
    s = []
    for ch in seq:
        if ch == "A":
            s.append("T")
        elif ch == "T":
            s.append("A")
        elif ch == "C":
            s.append("G")
        elif ch == "G":
            s.append("C")
    return ''.join(s)


def reverse_complement(seq: str) -> str:
    """
    Generate the reverse complement sequence
    """
    return complement(seq[::-1])
