#!/usr/bin/env python3

import argparse
import logging
import os
from collections import namedtuple, Counter

import regex

from utils.trie import Trie

logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument('pairs_tsv', type=argparse.FileType('r'), help='TSV file with pairing information from linktags')
parser.add_argument('readtag_dir', type=str, help="Folder with readtag sequences")
parser.add_argument('clusters_tsv', type=argparse.FileType('r'), help='TSV file with clusters from readtags')
parser.add_argument('output_base', type=str, help='Folder to write output to')
args = parser.parse_args()
try_fuzzy = False

SequenceAndQuality = namedtuple('SequenceAndQuality', ['seq', 'qual'])


def read_fastq(fq_file):
    reads = {}

    with open(fq_file, "r") as file:
        for line in file:
            header = line[1:].strip().split()[0]
            sequence = file.readline().strip()
            file.readline()
            quality = file.readline().strip()
            reads[header] = SequenceAndQuality(seq=sequence, qual=quality)

    return reads


def process_paired(r1, r2):
    return read_fastq(r1), read_fastq(r2)


# remove old output
if os.path.exists(args.output_base):
    logging.error("Output path exists. Refusing to overwrite")
    exit(1)

sample_name = os.path.basename(args.clusters_tsv.name).replace('.15N.clusters.tsv', '')
sample_name_r1 = os.path.join(args.readtag_dir, f"{sample_name}_R1_001.fastq")
sample_name_r2 = os.path.join(args.readtag_dir, f"{sample_name}_R2_001.fastq")

fastq_reads = process_paired(sample_name_r1, sample_name_r2)

# read in pairs
line: str
pairs = {}
clusters = {}
t = Trie()
for i, line in enumerate(args.pairs_tsv):
    if line.startswith('#'):
        continue
    fields = line.strip().split('\t')
    _, *sequences = fields

    for sequence in sequences:
        node = t.insert(sequence, i)
        if len(node.values) != 1:
            logging.warning(f"Ambiguous sequence {sequence}")

        v = pairs.get(sequence, [])
        v.append(i)
        pairs[sequence] = v

    clusters[i] = sequences

trivial_found = 0
trivial_found_comp = 0
non_trivial_found = 0
current_rt_cluster_id = None
current_lt_cluster = None

cluster_sizes = {}

with open('matches.tsv', 'w') as matches:
    for line in args.clusters_tsv:
        if line.startswith('#'):
            if line[1:].startswith('Cluster'):
                current_rt_cluster_id = int(line[9:])
                current_lt_cluster = None
            continue
        else:
            fields = line.strip().split('\t')
            ident, sequence, origin, cluster_id, cluster_sequence = fields
            if current_lt_cluster is None:
                logging.info(f"Search match of cluster {ident}")
                exact_search_target = cluster_sequence[:15]
                node = t.search(exact_search_target)
                if node:
                    trivial_found += 1
                    logging.debug(f"Trivial match: n={trivial_found}")
                    current_lt_cluster = node.values
                    matches.write(f"{exact_search_target}\t{cluster_sequence}\texact\n")
                else:
                    if try_fuzzy:
                        # try to identify with fuzzy regex matching
                        for pairs_seq in pairs.keys():
                            allowed_mm = round(len(pairs_seq) * 0.1) + abs(len(pairs_seq) - len(cluster_sequence))
                            match = regex.match(f"{pairs_seq}{{e<={allowed_mm}}}", cluster_sequence)
                            if match:
                                non_trivial_found += 1
                                logging.debug(f"Non trivial match: n={non_trivial_found}")
                                current_lt_cluster = pairs.get(pairs_seq)
                                matches.write(f"{pairs_seq}\t{cluster_sequence}\tfuzzy\n")
                                break
                if current_lt_cluster is None:
                    current_lt_cluster = False
                    continue
            elif current_lt_cluster is False:
                continue

            for lt_cluster in current_lt_cluster:
                lt_cluster_seqs = clusters[lt_cluster]
                cluster_name = '-'.join(lt_cluster_seqs)
                path = os.path.join(args.output_base, cluster_name)
                os.makedirs(path, exist_ok=True)
                with open(os.path.join(path, f"{sample_name}_R1_001.matching.fastq"), 'a') as r1, open(
                        os.path.join(path, f"{sample_name}_R2_001.matching.fastq"), 'a') as r2:

                    r1_entry = fastq_reads[0].get(ident[1:])
                    r2_entry = fastq_reads[1].get(ident[1:])

                    # write out clipped sequences
                    if origin == "R1":
                        r1.write(f"@{ident[1:]}\n{r1_entry.seq[15+len(sequence)+20:]}\n+\n{r1_entry.qual[15+len(sequence)+20:]}\n")
                        r2.write(f"@{ident[1:]}\n{r2_entry.seq}\n+\n{r2_entry.qual}\n")
                    else:
                        r1.write(f"@{ident[1:]}\n{r1_entry.seq}\n+\n{r1_entry.qual}\n")
                        r2.write(f"@{ident[1:]}\n{r2_entry.seq[15+len(sequence)+18:]}\n+\n{r2_entry.qual[15+len(sequence)+18:]}\n")

                    v = cluster_sizes.get(cluster_name, [0, 0])
                    if lt_cluster_seqs[0] == cluster_sequence[:15]:
                        v[0] += 1
                    elif lt_cluster_seqs[1] == cluster_sequence[:15]:
                        v[1] += 1
                    cluster_sizes[cluster_name] = v

            logging.debug(cluster_id)
            matches.flush()


with open(os.path.join(f"{sample_name}.stats"), 'w') as stats:
    stats.write(f"#cluster name\tmatches seqA\tmatches seqB\tmatches total\n")
    for key, value in cluster_sizes.items():
        stats.write(f"{key}\t{value[0]}\t{value[1]}\t{value[0]+value[1]}\n")

print(non_trivial_found)
print(trivial_found)
print(cluster_id)
