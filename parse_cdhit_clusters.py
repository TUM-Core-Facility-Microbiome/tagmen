#!/usr/bin/env python3

import argparse
import logging
from collections import namedtuple, Counter

logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument('readtag_tsv', type=argparse.FileType('r'), help="TSV output from readtag script")
parser.add_argument('cd_hit_est_fasta', type=argparse.FileType('r'), help="FASTA output from CD-HIT-EST")
parser.add_argument('cd_hit_est_clstr', type=argparse.FileType('r'), help="CLSTR output from CD-HIT-EST")
args = parser.parse_args()

# read in clstr file
logging.info("Read clusters from clstr file")

clusters = {}
tmp_cluster = ""
tmp_ids = []

# read in tsv file
logging.info("Read in TSV file")

ReadtagInfo = namedtuple('ReadtagInfo', ['sequence', 'origin'])

tsv = {}
line: str
for line in args.readtag_tsv:
    if line.startswith('#'):
        continue
    fields = line.strip().split('\t')
    ident, sequence, origin = fields
    tsv[ident] = ReadtagInfo(sequence=sequence, origin=origin)
logging.debug(tsv)

# read in fasta
id_to_cluster_sequence = {}
logging.info("Read in FASTA")
for line in args.cd_hit_est_fasta:
    if line.startswith(';'):
        continue
    if line.startswith('>'):
        header = line[1:].strip()
        sequence = args.cd_hit_est_fasta.readline().strip()
        id_to_cluster_sequence[header] = sequence
logging.debug(id_to_cluster_sequence)

for line in args.cd_hit_est_clstr:
    if line.startswith('>'):
        if tmp_cluster != "":
            clusters[tmp_cluster] = tmp_ids
            tmp_ids = []
        tmp_cluster = line[1:].strip()
    else:

        sequence_id = line.split('>')[1].split('...')[0]
        tmp_ids.append(sequence_id)
clusters[tmp_cluster] = tmp_ids

logging.debug(clusters)

# generate output
output_tsv_file = args.readtag_tsv.name.replace('.tsv', '.clusters.tsv')
logging.info(f"Writing enriched TSV file to {output_tsv_file}")
cluster_stats_file = args.readtag_tsv.name.replace('.tsv', '.clusters.stats')
logging.info(f"Writing stats file to {cluster_stats_file}")

histogram_counts = Counter()

with open(output_tsv_file, 'w') as out, open(cluster_stats_file, 'w') as stats:
    out.write("#id\t15N sequence\tfound in\tcluster\tcluster_seq\n")
    stats.write("#cluster id\tsequence\tsize\n")
    for cluster_id, cluster_sequence_ids in clusters.items():
        out.write(f"#{cluster_id}\n")

        # find cluster seq
        cluster_seq = ""
        for id_in_cluster in cluster_sequence_ids:
            if not cluster_seq == "":
                break
            cluster_seq = id_to_cluster_sequence.get(id_in_cluster, '')
        if cluster_seq == "":
            logging.warning(f"No cluster sequence found for {cluster_id}")

        cluster_size = len(cluster_sequence_ids)
        stats.write(f"{cluster_id}\t{cluster_seq}\t{cluster_size}\n")
        histogram_counts[cluster_size] += 1

        for sequence_id in cluster_sequence_ids:
            tsv_info = tsv[sequence_id]
            out.write(f"{sequence_id}\t{tsv_info.sequence}\t{tsv_info.origin}\t{cluster_id}\t{cluster_seq}\n")

cluster_hist_file = args.readtag_tsv.name.replace('.tsv', '.clusters.histdata')
logging.info(f"Writing stats file to {cluster_hist_file}")
with open(cluster_hist_file, 'w') as hist_data:
    hist_data.write(f"#cluster size\tcount\n")
    for cls, count in histogram_counts.items():
        hist_data.write(f"{cls}\t{count}\n")

args.readtag_tsv.close()
args.cd_hit_est_fasta.close()
args.cd_hit_est_clstr.close()
