#!/usr/bin/env python3

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('assembly_stats_file', default='aggregated.assembly_stats.csv', help="File with assembly stats")
parser.add_argument('--threshold', type=int, default=1200, help="Minimal genome size")
parser.add_argument('--output-file', default='aggregated.fasta', help="Name for output file with all scaffolds")
parser.add_argument('--output-file-scaffold0', default='aggregated_scaffold_0.fasta',
                    help="Name for output file with only scaffold_0 entries")
args = parser.parse_args()
print(args)

stats_file = args.assembly_stats_file
output_all_file = args.output_file
output_scaffold0_file = args.output_file_scaffold0
threshold = args.threshold
total = 0
match = 0
with open(stats_file, 'r') as fh, open(output_all_file, 'w') as out_all, open(output_scaffold0_file, 'w') as out_sc0:
    for line in fh:
        if not line.startswith('.'):
            continue
        total += 1
        fields = line.split('\t')

        origin_folder = fields[0]
        genome_size = int(fields[4])

        if genome_size >= threshold:
            match += 1
            with open(os.path.join(origin_folder, 'a5.contigs.fasta'), 'r') as fasta:
                scaffold_counter = 0
                for fa_line in fasta:
                    if fa_line.startswith('>'):
                        fa_line = f">{origin_folder[2:]}_{fa_line[1:]}"  # modify header
                        scaffold_counter += 1
                    out_all.write(fa_line)
                    if scaffold_counter > 1:
                        out_sc0.write(fa_line)


print(f"Wrote all sequences to {output_all_file!r}")
print(f"Wrote scaffold_0 sequences to {output_scaffold0_file!r}")
print(f"From {total} entries {match} matched the threshold ({threshold}). {round(100 * match/total, 2)}%")
