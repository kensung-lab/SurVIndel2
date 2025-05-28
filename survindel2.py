#!/usr/bin/env python

from __future__ import print_function
import argparse, os, pysam, pyfaidx, gc, sys
import numpy as np
from random_pos_generator import RandomPositionGenerator
import timeit

MAX_READS = 1000
GEN_DIST_SIZE = 100000
MAX_ACCEPTABLE_IS = 20000

VERSION = "1.1.4"

cmd_parser = argparse.ArgumentParser(description='SurVIndel2, a CNV caller.')
cmd_parser.add_argument('bam_file', help='Input bam file.')
cmd_parser.add_argument('workdir', help='Working directory for Surveyor to use.')
cmd_parser.add_argument('reference', help='Reference genome in FASTA format.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to be used.')
cmd_parser.add_argument('--seed', type=int, default=0, help='Seed for random sampling of genomic positions.')
cmd_parser.add_argument('--min_sv_size', type=int, default=50, help='Min SV size.')
cmd_parser.add_argument('--min_clip_len', type=int, default=15, help='Min length for a clip to be used.')
cmd_parser.add_argument('--max_seq_error', type=float, default=0.04, help='Max sequencing error admissible on the platform used.')
cmd_parser.add_argument('--max_clipped_pos_dist', type=int, default=3, help='Max clipped position distance.')
cmd_parser.add_argument('--min_size_for_depth_filtering', type=int, default=1000, help='Minimum size for depth filtering.')
cmd_parser.add_argument('--samplename', default='', help='Name of the sample to be used in the VCF output.'
                                                                 'If not provided, the basename of the bam/cram file will be used,'
                                                                 'up until the first \'.\'')
cmd_parser.add_argument('--sampling-regions', help='File in BED format containing a list of regions to be used to estimate'
                                                   'statistics such as depth.')
cmd_parser.add_argument('--log', action='store_true', help='Activate in-depth logging (can be very large and cryptic).')
cmd_parser.add_argument('--match_score', type=int, default=1, help='Match score used by the aligner that produced tha BAM/CRAM file (TODO: auto-determine).')
cmd_parser.add_argument('--min-diff-hsr', type=int, default=3, help='Minimum number of differences with the reference \
                        (considered as number of insertions, deletions and mismatches) for a read to be considered a hidden split read.')
cmd_parser.add_argument('--version', action='version', version="SurVIndel2 v%s" % VERSION, help='Print version number.')
cmd_args = cmd_parser.parse_args()

def run_cmd(cmd, error_msg=None):
    start_time = timeit.default_timer()
    print("Executing:", cmd)
    return_code = os.system(cmd)
    if return_code != 0:
        if error_msg:
            print(error_msg)
        else:
            print("Error executing:", cmd)
            print("Return code:", return_code)
        exit(1)
    elapsed = timeit.default_timer() - start_time
    print(cmd, "was run in %.2f seconds" % elapsed)

# Create config file in workdir

SURVINDEL_PATH = os.path.dirname(os.path.realpath(__file__))

with open(cmd_args.workdir + "/full_cmd.txt", "w") as full_cmd_out:
    print(" ".join(sys.argv[:]), file=full_cmd_out)

config_file = open(cmd_args.workdir + "/config.txt", "w")
config_file.write("threads %d\n" % cmd_args.threads)
config_file.write("seed %d\n" % cmd_args.seed)
config_file.write("min_sv_size %s\n" % cmd_args.min_sv_size)
config_file.write("min_clip_len %s\n" % cmd_args.min_clip_len)
config_file.write("max_seq_error %s\n" % cmd_args.max_seq_error)
config_file.write("match_score %s\n" % cmd_args.match_score)
config_file.write("max_clipped_pos_dist %s\n" % cmd_args.max_clipped_pos_dist)
config_file.write("min_size_for_depth_filtering %s\n" % cmd_args.min_size_for_depth_filtering)
config_file.write("min_diff_hsr %s\n" % cmd_args.min_diff_hsr)
config_file.write("log %d\n" % cmd_args.log)
if cmd_args.sampling_regions:
    config_file.write("sampling_regions %s\n" % cmd_args.sampling_regions)
config_file.write("version %s\n" % VERSION)

with pysam.AlignmentFile(cmd_args.bam_file) as bam_file, open(cmd_args.workdir + "/contig_map", "w") as contigs_file:
    for r in bam_file.references:
        contigs_file.write("%s\n" % r)

# Find read length

read_len = 0
bam_file = pysam.AlignmentFile(cmd_args.bam_file, reference_filename=cmd_args.reference)
for i, read in enumerate(bam_file.fetch(until_eof=True)):
    if i > MAX_READS: break
    read_len = max(read_len, read.query_length)
config_file.write("read_len %d\n" % read_len)

# Generate general distribution of insert sizes

reference_fa = pyfaidx.Fasta(cmd_args.reference)
rand_pos_gen = RandomPositionGenerator(reference_fa, cmd_args.seed, cmd_args.sampling_regions)
random_positions = []
for i in range(1,1000001):
    if i % 100000 == 0: print(i, "random positions generated.")
    chr, pos = rand_pos_gen.next()
    random_positions.append((chr, pos))
rand_pos_gen = None

with open("%s/random_pos.txt" % cmd_args.workdir, "w") as random_pos_file:
    for random_pos in random_positions:
        random_pos_file.write("%s %d\n" % random_pos)

general_dist = []
rnd_i = 0
while len(general_dist) < GEN_DIST_SIZE:
    chr, pos = random_positions[rnd_i]
    rnd_i += 1

    if pos > len(reference_fa[chr])-10000:
        continue

    i = 0
    for read in bam_file.fetch(contig=chr, start=pos, end=pos+10000):
        if read.is_proper_pair and not read.is_secondary and not read.is_supplementary and \
                0 < read.template_length < MAX_ACCEPTABLE_IS:
            if i > 100:
                break
            i += 1
            general_dist.append(read.template_length)

reference_fa = None

mean_is = np.mean(general_dist)
stddev_is = np.std(general_dist)

general_dist = [x for x in general_dist if abs(x-mean_is) < 5*stddev_is]

mean_is = int(np.mean(general_dist))
lower_stddev_is = int(np.sqrt(np.mean([(mean_is-x)**2 for x in general_dist if x < mean_is])))
higher_stddev_is = int(np.sqrt(np.mean([(x-mean_is)**2 for x in general_dist if x > mean_is])))
general_dist = None

min_is, max_is = mean_is-3*lower_stddev_is, mean_is+3.5*higher_stddev_is

config_file.write("min_is %d\n" % min_is)
config_file.write("avg_is %d\n" % mean_is)
config_file.write("max_is %d\n" % max_is)
config_file.close()

workspace = cmd_args.workdir + "/workspace"
if not os.path.exists(workspace):
    os.makedirs(workspace)

gc.collect()

read_categorizer_cmd = SURVINDEL_PATH + "/reads_categorizer %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference)
run_cmd(read_categorizer_cmd)

if cmd_args.samplename:
    sample_name = cmd_args.samplename
else:
    sample_name = os.path.basename(cmd_args.bam_file).split(".")[0]

clip_consensus_builder_cmd = SURVINDEL_PATH + "/clip_consensus_builder %s %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference, sample_name)
run_cmd(clip_consensus_builder_cmd)

normalise_cmd = SURVINDEL_PATH + "/normalise %s/sr.vcf.gz %s/sr.norm.vcf.gz %s" % \
                (cmd_args.workdir, cmd_args.workdir, cmd_args.reference)
run_cmd(normalise_cmd)

merge_identical_calls_cmd = SURVINDEL_PATH + "/merge_identical_calls %s/sr.norm.vcf.gz %s/sr.norm.dedup.vcf.gz" % \
                (cmd_args.workdir, cmd_args.workdir)
run_cmd(merge_identical_calls_cmd)

dp_clusterer = SURVINDEL_PATH + "/dp_clusterer %s %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference, sample_name)
run_cmd(dp_clusterer)
