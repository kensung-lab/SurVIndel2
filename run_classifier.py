import pysam, argparse, os
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from collections import defaultdict
import joblib

cmd_parser = argparse.ArgumentParser(description='Classify SVs using a built ML model.')
cmd_parser.add_argument('in_vcf', help='Input VCF file.')
cmd_parser.add_argument('out_vcf', help='Output VCF file.')
cmd_parser.add_argument('stats', help='Stats of the test VCF file.')
cmd_parser.add_argument('svtype', help='SV type to filter.', choices=['DEL', 'DUP', 'INS', 'ALL'])
cmd_parser.add_argument('model_dir', help='Directory containing the trained model.')
cmd_args = cmd_parser.parse_args()

source_to_int = { "2SR" : 0, "2HSR" : 1, "HSR-SR" : 2, "SR-HSR" : 3, "1SR_LC" : 4, "1SR_RC" : 5,
                 "1HSR_LC" : 6, "1HSR_RC" : 7, "DP" : 8 }

svtype_to_int = { "DEL" : 0, "DUP" : 1, "INS" : 2 }

feature_names = [
    'svtype',
    'source',
    'svlen',
    'median_depths[0]',
    'median_depths[1]',
    'median_depths[2]',
    'median_depths[3]',
    'median_depths_ratio[0]',
    'median_depths_ratio[1]',
    'median_depths_above_max[0]',
    'median_depths_above_max[1]',
    'median_depths_above_max[2]',
    'median_depths_above_max[3]',
    'median_depths_below_min[0]',
    'median_depths_below_min[1]',
    'cluster_depths_above_max[0]',
    'cluster_depths_above_max[1]',
    'disc_pairs',
    'disc_pairs_highmapq',
    'disc_pairs_highmapq_ratio',
    'disc_pairs_trusted',
    'disc_pairs_trusted_ratio',
    'conc_pairs',
    'disc_pairs_maxmapq',
    'ptn_ratio',
    'ks_pval',
    'disc_pair_surrounding[0]',
    'disc_pair_surrounding[1]',
    'max_size_diff',
    'lb_diff',
    'ub_diff',
    'clipped_reads[0]',
    'clipped_reads[1]',
    'max_mapq[0]',
    'max_mapq[1]',
    'full_to_split_junction_score_ratio',
    'full_to_split_junction_score_diff',
    'split2_to_split1_junction_score_ratio[0]',
    'split2_to_split1_junction_score_ratio[1]',
    'split_junction_size[0]',
    'split_junction_size[1]',
    'split_to_size_ratio[0]',
    'split_to_size_ratio[1]',
    'split2_to_split1_junction_score_diff_ratio[0]',
    'split2_to_split1_junction_score_diff_ratio[1]',
    'ext_1sr_reads[0]',
    'ext_1sr_reads[1]',
    'hq_ext_1sr_reads[0]',
    'hq_ext_1sr_reads[1]'
]

def get_value(info, key, default):
    if key in info:
        return info[key]
    else:
        return default
    
def read_false_ids(file_path, tolerate_no_fps = False):
    if not os.path.exists(file_path) and tolerate_no_fps:
        return set()
    with open(file_path, 'r') as file:
        false_ids = set([line.strip() for line in file])
    return false_ids

def record_to_features(record, min_depth, avg_depth, max_depth):
    info = record.info
    svtype = svtype_to_int[record.info['SVTYPE']]
    source = source_to_int[info['SOURCE']]
    svlen = abs(float(info['SVLEN']))
    median_depths = [float(x) for x in info['MEDIAN_DEPTHS']]
    median_depths_norm = [float(x)/avg_depth for x in info['MEDIAN_DEPTHS']]
    median_depths_ratio = [median_depths[0]/max(1, median_depths[1]), median_depths[3]/max(1, median_depths[2])]
    median_depths_above_max = [max(0, x-max_depth)/avg_depth for x in median_depths]
    median_depths_below_min = [max(0, min_depth-median_depths[0])/avg_depth, 
                               max(0, min_depth-median_depths[3])/avg_depth]
    if 'CLUSTER_DEPTHS' in info:
        cluster_depths_above_max = [max(0, float(x)-max_depth)/avg_depth for x in info['CLUSTER_DEPTHS']]
    else:
        cluster_depths_above_max = [0, 0]
    disc_pairs = float(info['DISC_PAIRS'])/avg_depth
    disc_pairs_highmapq = float(get_value(info, 'DISC_PAIRS_HIGHMAPQ', 0))/avg_depth
    disc_pairs_maxmapq_ratio = disc_pairs_highmapq/max(1, disc_pairs)
    disc_pairs_trusted = 0 #float(get_value(info, 'DISC_PAIRS_TRUSTED', 0))/avg_depth
    disc_pairs_trusted_ratio = disc_pairs_trusted/max(1, disc_pairs)
    disc_pairs_maxmapq = float(get_value(info, 'DISC_PAIRS_MAXMAPQ', 60))
    disc_pair_surrounding = [float(x)/avg_depth for x in info['DISC_PAIRS_SURROUNDING']]
    conc_pairs = float(get_value(info, 'CONC_PAIRS', 0))/avg_depth
    ptn_ratio = disc_pairs/max(1, disc_pairs+conc_pairs)
    ks_pval = 0
    if 'KS_PVAL' in info:
        ks_pval = max(0, float(info['KS_PVAL']))
    max_size_diff = 0
    if 'MAX_SIZE' in info:
        max_size = float(info['MAX_SIZE'])
        max_size_diff = max(0, svlen - 2*max_size)
    lb_diff, ub_diff = 0, 0
    if 'REMAP_LB' in info and 'REMAP_UB' in info:
        remap_lb = float(info['REMAP_LB'])
        remap_ub = float(info['REMAP_UB'])
        lb_diff = max(0, remap_lb-record.pos)
        ub_diff = max(0, record.stop-remap_ub)
    clipped_reads = [float(x)/avg_depth for x in info['CLIPPED_READS']]
    max_mapq = [float(x) for x in info['MAX_MAPQ']]
    full_junction_score = float(info['FULL_JUNCTION_SCORE'])
    split_junction_score1 = [float(x) for x in info['SPLIT_JUNCTION_SCORE']]
    split_junction_score2 = [float(x) for x in info['SPLIT_JUNCTION_SCORE2']]
    split_junction_size = [float(x) for x in info['SPLIT_JUNCTION_SIZE']]
    if sum(split_junction_score1) == 0:
        full_to_split_junction_score_ratio = 0
        full_to_split_junction_score_diff = 0
        split2_to_split1_junction_score_ratio = [0, 0]
        split2_to_split1_junction_score_diff_ratio = [0, 0]
        split_to_size_ratio = [0, 0]
    else:
        full_to_split_junction_score_ratio = full_junction_score/sum(split_junction_score1)
        full_to_split_junction_score_diff = full_junction_score-sum(split_junction_score1)
        split2_to_split1_junction_score_ratio = [split_junction_score2[0]/split_junction_score1[0], split_junction_score2[1]/split_junction_score1[1]]
        s1 = max(split_junction_size[0]-split_junction_score1[0], 0.01)/max(split_junction_size[0]-split_junction_score2[0], 0.01)
        s2 = max(split_junction_size[1]-split_junction_score1[1], 0.01)/max(split_junction_size[1]-split_junction_score2[1], 0.01)
        split2_to_split1_junction_score_diff_ratio = [s1, s2]
        split_to_size_ratio = [split_junction_score1[0]/split_junction_size[0], split_junction_score1[1]/split_junction_size[1]]
    if 'EXT_1SR_READS' in info:
        ext_1sr_reads = [float(x)/avg_depth for x in info['EXT_1SR_READS']]
        hq_ext_1sr_reads = [float(x)/avg_depth for x in info['HQ_EXT_1SR_READS']]
    else:
        ext_1sr_reads = [0, 0]
        hq_ext_1sr_reads = [0, 0]

    feature_values = [svtype, source, svlen] + median_depths_norm + median_depths_ratio
    feature_values += median_depths_above_max + median_depths_below_min + cluster_depths_above_max
    feature_values += [disc_pairs, disc_pairs_highmapq, disc_pairs_maxmapq_ratio, disc_pairs_trusted, disc_pairs_trusted_ratio]
    feature_values += [conc_pairs, disc_pairs_maxmapq, ptn_ratio, ks_pval] + disc_pair_surrounding + [max_size_diff] + [lb_diff, ub_diff]
    feature_values += clipped_reads + max_mapq + [full_to_split_junction_score_ratio, full_to_split_junction_score_diff]
    feature_values += split2_to_split1_junction_score_ratio + split_junction_size + split_to_size_ratio + split2_to_split1_junction_score_diff_ratio
    feature_values += ext_1sr_reads + hq_ext_1sr_reads
    return feature_values

# Function to parse the VCF file and extract relevant features using pysam
def parse_vcf(vcf_fname, stats_fname, fp_fname, tolerate_no_fps = False):
    false_ids = read_false_ids(fp_fname, tolerate_no_fps=tolerate_no_fps)
    vcf_reader = pysam.VariantFile(vcf_fname)
    stats_reader = open(stats_fname, 'r')
    data_by_source = defaultdict(list)
    labels_by_source = defaultdict(list)
    variant_ids_by_source = defaultdict(list)

    # read the stats file and extract the relevant values
    stats = dict()
    for line in stats_reader:
        sl = line.strip().split()
        stats[sl[0]] = int(sl[1])
    min_depth = stats['min_depth']
    avg_depth = stats['avg_depth']
    max_depth = stats['max_depth']

    for record in vcf_reader.fetch():
        if cmd_args.svtype != 'ALL' and record.info['SVTYPE'] != cmd_args.svtype:
            continue
        source = record.info['SVTYPE'] + "_" + record.info['SOURCE']
        feature_values = record_to_features(record, min_depth, avg_depth, max_depth)
        data_by_source[source].append(feature_values)
        labels_by_source[source].append(0 if record.id in false_ids else 1)
        variant_ids_by_source[source].append(record.id)
    for source in data_by_source:
        data_by_source[source] = np.array(data_by_source[source])
        labels_by_source[source] = np.array(labels_by_source[source])
    return data_by_source, labels_by_source, variant_ids_by_source

test_data, test_labels, test_variant_ids = parse_vcf(cmd_args.in_vcf, cmd_args.stats, "XXX", tolerate_no_fps = True)

def write_vcf(vcf_reader, vcf_header, pred_variant_ids, fname):
    pred_variant_ids = set(pred_variant_ids)
    vcf_writer = pysam.VariantFile(fname, 'w', header=vcf_header)
    for record in vcf_reader.fetch():
        if record.id in pred_variant_ids:
            record.info['HARD_FILTERS'] = ",".join(record.filter.keys())
            record.filter.clear()
            record.filter.add('PASS')
            vcf_writer.write(record)
    vcf_writer.close()

pred_labels = list()
pred_variant_ids = list()
for source in test_data:
    classifier = joblib.load(cmd_args.model_dir + '/' + source + '.model')
    predictions = classifier.predict(test_data[source])
    pred_labels.append(predictions)
    pred_variant_ids.append([test_variant_ids[source][i] for i in range(len(predictions)) if predictions[i] == 1])
pred_labels = np.concatenate(pred_labels)
pred_variant_ids = np.concatenate(pred_variant_ids)

# write the predictions to a VCF file
vcf_reader = pysam.VariantFile(cmd_args.in_vcf)
header = vcf_reader.header
header.add_line('##INFO=<ID=HARD_FILTERS,Number=.,Type=String,Description="PASS or not according to hard filters.">')
write_vcf(vcf_reader, header, pred_variant_ids, cmd_args.out_vcf)
