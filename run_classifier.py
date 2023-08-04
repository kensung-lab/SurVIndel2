import pysam, argparse, os
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from collections import defaultdict
import joblib
import features

cmd_parser = argparse.ArgumentParser(description='Classify SVs using a built ML model.')
cmd_parser.add_argument('in_vcf', help='Input VCF file.')
cmd_parser.add_argument('out_vcf', help='Output VCF file.')
cmd_parser.add_argument('stats', help='Stats of the test VCF file.')
cmd_parser.add_argument('svtype', help='SV type to filter.', choices=['DEL', 'DUP', 'INS', 'ALL'])
cmd_parser.add_argument('model_dir', help='Directory containing the trained model.')
cmd_args = cmd_parser.parse_args()
    
def read_false_ids(file_path, tolerate_no_fps = False):
    if not os.path.exists(file_path) and tolerate_no_fps:
        return set()
    with open(file_path, 'r') as file:
        false_ids = set([line.strip() for line in file])
    return false_ids

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
        feature_values = features.Features.record_to_features(record, min_depth, avg_depth, max_depth)
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
