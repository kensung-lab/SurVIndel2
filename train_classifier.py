import pysam, argparse, os
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from collections import defaultdict
import joblib
import features

cmd_parser = argparse.ArgumentParser(description='Train ML model.')
cmd_parser.add_argument('training_prefixes', help='Prefix of the training VCF and FP files.')
cmd_parser.add_argument('svtype', help='SV type to filter.', choices=['DEL', 'DUP', 'INS', 'ALL'])
cmd_parser.add_argument('outdir')
cmd_parser.add_argument('--n-trees', type=int, default=5000, help='Number of trees in the random forest.')
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

training_prefixes = cmd_args.training_prefixes.split(",")
training_data, training_labels, variant_ids = None, None, None
for training_prefix in training_prefixes:
    vcf_training_data, vcf_training_labels, vcf_variant_ids = parse_vcf(training_prefix + ".vcf.gz", training_prefix + ".stats",
                                                                        training_prefix + ".fps", tolerate_no_fps = False)
    if training_data is None:
        training_data = vcf_training_data
        training_labels = vcf_training_labels
    else:
        for source in vcf_training_data:
            training_data[source] = np.concatenate((training_data[source], vcf_training_data[source]))
            training_labels[source] = np.concatenate((training_labels[source], vcf_training_labels[source]))

for sources in training_data:
    classifier = RandomForestClassifier(n_estimators=cmd_args.n_trees, max_depth=15, n_jobs=-1)
    classifier.fit(training_data[sources], training_labels[sources])

    #save model to outdir
    model_fname = cmd_args.outdir + "/" + sources + ".model"
    joblib.dump(classifier, open(model_fname, 'wb'))
