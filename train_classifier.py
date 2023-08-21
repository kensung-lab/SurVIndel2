import argparse
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

training_prefixes = cmd_args.training_prefixes.split(",")
training_data, training_labels, variant_ids = None, None, None
for training_prefix in training_prefixes:
    vcf_training_data, vcf_training_labels, vcf_variant_ids = features.parse_vcf(training_prefix + ".vcf.gz", training_prefix + ".stats",
                                                                        training_prefix + ".fps", cmd_args.svtype, tolerate_no_fps = False)
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
