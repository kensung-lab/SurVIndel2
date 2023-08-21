import pysam, argparse
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

test_data, test_labels, test_variant_ids = features.parse_vcf(cmd_args.in_vcf, cmd_args.stats, "XXX", cmd_args.svtype, tolerate_no_fps = True)

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
