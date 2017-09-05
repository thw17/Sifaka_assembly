from __future__ import print_function
import argparse
import cyvcf2
import numpy as np


def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		description="This program filters VCF files based on user defined "
		"thresholds")

	parser.add_argument(
		"--vcf", required=True,
		help="Input VCF file.  Must be bgzipped and tabix-indexed.")

	parser.add_argument(
		"--output_vcf", required=True,
		help="Output VCF file.")

	parser.add_argument(
		"--variant_caller", type=str.lower,
		default="freebayes", choices=["freebayes"],
		help="Variant caller used.  Currently only supports freebayes.  Default "
		"is freebayes.")

	parser.add_argument(
		"--QUAL", type=int, default=0,
		help="Minimum phred-scaled *site* quality (QUAL column of vcf). "
		"Default is 0.")

	parser.add_argument(
		"--sample_depth", type=int, default=0,
		help="Minimum depth per sample required to retain site.  Default is 0.")

	parser.add_argument(
		"--var_type", type=str.upper, nargs="*", default="ALL",
		choices=["ALL", "SNP", "INDEL", "MNP", "COMPLEX", "INS", "DEL"],
		help="Type of variant to retain. Default is to retain all.  Choices "
		"include (can be either all uppcase or all lower case): all, snp, "
		"indel, mnp.  Multiple options can be selected, e.g. '--type snp indel'")

	parser.add_argument(
		"--min_samples", type=int, default=0,
		help="Minimum number of samples passing all filters required to retain "
		"sites. Integer. Default is 0.")

	parser.add_argument(
		"--min_support", type=int, default=0,
		help="Minimum number of reads supporting an allele. Integer. Default is 0.")

	parser.add_argument(
		"--genotype_quality", type=int, default=0,
		help="Minimum phred-scaled genotype quality per sample. "
		"Default is 0.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	vcf = cyvcf2.VCF(args.vcf)

	vcf.add_to_header(
		"##Filter_vcf_CMD=python Filter_vcf.py "
		"--vcf {} "
		"--output_vcf {} "
		"--variant_caller {} "
		"--min_samples {} "
		"--QUAL {} "
		"--sample_depth {} "
		"--min_support {} "
		"--genotype_quality {} "
		"--type {}".format(
			str(args.vcf),
			str(args.output_vcf),
			str(args.variant_caller),
			str(args.min_samples),
			str(args.QUAL),
			str(args.sample_depth),
			str(args.min_support),
			str(args.genotype_quality),
			str(args.var_type)))

	out_vcf = cyvcf2.Writer(args.output_vcf, vcf)

	for variant in vcf:
		if variant.QUAL < args.QUAL:
			continue
		if variant.INFO.get("DP") < args.min_samples * args.sample_depth:
			continue
		if args.var_type != "ALL":
			var_type = variant.INFO.get("type")
			if args.var_type == "INDEL":
				if var_type != "ins":
					if var_type != "del":
						continue
			else:
				if var_type != args.var_type.lower():
					continue
		dp = variant.format('DP')
		dp = dp[np.where(dp >= args.sample_depth)]
		if len(dp) < args.min_samples:
			continue
		gq = variant.format('GQ')
		gq = gq[np.where(gq >= args.genotype_quality)]
		if len(gq) < args.min_samples:
			continue
		gt = variant.gt_types
		gt_ref = variant.gt_ref_depths
		gt_alt = variant.gt_alt_depths
		if args.min_support > 0:
			passing = 0
			for idx, i in enumerate(gt):
				if i == 0:
					if gt_ref[idx] >= args.min_support:
						passing += 1
				elif i == 1:
					if gt_ref[idx] >= args.min_support and gt_alt[idx] >= args.min_support:
						passing += 1
				elif i == 2:
					if gt_alt[idx] >= args.min_support:
						passing += 1
			if passing < args.min_samples:
				continue
		out_vcf.write_record(variant)

	out_vcf.close()
	# vcf.close()

if __name__ == "__main__":
	main()
