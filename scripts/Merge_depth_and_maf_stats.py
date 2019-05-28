import pandas as pd
import argparse

def parse_args():

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--input_depth", required=True,
		help="REQUIRED. Input per target BED file with depth stats.")

	parser.add_argument(
		"--input_maf_stats", required=True,
		help="REQUIRED. Input per target BED file with maf stats.")

	parser.add_argument(
		"--output", required=True,
		help="REQUIRED. Output BED file.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	input_depth = pd.read_csv(args.input_depth, sep="\t")
	input_maf = pd.read_csv(args.input_maf_stats, sep="\t")

	input_depth = input_depth.rename(columns={'length': 'sites_covered'})

	df = pd.merge(
		input_maf, input_depth[['sites_covered', 'mean_depth', 'chrom_start', 'id']], on='chrom_start')

	df.to_csv(args.output, sep="\t")


if __name__ == "__main__":
		main()
