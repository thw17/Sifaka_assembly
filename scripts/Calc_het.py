from __future__ import print_function
import argparse
import csv
import cyvcf2
import numpy as np


def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		description="This program calculates the ratio of heterozygous sites"
		" to all callable sites. It requires a bed file of callable sites"
		" and a vcf. Both are expected to be filtered in the same way - no"
		"filtering is done within this program. It also does not parse sex"
		"linked sequence.")

	parser.add_argument(
		"--vcf", nargs='*', required=True,
		help="Input VCF file.  Must be bgzipped and tabix-indexed.")

	parser.add_argument(
		"--callable_regions", required=True,
		help="REQUIRED. Bed file (no header) containing callable "
		"regions for population (to accurately calculate pi).")

	args = parser.parse_args()


def calc_callable_seq(callable_regions_bed):
	""" """
	seq = 0

	with open(callable_regions_bed, "r") as f:
		for i in csv.reader(f, delimiter="\t"):
			seq += int(i[2]) - int(i[1])

	return seq


def main():
	args = parse_args()
	vcf = cyvcf2.VCF(args.vcf)
