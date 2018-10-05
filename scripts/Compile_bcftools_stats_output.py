import argparse
import collections
import os
import pandas as pd


def parse_args():

	""" Parse command line arguments """
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--input_files", nargs="+",
		help="Space separated list of input files "
		"(output of 'bcftools stats <file.vcf.gz> | grep ^SN)'")

	parser.add_argument(
		"--output_file", required=True,
		help="Path to and name of output file")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	data_dict = collections.OrderedDict()
	data_dict['species'] = []
	data_dict['assembly'] = []
	data_dict['caller'] = []
	data_dict['region'] = []

	ex_counter = 0
	in_counter = 0
	for i in args.input_files:
		print(i)
		base = os.path.basename(i)
		base_split = base.split(".")
		data_dict['species'].append(base_split[0])
		data_dict['assembly'].append(base_split[1])
		data_dict['caller'].append(base_split[2])
		data_dict['region'].append(base_split[4])
		ex_counter += 1
		print(ex_counter)

		with open(i, "r") as infile:
			for line in infile:
				line_split = line.split(" of ")
				data_split = [x.strip() for x in line_split[1].split(":")]

				id = data_split[0].replace(" ", "_")
				if id in data_dict:
					data_dict[id].append(data_split[1])
				else:
					data_dict[id] = [data_split[1]]
			in_counter += 1
			print(in_counter)
	print(len(args.input_files))
	for i in data_dict:
		print(i, len(data_dict[i]))
	print(data_dict)
	df = pd.DataFrame.from_dict(data_dict)
	df.to_csv(args.output_file, index=False)


if __name__ == "__main__":
	main()
