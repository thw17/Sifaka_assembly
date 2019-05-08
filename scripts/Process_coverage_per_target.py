from __future__ import division
import argparse
import collections

def parse_args():

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--coverage_bed", required=True,
		help="REQUIRED. Input per-site BED file with coverage and target info.")

	parser.add_argument(
		"--depth_idx", type=int, default=3,
		help="Index (0-based Python index) of depth data. Default is 3.")

	parser.add_argument(
		"--target_idx", type=int, default=4,
		help="Index (0-based Python index) at which target data starts on each line. Default is 4.")

	parser.add_argument(
		"--output_distance", required=True,
		help="REQUIRED. Output per-site BED file to write distance and depth data")

	parser.add_argument(
		"--output_target_depth", required=True,
		help="REQUIRED. Output BED file to write target data")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	depth_idx = args.depth_idx
	target_idx = args.target_idx

	target_dict = collections.OrderedDict()
	current_target = None
	site_collector = []

	with open(args.coverage_bed, "r") as f:
		with open(args.output_distance, "w") as od:
			with open(args.output_target_depth, "w") as ot:
				for line in f:
					line1 = line.strip().split()
					if line1 == []:
						continue
					if line1[0][0] == "#":
						continue
					target = line1[target_idx:]
					target_tuple = (target[0], target[1], target[2])
					if target_tuple != current_target:
						# First iteration: set current_target and start site_collector
						# then continue
						if current_target is None:
							current_target = target_tuple
							site_collector = [line1[:target_idx]]
						# Done with previous target; process then reset current_target
						# and site_collector
						else:
							# process previous target
							start = int(current_target[1])
							stop = int(current_target[2])
							total_length = stop - start
							site_counter = 0
							total_depth = 0
							for site in site_collector:
								total_depth += site[depth_idx]
								site_counter += 1
							mean_depth = float(total_depth) / float(total_length)
							ot.write("{}\t{}\t{}\t{}\t{}\n".format(
								current_target[0], start, stop, site_counter, mean_depth))
							for site in site_collector:
								rel_depth = float(site[depth_idx]) / mean_depth
								od.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
									site[0], site[1], site[2], site[depth_idx],
									rel_depth, mean_depth, site[1] - start, stop - site[2]))

							# reset with new target
							current_target = target_tuple
							site_collector = [line1[:target_idx]]
					# Still in same target; add to site_collector then continue
					else:
						site_collector.append(line1[:target_idx])


if __name__ == "__main__":
	main()
