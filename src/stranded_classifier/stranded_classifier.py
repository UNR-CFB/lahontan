#!/usr/bin/python

import csv
from optparse import OptionParser

def main():
	parser = OptionParser()
	parser.add_option('-1', '--read1', dest='read1', type='string', action='store', help='BLAST results file name from read 1')
	parser.add_option('-2', '--read2', dest='read2', type='string', action='store', help='BLAST results file name from read 2')
	parser.add_option('-p', '--percthresh', dest='perc_thresh', type='float', action='store', help='percentage threshold to be considered stranded', default=0.85)
	parser.add_option('-d', '--percdelta', dest='perc_delta', type='float', action='store', help='maximum difference in read1 and read2 plus percentages to be considered stranded', default=0.05)

	options, args = parser.parse_args()
	print classify_stranded(options.perc_thresh, options.perc_delta, options.read1, options.read2)

def classify_stranded(perc_thresh, perc_delta, read1, read2):
	read1_plus = percent_plus(accum_directions(read1))
	read2_plus = percent_plus(accum_directions(read2))

	print 'read1 plus percent: ', read1_plus
	print 'read2 plus percent: ', read2_plus
	print 'combined: ', read1_plus + read2_plus

	if abs(1 - read1_plus - read2_plus) < perc_delta:
		# TODO: Check for FF and RR conditions. Where both read1 and read2 >=perc_thresh.
		if read1_plus >= perc_thresh or read2_plus >= perc_thresh:
			return 'True'

	return 'False'

def accum_directions(filename):
	with open(filename, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		directions = { 'minus': 0,
				'plus': 0
			}
		for row in reader:
			directions[row[12]] += 1

	return directions

def percent_plus(directions_dict):
	return directions_dict['plus'] / float(directions_dict['plus'] + directions_dict['minus'])

if __name__ == '__main__':
	main()
