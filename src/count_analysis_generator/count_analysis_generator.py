#!/usr/bin/python

import csv
import os
from optparse import OptionParser

def main():
	parser = OptionParser()
	parser.add_option('-f', '--in', dest='in_files', type='string', action='callback', callback=multi_file_callback, help='input files, separate multiple files with comma')
	parser.add_option('-o', '--out', dest='out_file', type='string', action='store', help='output file name')

	option, args = parser.parse_args()

	aggregate_gene_counts = {}
	for file_name in option.in_files:
		# use filename as header
		header = os.path.splitext(os.path.basename(file_name))[0]
		gene_counts = get_gene_counts(file_name)
		for k, v in gene_counts.iteritems():
			if k in aggregate_gene_counts:
				aggregate_gene_counts[k].gene_counts[header] = v
			else:
				summary = GeneSummary()
				summary.length = v['length']
				summary.gene_counts = { header: v['counts'] }
				aggregate_gene_counts[k] = summary

	print aggregate_gene_counts

def get_gene_counts(filename):
	gene_counts = {}
	with open(filename) as f:
		# throw away first line, contains program invocation
		f.readline()

		# read counts file
		reader = csv.reader(f, delimiter='\t')
		for row in reader:
			gene_counts[row[0]] = { 'length': row[5],
				'counts': row[6] }

	return gene_counts

def multi_file_callback(option, opt, value, parser):
	setattr(parser.values, option.dest, value.split(','))


class GeneCounts(object):
	def __init__(self):
		self.header = [ 'GeneID', 'Length' ]
		self.sample_names = []
		self.gene_length = {}
		self.gene_counts = {}

	def add_counts(sample_name, geneID, length, counts):
		if sample_name not in self.sample_names:
			self.sample_names.append(sample_name)

		if gene not in self.gene_length:
			self.gene_length[gene] = length

		self.sample_counts.setdefault(gene, { sample_name: counts }).append(counts)

	def __str__(self):
		pass

class GeneSummary(object):
	def __init__(self):
		self.length = 0
		self.gene_counts = {}


if __name__ == '__main__':
	main()