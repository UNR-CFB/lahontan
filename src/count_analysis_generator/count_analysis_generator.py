#!/usr/bin/python

import csv
import os
from optparse import OptionParser

class GeneCounts(object):
	def __init__(self):
		self.header = [ 'GeneID', 'Length' ]
		self.sample_names = []
		self.gene_lengths = {}
		self.gene_counts = {}

	def add_counts(self, sample_name, geneID, length, counts):
		if sample_name not in self.sample_names:
			self.sample_names.append(sample_name)

		if geneID not in self.gene_lengths:
			self.gene_lengths[geneID] = length

		if geneID in self.gene_counts:
			self.gene_counts[geneID][sample_name] = counts
		else:
			self.gene_counts[geneID] = { sample_name: counts }

	def get_header(self):
		return self.header + self.sample_names

	def __iter__(self):
		for k, v in self.gene_counts.iteritems():
			yield [k, self.gene_lengths[k]] + [v[sample] for sample in self.sample_names]


def main():
	parser = OptionParser()
	parser.add_option('-f', '--in', dest='in_files', type='string', action='callback', callback=multi_file_callback, help='input files, separate multiple files with comma')
	parser.add_option('-o', '--out', dest='out_file', type='string', action='store', help='output file name')

	option, args = parser.parse_args()

	counts = GeneCounts()
	for file_name in option.in_files:
		# use filename as header
		sample = os.path.splitext(os.path.basename(file_name))[0]

		for count in read_gene_counts(file_name):
			counts.add_counts(sample, count['geneid'], count['length'], count['counts'])

	write_gene_counts(option.out_file, counts)

def read_gene_counts(filename):
	with open(filename) as f:
		# throw away first line, contains program invocation
		f.readline()
		# throw away second line, contains header
		f.readline()

		# read counts file
		reader = csv.reader(f, delimiter='\t')
		for row in reader:
			yield { 'geneid': row[0], 'length': int(row[5]), 'counts': int(row[6]) }

def write_gene_counts(filename, genes):
	with open(filename, 'wb') as f:
		writer = csv.writer(f, delimiter='\t')

		# write header
		writer.writerow(genes.get_header())

		# write data
		for data in genes:
			writer.writerow(data)

def multi_file_callback(option, opt, value, parser):
	setattr(parser.values, option.dest, value.split(','))

if __name__ == '__main__':
	main()