#!/usr/bin/env python

# Compare gene-ID between expression-file and association-file
# and make new expression-file and association-file which are consisted from same gene ID.
# The gene IDs are based on the input expression-file.

import sys, os
from optparse import OptionParser

__author__ = "Aika TERADA"

SEPARATOR_EXP = "\t"
SEPARATOR_CSV = ","

##
# read expression-file and return the set of gene IDs in expression-file
# exp_file: The expression filename ( inputted filename )
##
def readEXPFile( exp_file ):
	gene_list = []
	try:
		f = open( exp_file, 'r' ); line = ""
		for line in f:
			s = line[:-1].split(SEPARATOR_EXP)
			gene_list.append(s[0] )
		f.close()
		return gene_list
	except IOError, e:
		sys.stderr.write( "Error in read expression-file.\n" )
		sys.exit()

def readCSVFile( csv_file ):
	association_dict = {}
	try:
		f = open( csv_file, 'r' )
		line = ""; column_line = "";
		for line in f:
			line = line[:-1]
			if line.startswith("#"):
				column_line = line
				continue
			gene = line.split( SEPARATOR_CSV )[0]
			association_dict[ gene ] = line
		f.close()
		return column_line, association_dict
	except IOError, e:
		sys.stderr.write( "Error in read csv-file.\n" )
		sys.exit()
		
##
# read CSV file and output a new file when the gene included gene_set
##		
def makeCSVFile( out_csv_file, gene_list, association_dict, column_line ):
	tf_size = len( column_line.split(SEPARATOR_CSV) ) - 1
	try:
		fo = open( out_csv_file, 'w' )
		fo.write("%s\n" % column_line)
		for gene in gene_list:
			if gene in association_dict:
				fo.write("%s\n" % association_dict[gene])
			else:
				fo.write("%s" % gene)
				for i in xrange(0, tf_size):
					fo.write(",0")
				fo.write("\n")
		fo.close()
	except IOError, e:
		sys.stderr.write( "Error in make a new csv file.\n" )
		sys.exit()

def run( exp_file, csv_file,  out_csv_file):
	gene_list = readEXPFile( exp_file )
	column_line, association_dict = readCSVFile( csv_file )
	makeCSVFile( out_csv_file, gene_list, association_dict, column_line )
	sys.stdout.write("# of TFs: %s, # of genes : %s\n" % (len(column_line.split(","))-1, len(association_dict) - 1))

if __name__ == "__main__":
	usage = "usage: %prog expression-file csv-file output-csv-file"
	p = OptionParser(usage = usage)
	opts, args = p.parse_args()

	# Check argments
	if (len(args) < 3):
		sys.stderr.write("Error: input expression-file, csv-file and ontput-csv-file.\n")
		sys.exit()
		
	# Check the file existence.
	if not os.path.isfile(args[0]):
		sys.stderr.write("IOError: No such file: \'" + args[0] + "\'\n")
		sys.exit()
	if not os.path.isfile(args[1]):
		sys.stderr.write("IOError: No such file: \'" + args[1] + "\'\n")
		sys.exit()

	run( args[0], args[1], args[2] )
