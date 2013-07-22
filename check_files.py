#!/usr/bin/env python

# Check gene-IDs between target-file expression-file.

import sys, os
from optparse import OptionParser

__author__ = "Aika TERADA"

SEPARATOR_EXP = "\t"
SEPARATOR_CSV = ","

def readFile( input_file, separator ):
	try:
		f = open( input_file, 'r' )
		line = ""; genes_set = set()
		for line in f:
			if line.startswith("#"):
				continue
			gene = line.split(separator)[0]
			if gene in genes_set:
				sys.stderr.write("%s is contained two or more times in %s.\n" \
								 % (gene, input_file))
				sys.exit()
			genes_set.add( line.split(separator)[0] )
		f.close()
		return genes_set
	except IOError, e:
		sys.stderr.write("Error in read %s\n" % input_file)
		sys.exit()

def compareSet( set1, set2 ):
	diff_12 = set1 - set2
	diff_21 = set2 - set1
	# If set 1 and set2 are identical, finish this program
	if len( diff_12 ) == 0 and len( diff_21 ) == 0:
		sys.stdout.write("Both file's genes are identical.\n")
		sys.exit()
	# If set1 and set2 are different, output the different genes.
	if len( diff_12 ) > 0:
		sys.stdout.write("Only included expression-file:\n")
		for i in diff_12:
			sys.stdout.write("    %s\n" % i)
	if len( diff_21 ) > 0:
		sys.stdout.write("Only included target-file:\n")
		for i in diff_21:
			sys.stdout.write("    %s\n" % i)

def run( target_file, exp_file ):
	targeted_genes_set = readFile( target_file, SEPARATOR_CSV )
	exp_genes_set = readFile( exp_file, SEPARATOR_EXP )
	compareSet( targeted_genes_set, exp_genes_set )
	

if __name__ == "__main__":
	usage = "usage: %prog target-file expression-file"
	p = OptionParser(usage = usage)

	opts, args = p.parse_args()
	
	# Check arguments
	if ( len(args) < 2 ):
		sys.stderr.write("Error: input target-file and expression-file.\n")
		sys.exit()

	# check the file exist.
	if not os.path.isfile(args[0]):
		sys.stderr.write("IOError: No such file: \'" + args[0] + "\'\n")
		sys.exit()
	if not os.path.isfile(args[1]):
		sys.stderr.write("IOError: No such file: \'" + args[1] + "\'\n")
		sys.exit()

	run( args[0], args[1] )
