#!/usr/bin/env python

# Convert the expression to 1 or 0 according to the gene expression profile.
# If the gene profile >= threshold, the expression is 1, otherwise 0.
# @author A. Terada

__author__ = "Aika TERADA"

import sys, os, math
from optparse import OptionParser

NAME_COLUMN = 0 # gene name column number
EXPRESSION_COLUMN = 1 # gene expression profile number
SEPARATOR = '\t' # file separator


##
# Read expression file
# exp_file: The expression-file name
# return: the list of tuple:
#    tuple[0]: gene name
#    tuple[2]: expression profile
##
def readExpFile( exp_file ):
	exp_list = []
	try:
		f = open( exp_file, 'r' )
		line = None
		for line in f:
			s = line[:-1].split(SEPARATOR)
			exp_list.append( tuple( [s[NAME_COLUMN], float(s[EXPRESSION_COLUMN])]) )
		f.close()
		return exp_list
	except IOError, e:
		sys.stderr.write("Error in read %s\n" % exp_file)
		sys.exit()

##
# Classify the gene expression profile to 1 or 0.
# exp_list: the list of tuple (obtained by readExpFile)
# threshold: threshold to classify the expression level.
#    tuple[0]: gene name
#    tuple[1]: expression profile
# return: List of tuple
#    tuple[0]: gene name
#    tuple[1]: 1 or 0 according to gene expression profile
##
def deriveFlag( exp_list, threshold ):
	flag_list = []; flag1 = 0;
	for t in exp_list:
		# If exp >= threshold the expressions to 1, otherwise 0
		if ( t[1] >= threshold ):
			flag_list.append(tuple([t[0], 1]))
			flag1 = flag1 + 1
		else:
			flag_list.append(tuple([t[0], 0]))
	return flag_list, flag1

def output( flag_list, output_file ):
	try:
		f = open( output_file, 'w' )
		for t in flag_list:
			f.write("%s%s%s\n" % (t[0], SEPARATOR, t[1]))
		f.close()
	except IOError, e:
		sys.stderr.write("Error during output the result to %s" % output_file)
		sys.exit()
	
def run( exp_file, threshold, output_file ):
	exp_list = readExpFile( exp_file )
	flag_list, flag1 = deriveFlag( exp_list, threshold )
	output( flag_list, output_file )
	sys.stdout.write("%s genes are included.\n" % len(flag_list))
	sys.stdout.write("%s genes are classified to expression level 1.\n" % flag1)

if __name__ == "__main__":
	usage = "usage: %prog expression-file threshold output-file"
	p = OptionParser(usage = usage)
		
	opts, args = p.parse_args()
	
	# check arguments
	if len(args) < 3:
		sys.stderr.write("Error: input expression-file, threshold and output-file\n")
		sys.exit()

	threshold = None
	# Convert threshold float.
	try:
		threshold = float(args[1])
	except ValueError:
		sys.stderr.write("Error: threshold must be float.\n")
	
	# Check the file exist.
	if not os.path.isfile(args[0]):
		sys.stderr.write("IOError: No such file: \'" + args[0] + "\'\n")
		sys.exit()
		
	run( args[0], threshold, args[2] )
