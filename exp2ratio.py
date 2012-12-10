#!/usr/bin/env python

# Convert expression to log2 ratio. 
# The expression value devide
# @author A. Terada

__author__ = "Aika TERADA"

import sys, os, math
from optparse import OptionParser

NAME_COLUMN = 0 # gene name column number
SEPARATOR = '\t' # file separator
DEFAULT_CONTROL = "1" # default number of control expression

def log2( value ):
	if value == 1:
		return 0
	else:
		return math.log( value, 2 )

##
# Read expression file
# exp_file: The expression-file name
# control_column: List of control expression column
# target_column: Integer of target expression
# return: the list of tuple:
#    tuple[0]: gene name
#    tuple[1]: List of control expressions vlaues
#    tuple[2]: Float value of target expression
##
def readExpFile( exp_file, control_column, target_column ):
	exp_list = []
	try:
		f = open( exp_file, 'r' )
		line = None
		for line in f:
			s = line[:-1].split(SEPARATOR)
			control_list = []
			for i in control_column:
				control_list.append( float(s[i]) )
			exp_list.append( tuple( [s[NAME_COLUMN], control_list, float(s[target_column])]) )
		f.close()
		return exp_list
	except IOError, e:
		sys.stderr.write("Error in read %s\n" % exp_file)
		sys.exit()

##
# Derive target-expression/mean(control-expression)
# exp_list: the list of tuple (obtained by readExpFile)
#    tuple[0]: gene name
#    tuple[1]: List of control expressions vlaues
#    tuple[2]: Float value of target expression
# return: List of tuple
#    tuple[0]: gene name
#    tuple[1]: ratio
##
def changeRatio( exp_list, min_threshold ):
	ratio_list = []
	for t in exp_list:
		control_log_list = map(lambda x: log2(x), t[1]) # control exp log2
		target_log = log2(t[2]) # target exp log2
		# If exp >= min_threshold of controls and target expressions at least,
		# calculate ratio and save ratio_list
		if (len([x for x  in control_log_list if x >= min_threshold]) > 0) or (target_log >= min_threshold):
			ratio = target_log - sum(control_log_list)/len(control_log_list)
			ratio_list.append(tuple([t[0], ratio]))
	return ratio_list

def mergeID( ratio_list ):
	merge_list = []
	output_set = set()
	for t in ratio_list:
		if not t[0] in output_set:
			merge_list.append( t )
			output_set.add( t[0] )
	return merge_list

def output( ratio_list, output_file ):
	try:
		f = open( output_file, 'w' )
		for t in ratio_list:
			f.write("%s%s%s\n" % (t[0], SEPARATOR, t[1]))
		f.close()
	except IOError, e:
		sys.stderr.write("Error during output the result to %s" % output_file)
		sys.exit()
	
def run( exp_file, output_file, control_column, target_column ):
	exp_list = readExpFile( exp_file, control_column, target_column )
	ratio_list = changeRatio( exp_list, 4 )
	ratio_list = mergeID( ratio_list ) # reduce the duplicate gene names
	output( ratio_list, output_file )
	sys.stdout.write("%s genes are included.\n" % len(ratio_list))

if __name__ == "__main__":
	usage = "usage: %prog expression-file output-file"
	p = OptionParser(usage = usage)
	# Option to set control column
	p.add_option('-c', '--control', dest = "control_column", default = DEFAULT_CONTROL,
				 help = "Control expression column number. When there are several controls, the columns are delimited by ','. The default is 1.")
	p.add_option('-t', '--target', dest = "target_column",
				 help = "Analyzed expression column number.")
		
	opts, args = p.parse_args()
	
	# check arguments
	if len(args) < 2:
		sys.stderr.write("Error: input expression-file and output-file\n")
		sys.exit()

	# Convert control column number to integer.
	control_column = map(int, opts.control_column.split(','))
	# Convert target column number to integer.
	target_column = None
	try:
		target_column = int(opts.target_column)
	except TypeError:
		sys.stderr.write("Error: select column number using -t option.\n")
		sys.exit()
	except ValueError:
		sys.stderr.write("Error: select column number must be integer.\n")
	
	# Check the file exist.
	if not os.path.isfile(args[0]):
		sys.stderr.write("IOError: No such file: \'" + args[0] + "\'\n")
		sys.exit()
		
	run( args[0], args[1], control_column, target_column )
