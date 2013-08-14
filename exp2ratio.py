#!/usr/bin/env python

"""
Copyright (c) 2013, LAMP development team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the LAMP development team nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL LAMP DEVELOPMENT TEAM BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# Convert expression to log2 ratio. 
# The expression value devide
# @author A. Terada

__author__ = "Aika TERADA"

import sys, os, math
from optparse import OptionParser

EGF_COLS = range(3,30) # column numbers of the EGF-induced condition. 
HRG_COLS = range(30,58) # column numbers of the HRG-induced condition. 
CONTROL_COLS = (1, 2) # column numbers of the control condition.
MIN_THRESHOLD = 4 # If the gene expresses less than this value in all condition, the gene does not use for the analysis.
#DEFAULT_CONTROL = "1" # default number of control expression

NAME_COLUMN = 0 # gene name column number
SEPARATOR = '\t' # file separator

def log2( value ):
	if value == 1:
		return 0
	else:
		return math.log( value, 2 )

##
# Read expression file
# exp_file: The expression-file name
# control_column: Integer of control expression column
# target_column: Integer of target expression
# obj_cols: Integer list of experimented conditions
# return: the list of tuple:
#    tuple[0]: gene name
#    tuple[1]: Float value of control expressions
#    tuple[2]: Float value of target expression
##
def readExpFile( exp_file, control_column, target_column, obj_cols ):
	exp_list = []
	try:
		f = open( exp_file, 'r' )
		line = None
		for line in f:
			s = line[:-1].split(SEPARATOR)
			log_values = map(lambda x: log2(float(x)), s[1:]) # log2 expression values
			# If the gene expresses lower than the given threshold, then continue.
			if len([x for x in obj_cols if log_values[x-1] > MIN_THRESHOLD]) == 0:
				continue
			exp_list.append( tuple( [s[NAME_COLUMN], log_values[ control_column - 1 ], log_values[target_column - 1]] ) )
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
def changeRatio( exp_list ):
	ratio_list = []
	for t in exp_list:
		ratio = t[2] - t[1]
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

##
# return columun numbers for the control column and the experimented condition
# target_column: Integer of the analyzed condition
##
def getColumnNumbers( target_column ):
	base_col = -1; obj_cols = None;
	if target_column in EGF_COLS: # if control is EGS
		base_col = 1
		obj_cols = EGF_COLS
	else: # if control is HRG
		base_col = 2
		obj_cols = HRG_COLS
	return base_col, obj_cols

def run( exp_file, output_file, target_column ):
	control_column, obj_cols = getColumnNumbers( target_column ) # get columun numbers for the control column and the experimented condition
	exp_list = readExpFile( exp_file, control_column, target_column, obj_cols )
	ratio_list = changeRatio( exp_list )
	ratio_list = mergeID( ratio_list ) # reduce the duplicate gene names
	output( ratio_list, output_file )
	sys.stdout.write("%s genes are included.\n" % len(ratio_list))

if __name__ == "__main__":
	usage = "usage: %prog expression-file output-file"
	p = OptionParser(usage = usage)
	# Option to set control column
	p.add_option('-t', '--target', dest = "target_column",
				 help = "Analyzed expression column number.")
		
	opts, args = p.parse_args()
	
	# check arguments
	if len(args) < 2:
		sys.stderr.write("Error: input expression-file and output-file\n")
		sys.exit()
		
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
		
	run( args[0], args[1], target_column )
