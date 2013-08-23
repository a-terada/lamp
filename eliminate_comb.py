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

# Remove the redundant combinations in the result of LAMP.
# @author Terada 3, Aug., 2013


import sys, os
from optparse import OptionParser
from operator import itemgetter

# read the result of LAMP.
# filename: the result filename of LAMP
def readResult( filename ):
	try:
		detections_list = [] # Store the combinations 
		meta_line_list = [] # Store the meta lines
		time_line = "" # Store the line containing running time
		f = open(filename, 'r')
		flag_broken = True # a flag to decide whether the file is broken or not.
		flag_comb = False
		detection_size = 0 # number of detection sets.
		for line in f:
			line = line[:-1]
			if ( line.startswith("Time (sec.): ") ):
				time_line = line
				flag_broken = False
				continue
			if not flag_comb:
				meta_line_list.append( line )
				if ( line.startswith("Rank") ):
					flag_comb = True
				continue
			# [0]: Rank, [1]: Raw p-value, [2]: Adjusted P-value,
			# [3]: detections, [4]: arity, [5]: support, [6]: statistic_score.
			detections = line.split('\t')
			detections[1] = float(detections[1])
			detections[2] = float(detections[2])
			detections[4] = int(detections[4])
			detections_list.append( tuple( detections ) )
		# if the file does not have results, output error.
		if (flag_broken):
			sys.stderr.write("%s is broken.\n" % filename)
		f.close()
		return detections_list, meta_line_list, time_line
	except IOError, e:
		sys.stderr.write("'%s' cannot be opened.\n" % filename)
		sys.exit()

# Return whether set_j is subset of set_i.
def isSubset(set_i, set_j):
#	print set_i
#	print set_j
	diff_set = set_i - set_j
	if (len(diff_set) == (len(set_i) - len(set_j))):
#		print "subset"
		return True
	else:
		return False

# Eliminate the redundant combinations.
# detections_list: a list of combinations made by readResult function.
def mergeResult( detections_list ):
	upper_list = []; merged_list = []
	for detection in detections_list:
		detect_set = set(detection[3].split(','))
#		print detect_set
		flag = True
		for i in upper_list:
			if ( isSubset(detect_set, i) ) or ( isSubset(i, detect_set) ):
				flag = False
				break
		if flag:
			merged_list.append( detection )
		upper_list.append( detect_set )
	return merged_list

# Sort detections_list by P-value.
# If comb1 and comb2 have equal P-values, the larger combination gives high rank.
def sortComb( detections_list ):
	detections_list = sorted( detections_list, key = itemgetter( 1, 4 ) )
	return detections_list

# output result
def output( out_filename, detections_list, meta_line_list, time_line ):
	if len( out_filename ) > 0:
		try:
			sys.stdout = open(out_filename, 'w')
		except IOError, e:
			sys.stderr.write("Error: Cannot output to %s'\n" % out_filename)
			sys.exit()
	sys.stdout.write("# Non-redundant combinations\n")
	for i in meta_line_list:
		sys.stdout.write("%s" % i)
		if i.startswith("# # of significant combinations: "):
			sys.stdout.write(" -> %d" % len(detections_list))
		sys.stdout.write("\n")
	rank = 0
	for detect in detections_list:
		rank = rank + 1
		sys.stdout.write("%d" % rank)
		for i in detect[1:]:
			sys.stdout.write("\t%s" % i)
		sys.stdout.write("\n")
	sys.stdout.write("%s\n" % time_line)
	sys.stdout.close()
	

# in_filename: the filename of LAMO
# out_filename: the filename to output the eliminated result.
#               When out_filename is "", the result is printed as standard output.
def run( in_filename, out_filename ):
	detections_list, meta_line_list, time_line = readResult( in_filename )
	detections_list = sortComb( detections_list )
	merged_list = mergeResult( detections_list )
	output( out_filename, merged_list, meta_line_list, time_line )
	

if __name__ == "__main__":
	usage = "usage: %prog filename"
	p = OptionParser(usage = usage)

	p.add_option('-o', dest = "output_filename", default = "",\
				 help = "The file name to output the eliminated result.\n")

	opts, args = p.parse_args()
	
	# check the arguments
	if not (len(args) == 1):
		sys.stderr.write("Error: %s filename\n" % __file__)
		sys.exit()
	
	run( args[0], opts.output_filename )
