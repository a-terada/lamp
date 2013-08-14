#!/usr/bin/env python

# Check gene-IDs between target-file expression-file.

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
