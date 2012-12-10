#!/usr/bin/env python

# Convert ProbeID to OtherID using Affymetrix annotation file.
# This code convert to EntrezGene (column number is set by TO_COLUMN_ID).
# If a probe maps several IDs, this code selects the first ID in the annotation file.
# @author A. Terada

__author__ = "Aika TERADA"

import sys, os
from optparse import OptionParser

# Definitions for map-file
ORIG_COLUMN_ID = 0 # Column index to converted ID
TO_COLUMN_ID = 18 # Column index to new ID
MAP_SEPARATOR = '","' # Separator for the line
MAP_SEPARATOR_COLUMN = "///" # Separator for a column.
# Definitions for converted-file
CONV_ORIG_COLUMN_ID = 0 # Prob ID column number
CONV_SEPARATOR = "\t" # Fils separator

def readMapFile(filename):
	old2new = {} # the mapping of prob ID for EntrezGene ID
	try:
		f = open(filename, 'r')
		line = ""
		for line in f:
			line = line[:-1]
			# If line is the comment line, skip
			if (line.startswith('#')):
				continue
			else:
				s = line.split(MAP_SEPARATOR)
				old_id = s[ORIG_COLUMN_ID].strip('"') # String of an old ID
				new_ids = map(str.strip, s[TO_COLUMN_ID].strip('"').split('///')) # List of new IDs
				if not (new_ids[0] == "---"):
					old2new[old_id] = new_ids[0]
		f.close()
		return old2new
	except IOError, e:
		sys.stderr.write("Error in read %s\n" % filename)
		sys.exit()

def convertID( old2new, converted_file, output_file ):
	try:
		fr = open( converted_file, 'r' )
		fw = open( output_file, 'w' )
		line = ""
		total_size = 0; converted_size = 0;
		for line in fr:
			if not ( line.startswith("!") ):
				total_size += 1
				line = line[:-1]
				s = line.split(CONV_SEPARATOR)
				old_id = s[CONV_ORIG_COLUMN_ID].strip('"')
				if old_id in old2new:
					converted_size += 1
					new_id = old2new[old_id]
					fw.write("%s" % new_id)
					for i in s[1:]:
						fw.write("%s%s" % (CONV_SEPARATOR, i))
					fw.write("\n")
		fr.close()
		fw.close()
		return converted_size, total_size
	except IOError, e:
		sys.stderr.write("Error in convert %s\n" % filename)
		sys.exit()

def run( map_file, converted_file, output_file ):
	sys.stderr.write("read files ...\n")
	old2new = readMapFile(map_file)
	sys.stderr.write("convert ID...\n")
	converted_size, total_size = convertID( old2new, converted_file, output_file )
	sys.stderr.write("%s/%s IDs were converted\n" % (converted_size, total_size))

if __name__ == "__main__":
	usage = "usage: %prog map-file converted-file output-file"
	p = OptionParser(usage = usage)
	
	if (len(sys.argv) < 4):
		sys.stderr.write("Error: input map-file, converted-file and ontput-file.\n")
		sys.exit()
	
	opts, args = p.parse_args()
	
	# check the file exist.
	if not os.path.isfile(args[0]):
		sys.stderr.write("IOError: No such file: \'" + args[0] + "\'\n")
		sys.exit()
	if not os.path.isfile(args[1]):
		sys.stderr.write("IOError: No such file: \'" + args[1] + "\'\n")
		sys.exit()
	
	run(args[0], args[1], args[2])
