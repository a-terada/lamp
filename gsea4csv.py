#!/usr/bin/env python

# Convert the C3 of GSEA to csv format.
# The made file represents the association between gene and TF.
# "," is deliminated in the output file.
# The first column is a gene ID of Entrez gene.
# From the second columns,  the value indicates whether the gene is regulated the TF.
# If the gene is regulated the TF, then the value is 1, otherwise 0.
# @author aika, 16, Nov., 2011

__author__ = "Aika TERADA"

import sys
from optparse import OptionParser

##
# read GMT file and return gene dictionary.
# The dictionary map gene name and interacted miRNAs list.
##
def readGmtFile( gmt_file ):
	gene_dict = {} # key: gene, value: miRNA list
	all_tf_list = []
	try:
		f = open( gmt_file, 'r' ); line = ""
		for line in f:
			s = line[:-1].split('\t')
			tf = s[0] # TF name
			tf = tf.replace( ',', ' ' ) # Replace comma in TF name to space
			if tf in all_tf_list:
				sys.stderr.write("Error: %s is duplicate\n" % tf)
				sys.exit()
			# mapping gene and the TF
			all_tf_list.append( tf )
			regulated_genes = s[2:] # regulated genes
			for gene in regulated_genes:
				if gene in gene_dict:
					tf_set = gene_dict[ gene ]
					tf_set.append( tf )
				else:
					gene_dict[ gene ] = [ tf ]
		f.close()
		return all_tf_list, gene_dict
	except IOError, e:
		sys.stderr.write("Error in read gmt-file\n")
		sys.exit()

##
# output CSV format.
# all_tf_list: List of all TFs included GMT file.
# gene_dict: [key] gene, [value] regulate TF list.
# output_file: The filename to print CSV format.
##
def outCSVFormat( all_tf_list, gene_dict, output_file ):
	try:
		f = open( output_file, 'w' )
		# output header
		f.write( "#EntrezGene" )
		for tf in all_tf_list:
			f.write( ",%s" % tf )
		f.write( "\n" )
		
		for gene in gene_dict:
			f.write( "%s" % gene )
			regulated_tfs = gene_dict[ gene ]
			for tf in all_tf_list:
				if tf in regulated_tfs:
					f.write( ",1" )
				else:
					f.write( ",0" )
			f.write( "\n" )
	except IOError, e:
		sys.stderr.write("Error in output CSV file.\n")
		sys.exit()

	
def run( gmt_file, output_file ):
	sys.stderr.write( "Read GMT file...\n" )
	all_tf_list, gene_dict = readGmtFile( gmt_file )
	sys.stderr.write( "Make GSV file...\n" )
	outCSVFormat(all_tf_list, gene_dict, output_file)
	sys.stdout.write( "# of TFs: %s, # of genes: %s\n" % ( len(all_tf_list), len(gene_dict)) )

if __name__ == "__main__":
	usage = "usage: %prog gmt_file"
	p = OptionParser(usage = usage)

	opts, args = p.parse_args()
	
	# check arguments
	if (len(args) < 2):
		print "Error: input gmt-file output_file"
		sys.exit()

	run( args[0], args[1] )
