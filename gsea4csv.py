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
# The dictionary map gene name and interacted TFs list.
##
def readGmtFile( gmt_file ):
	gene_dict = {} # key: gene, value: TFs list
	all_motif_lst = [] # TFs list which appear in gme file
	for line in open(gmt_file, 'r'):
		line = line[:-1]
		s = line.split('\t')
		tf_set = s[0] # TFs set include represent motif
		motif = ""
		# If the TF does not known
		if (tf_set.find("$") < 0):
			motif = tf_set
		else:
			s_tf_set = tf_set.split('$')
			motif = s_tf_set[1].split("_")[0]
			
		if not motif in all_motif_lst:
			all_motif_lst.append(motif)
		
		# mapping gene and the tf motif
		interacted_genes = s[2:] # interacted genes
		for gene in interacted_genes:
			try:
				motif_set = gene_dict[gene]
				if not motif in motif_set:
					motif_set.append(motif)
					gene_dict[gene] = motif_set
			except KeyError, e:
				gene_dict[gene] = [motif]
	return all_motif_lst, gene_dict
	
##
# output CSV format.
# all_tf_list: List of all TFs included GMT file.
# gene_dict: [key] gene, [value] regulate TF list.
# output_file: The filename to print CSV format.
##
def outCSVFormat( all_tf_lst, gene_dict, output_file ):
	try:
		f = open( output_file, 'w' )
		# output header
		f.write( "#EntrezGene" )
		for tf in all_tf_lst:
			f.write( ",%s" % tf )
		f.write("\n")
		
		for gene in gene_dict:
			f.write( "%s" % gene )
			regulated_tfs = gene_dict[ gene ]
			for tf in all_tf_lst:
				if tf in regulated_tfs:
					f.write( ",1" )
				else:
					f.write( ",0" )
			f.write("\n")
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
