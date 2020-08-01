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

# Define methods to read transaction and flag files
# @author Terada 26, June, 2011
# @editor aika 1, Nov. 2011
#    Developing the readValueFile for handling the value.
# @editor aika 10, Nov. 2011
#    Change readTransactionFile. If file has space the front and the end, remove them.
# @editor aika, 11, Mar. 2014
#    Change readFiles for keeping transaction ID.

import sys, transaction, csv

##
# Read transaction and flag file and return transaction matrix.
# transaction_file: 
##
def readFiles(transaction_file, value_file, delm):
	transaction_list, gene2id, columnid2name = readTransactionFile(transaction_file, delm)
	transaction_list = readValueFile(value_file, transaction_list, gene2id, delm)
	transaction_list.sort() # sort transaction_list according to transaction_value
	# Generate IDs to transactions
	for i in xrange( 0, len(transaction_list) ):
		t = transaction_list[i]
		t.id = i
	# check transaction names two
	checkTransName(transaction_list, transaction_file) # check transaction names two
	return transaction_list, columnid2name

##
# Read transaction file and return list of transactions.
# item ID is integer value and begin from 0.
##
def readTransactionFile(transaction_file, delm):
	transaction_list = []
	gene2id = {} # dictionary that gene name -> transaction ID
	columnid2name = [] # list about mapping column id to column name
	gene_set = set([])
	line_num = 0; col_size = -1
	try:
		f = open( transaction_file, 'rU' )
		for row_list in csv.reader( f, delimiter = delm ):
			line_num = line_num + 1
			# If line is header line, read column name
			if line_num == 1:
				col_size = len( row_list )
				for i in range(1, col_size):
					colname = row_list[i]
					columnid2name.append(colname)
				continue
			
			t_name = row_list[0]
			if t_name in gene_set:
				sys.stderr.write("Error: %s is contained two or more times in %s.\n" \
								 % (t_name, transaction_file))
				sys.exit()
			# check the number of columns
			if len( row_list ) != col_size:
				sys.stderr.write("Error in %s\n" % transaction_file)
				sys.stderr.write("    The header line contains %s columns, while line %s contains %s columns.\n" \
								 % (col_size, line_num, len( row_list )))
				sys.exit()
				
			gene_set.add(t_name)
			t = transaction.Transaction(t_name)
			gene2id[t_name] = len(transaction_list)
			for i in range(1, len(row_list)):
				flag = int(row_list[i])
				if flag == 1:
					t.addItem(i)
				elif flag == 0:
					continue
				else:
					sys.stderr.write("line %s in \'%s\' contains the value neither 0 or 1.\n" \
									 % (line_num, transaction_file) )
					sys.exit()
			transaction_list.append(t)
	except IOError as e:
		sys.stderr.write("Error: %s\n" % e)
		sys.exit()
	return transaction_list, gene2id, columnid2name


##
# Read flag file and add information about flags to transaction list.
# value_file: Read flag file.
#     The column1 is gene name and the column2 is value (for example, gene expression).
# transaction_list: List of transactions. This is made of readTransactionFile method.
# gene2id: Dictionary that key indicates gene name and value is transaction ID(location of list)
##
def readValueFile(value_file, transaction_list, gene2id, delm):
	line_num = 0; gene_set = set([])
	try:
		f = open( value_file, 'rU' )
		for row_list in csv.reader( f, delimiter = delm ):
			line_num = line_num + 1
			if (row_list[0].startswith("#")):
				continue

			# This error raises if value file contains more than two columns.
			if not len(row_list) == 2:
				sys.stderr.write("Error: line %s in %s.\n" % (line_num, value_file) )
				sys.stderr.write("       value-file should contain two columns.\n")
				sys.exit()

			genename = row_list[0].strip()
			exp_value = row_list[1].strip()

			# This error raises if value cannot be converted to float.
			if not isFloat( exp_value ):
				sys.stderr.write("Error: line %s in %s.\n" % (line_num, value_file))
				sys.stderr.write("       \'%s\' could not be converted string to float.\n" % exp_value)
				sys.exit()
			# This error raises if the identical keys are contained more than two times.
			if genename in gene_set:
				sys.stderr.write("Error: %s is contained two or more times in %s.\n" \
								 % (genename, value_file))
				sys.exit()
			# This error raise if gene does not include in itemset file.
			if not genename in gene2id:
				sys.stderr.write("Error:line %s in %s.\n" % (line_num, value_file))
				sys.stderr.write("      \'%s\' is not contained in itemset file.\n" % genename)
				sys.exit()
				
			gene_set.add(genename)
			exp_value = float( exp_value )
			geneid = gene2id[genename]
			t = transaction_list[geneid]
			t.value = exp_value
	except IOError as e:
		sys.stderr.write("Error: %s cannot be found.\n" % value_file)
		sys.exit()
	return transaction_list

##
# Check two files transaction name are same
# If transaction contains only one file, flag is -1.
##
def checkTransName(transaction_list, transaction_file):
	for t in transaction_list:
		if t.value == None:
			sys.stderr.write("\"%s\" only appears in %s\n" % (t.name, transaction_file))
			sys.exit()


##
#
##
def colname2id(columnid2name):
	colname2id_dict = {}
	index = 0
	for i in columnid2name:
		colname2id_dict[i] = index
		index = index + 1
	return colname2id_dict

def isFloat( value_str ):
	try:
		float( value_str )
		return True
	except ValueError:
		return False
