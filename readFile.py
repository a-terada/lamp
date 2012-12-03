#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Define methods to read transaction and flag files
# @author Terada 26, June, 2011
# @editor aika 1, Nov. 2011
#    Developing the readValueFile for handling the value.
# @editor aika 10, Nov. 2011
#    Change readTransactionFile. If file has space the front and the end, remove them.

import sys, transaction

class ReadFileError(Exception):
	def __init__(self, e):
		print "ReadFileError: " + e
		
##
# Read transaction and flag file and return transaction matrix.
# transaction_file: 
##
def readFiles(transaction_file, value_file):
	transaction_list, gene2id, columnid2name = readTransactionFile(transaction_file)
	transaction_list = readValueFile(value_file, transaction_list, gene2id)
	transaction_list.sort() # sort transaction_list according to transaction_value
	# check transaction names two
	checkTransName(transaction_list, transaction_file) # check transaction names two
	lcm2transaction_id = makeLCM2TransactionList(transaction_list) # make transaction id list to convert LCM result to transaction_list
	return transaction_list, columnid2name, lcm2transaction_id

##
# Read transaction file and return list of transactions.
# item ID is integer value and begin from 0.
##
def readTransactionFile(transaction_file):
	transaction_list = []
	gene2id = {} # dictionary that gene name -> transaction ID
	columnid2name = [] # list about mapping column id to column name
	for line in open(transaction_file, 'r'):
		# If line is header line, read column name
		if line.startswith("#"):
			s = line[:-1].split(',')
			for i in range(1, len(s)):
				colname = s[i]				
				colname = colname.strip() # If column name starts or finishes spaces, remove them.
				columnid2name.append(colname)
			continue
		s = line[:-1].split(',')
		t_name = s[0].strip()
		t = transaction.Transaction(t_name)
		gene2id[t_name] = len(transaction_list)
		for i in range(1, len(s)):
			flag = s[i]
			flag = flag.strip() # If flag includes spaces, remove them.
			if flag == "1":
				t.addItem(i)
		transaction_list.append(t)
	return transaction_list, gene2id, columnid2name


##
# Read flag file and add information about flags to transaction list.
# value_file: Read flag file.
#     The column1 is gene name and the column2 is vlue (for example, gene expression).
# transaction_list: List of transactions. This is made of readTransactionFile method.
# gene2id: Dictionary that key indicates gene name and value is transaction ID(location of list)
##
def readValueFile(value_file, transaction_list, gene2id):
	line_num = 0
	for line in open(value_file, 'r'):
		line_num = line_num + 1
		if (line.startswith("#")):
			continue
		s = line[:-1].split('\t')
		try:
			genename = s[0].strip()
			value = float(s[1].strip())
			geneid = gene2id[genename]
			t = transaction_list[geneid]
			t.value = value
#			print "---"
#			print t.name
#			print t.itemset
#			print t.value
			
		# This error raise if gene does not include in itemset file.
		except KeyError, e:
			e_out =  "line " + str(line_num) + ", \"" + line[:-1] + "\"" + "\n       " + str(e) + " does not exist in item set file."
			raise ReadFileError, e_out
		except ValueError, e:
			e_out = "line " + str(line_num) + ", \"" + line[:-1] + "\"" + "\n       " + str(e) + " does not numerical value."
			raise ReadFileError, e_out
	return transaction_list

##
# Check two files transaction name are same
# If transaction contains only one file, flag is -1.
##
def checkTransName(transaction_list, transaction_file):
	for t in transaction_list:
		if t.value == None:
			e_out = "\"" + t.name + "\" only appears in " + transaction_file
			raise ReadFileError, e_out


##
# Make list to convert LCM result to transaction_list index
##
def makeLCM2TransactionList(transaction_list):
	lcm2transaction_id = []
	for i in range(0, len(transaction_list)):
		t = transaction_list[i]
		if ( len(t.itemset) > 0 ):
			lcm2transaction_id.append(i)
	return lcm2transaction_id

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
