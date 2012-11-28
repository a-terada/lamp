#!/usr/local/bin/python
# -*- coding: utf-8 -*-

# Define methods to read transaction and flag files
# @author aika 26, June, 2011
# @editor aika 1, Nov. 2011
#   developing the readValueFile for handling the value.
# @editor aika 10, Nov. 2011
#   change readTransactionFile. If file has space the front and the rear, remove them.

import sys, transaction

class ReadFileError(Exception):
	def __init__(self, e):
		print "ReadFileError: " + e
		
##
# read transaction and flag file
# and return transaction matrix.
# transaction_file: 
##
def readFiles(transaction_file, value_file):
	transaction_list, gene2id, columnid2name, lcm2transaction_id = readTransactionFile(transaction_file)
	transaction_list = readValueFile(value_file, transaction_list, gene2id)
	# check transaction names two
	checkTransName(transaction_list, transaction_file)
	return transaction_list, columnid2name, lcm2transaction_id

##
# read transaction file and return list of transactions.
# item ID is integer value and begin from 0.
##
def readTransactionFile(transaction_file):
	transaction_list = []
	gene2id = {} # dictionary that gene name -> transaction ID
	columnid2name = [] # list about mapping column id to column name
	lcm2transaction_id = [] # line number in LCM arguments file.
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
#				t.addItem(i-1)
		if ( len(t.itemset) > 0 ):
			lcm2transaction_id.append(len(transaction_list))
		transaction_list.append(t)
	return transaction_list, gene2id, columnid2name, lcm2transaction_id


##
# read flag file and add information about flags to transaction list.
# value_file: read flag file. The column1 is gene name and the column2 is vlue (for example, gene expression).
# transaction_list: list of transactions. This is made of readTransactionFile method.
# gene2id: dictionary that key indicates gene name and value is transaction ID(location of list)
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
			
		# this error raise if gene does not include in itemset file.
		except KeyError, e:
			e_out =  "line " + str(line_num) + ", \"" + line[:-1] + "\"" + "\n       " + str(e) + " does not exist in item set file."
			raise ReadFileError, e_out
		except ValueError, e:
			e_out = "line " + str(line_num) + ", \"" + line[:-1] + "\"" + "\n       " + str(e) + " does not numerical value."
			raise ReadFileError, e_out
	return transaction_list

##
# check two files transaction name are same
# If transaction contains only one file, flag is -1
##
def checkTransName(transaction_list, transaction_file):
	for t in transaction_list:
		if t.value == None:
			e_out = "\"" + t.name + "\" only appears in " + transaction_file
			raise ReadFileError, e_out

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


# readFiles(sys.argv[1], sys.argv[2])
