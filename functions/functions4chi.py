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

# Define fuctions that is used in multiple_test.py
# This source includes calculate P-value and MASL of the chi-square test.
# @author Terada, 16, Apr, 2013
# @editor Terada, 11, Mar, 2015,
#     Implement computation of the 'less' and 'two-sided' Chi-square test. 

from __future__ import division
import sys, os
import functionsSuper as fs
import pvalTable

pardir = os.path.dirname(os.path.dirname(os.path.abspath( __file__ )))
sys.path.append(pardir)
import readFile

##
# Define class
# This class calculate function f that define in paper. (C(n1, x)/C(n0+n1, x))
# transaction_list: list of transactions
# alternative: alternative hypothesis, 1 -> "greater" or "less", 0 -> two.sided
##
class FunctionOfX(fs.FunctionsSuper):
	def __init__(self, transaction_list, row_size, alternative):
		fs.FunctionsSuper.__init__(self)
		self.__t_size = len(transaction_list) # all transaction size
		self.__f_size = self.sumValue(transaction_list) # transaction size which have flag = 1 (n1)
		self.alternative = alternative # alternative hypothesis. greater or less -> 1, two.sided -> 0.
		self.__pvalTable = pvalTable.PvalTable( row_size ) # P-value table
		self.__chiTable = pvalTable.PvalTable( row_size ) # P-value table
		if self.__f_size == 0:
			sys.stdout.write("Error: There is no up-regulate gene.\n")
			sys.exit()
		# Check the transaction value.
		# If the value is not 1 or 0, raise error.
		# Because fisher's exact test does not handle numerical value.
		for t in transaction_list:
			if not (t.value == 1.0 or t.value == 0.0):
				sys.stderr.write("Error: \"" + t.name + "\" value is " + str(t.value)+".\n")
				sys.stderr.write("       But value is 1 or 0 if you test by fisher's exact test.\n")
				sys.exit()
				
	
	def getN1(self):
		return self.__f_size

	def getAllSize(self):
		return self.__t_size

	##
	# calclate MASL
	##
	def funcF(self, x):
		p1 = p2 = 1.0
		chi1 = chi2 = 0.0
		total_row1 = self.__f_size
		total = self.__t_size
		# when x < n_u
		if x < total_row1:
			ovalues = [[x, 0], [total_row1 - x, total - total_row1]]
			p1, chi1 = self.__probabilityTable( ovalues )
			ovalues = [[0, x], [total_row1, total - total_row1 - x]]
			p2, chi2 = self.__probabilityTable( ovalues )
		# when x >= n_u
		else:
			ovalues = [[total_row1, x-total_row1], [0, total - x]]
			p1, chi1 = self.__probabilityTable( ovalues )
			ovalues = [[0, x], [total_row1, total - total_row1 - x]]
			p2, chi2 = self.__probabilityTable( ovalues )
		if self.alternative == 0:
			p1 = min( p1 * 2., 1.0 )
			p2 = min( p2 * 2., 1.0 )
		if p1 < p2:
			return p1
		else:
			return p2
					
	##
	# Calculate p-value by using chi-square test.
	# transaction_list: List of transactions
	# flag_itemset_id: Transactions which have items
	# alternative: alternative hypothesis, 1 -> greater, 0 -> two-sided, -1 -> less.
	##
	def calPValue(self, transaction_list, flag_transactions_id):
		ovalues = self.contingencyTable( transaction_list, flag_transactions_id, self.__t_size, self.__f_size )
		total_col1 = self.__f_size
		total_row1 = sum( ovalues[0] )
		p = self.__pvalTable.getValue( total_row1, ovalues[0][0] )
		chi = self.__chiTable.getValue( total_row1, ovalues[0][0] )
		if p < 0: # calculate P-value and save to the table
			p, chi = self.__probabilityTable(ovalues)
			if (self.alternative > 0): 
				if (ovalues[0][0] < min(self.__f_size, total_row1)/2):
					p = 1. - p
			# when the alternative hypothesis is "two.sided", 
			# the P-value is doubled. 
			else:
				p = min( p * 2., 1.0 )
			self.__pvalTable.putValue( total_row1, ovalues[0][0], p )
			self.__chiTable.putValue( total_row1, ovalues[0][0], chi )
#		sys.stdout.write( "x: %d, a:%d, chi: %s, p: %s\n" % (total_row1, ovalues[0][0], chi, p) )
		return p, ovalues[0][0]
	
	def __calMeans(self, ovalues):
		total = self.__t_size
		total_col1 = self.__f_size # the number of all flag 1 transaction (n1)
		total_col2 = total - total_col1 # the number of all flag 0 transactio (n0)
		total_row1 = ovalues[0][0] + ovalues[0][1]
		total_row2 = ovalues[1][0] + ovalues[1][1]
		means = []; means.append([0]*2); means.append([0]*2)
		means[0][0] = float(total_row1 * total_col1) / total
		means[0][1] = float(total_row1 * total_col2) / total
		means[1][0] = float(total_row2 * total_col1) / total
		means[1][1] = float(total_row2 * total_col2) / total
		return means
	
	##
	# Calculate probability of occurrence probability about table.
	# a: Top left of table
	# b: Top right of table
	# n1: Sum of top and bottom left (a + c)
	# n0: Sum of top and bottom right (b + d)
	##
	def __probabilityTable(self, ovalues):
		means = self.__calMeans(ovalues) # calculate the exception value
		
		# Yate continuity correction
		yate_corr = 0
		for i in means:
			for j in i:
				if j < 5:
					yate_corr = 0.5
					break
		
		chi = 0
		for i in xrange(0, len(ovalues)):
			row = ovalues[i]
			for j in xrange(0, len(row)):
				chi = chi + (abs(row[j] - means[i][j]) - yate_corr)**2/means[i][j]
		return self.__chi2pval( chi ), chi

	def __chi2pval(self, chi):
		if (chi == 0.0):
			return 1.0
		else: # dimension = 1
			return (self.stdNorDistribution(chi**0.5))

"""
def maxLambda(transaction_list):
	# Count each item size
	item_sizes = {}
	for t in transaction_list:
#		print t.itemset
		for item in t.itemset:
#			print item
			# If item does not exist in item_size, then make mapping to 0
			if not item_sizes.has_key(item):
				item_sizes[item] = 0
			item_sizes[item] = item_sizes[item] + 1
	# Get max value in item_sizes
	max_value = 0
	for i in item_sizes.itervalues():
		if i > max_value:
			max_value = i
	return max_value


def run(xls_file, value_file, itemset_str_lst, delimiter):
	transaction_list, columnid2name = readFile.readFiles(xls_file, value_file, delimiter)
	max_lambda = maxLambda(transaction_list)
	func = FunctionOfX(transaction_list, max_lambda)
	colname2id_dict = readFile.colname2id(columnid2name)

	itemset = set()
	for i in itemset_str_lst:
		item_id = colname2id_dict[i]
		itemset.add(item_id + 1)
		
	flag_transactions_id = []
	for i in xrange( len(transaction_list) ):
		t = transaction_list[i]
		if len( itemset & t.itemset ) == len(itemset):
			flag_transactions_id.append( i )
	p_value, stat_value = func.calPValue(transaction_list, flag_transactions_id)
	n = len(transaction_list)
	n1 = func.getN1()
	sys.stdout.write("p-value: %s (N: %s, n1: %s, x: %s, chi: %s)\n"
					 % (p_value, n, n1, len(flag_transactions_id), stat_value))
	return p_value, len(flag_transactions_id)

if __name__ == "__main__":
	if (len(sys.argv) < 4):
		sys.stderr.write("Error: functions4chi.py [item-file] [value_file] [itemset]\n")
		sys.exit()
	
	xls_file = sys.argv[1]
	value_file = sys.argv[2]
	itemset_str_lst = sys.argv[3].split(',')
	delimiter = ','
	p_value, down_size = run(xls_file, value_file, itemset_str_lst, delimiter)
"""
