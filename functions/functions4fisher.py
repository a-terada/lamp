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
# Calculate the P-value and lower bound of the combination. 
# @author Terada, 26, June, 2011
# @editor Terada, 24, Feb, 2012,
#     If x > n1 in caluclation MASL, then raise error and exit (This case does not treat the case).
# @editor Terada, 16, Apr, 2013,
#     Acceletate of the calculation P-value by storing the calculated P-value.
# @editor Terada, 11, Mar, 2015,
#     Implement computation of the 'less' and 'two-sided' Fisher's exact test. 
from __future__ import division
import sys, os
from . import functionsSuper as fs
from . import pvalTable

pardir = os.path.dirname(os.path.dirname(os.path.abspath( __file__ )))
sys.path.append(pardir)


##
# Define class
# This class calculate function f that define in paper. (C(n1, x)/C(n0+n1, x))
# transaction_list: list of transactions
##
class FunctionOfX(fs.FunctionsSuper):
	def __init__(self, transaction_list, row_size, alternative):
		fs.FunctionsSuper.__init__(self)
		self.__t_size = len(transaction_list) # all transaction size
		self.__f_size = self.sumValue(transaction_list) # transaction size which have flag = 1 (n1)
		self.__pvalTable = pvalTable.PvalTable( row_size ) # P-value table
		self.__occrTable = pvalTable.PvalTable( row_size ) # occurence table for calculate P-value
		self.calTime = 0 # Total number of calculate P-value
		self.alternative = alternative # alternative hypothesis. greater or less -> 1, two.sided -> 0.
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

	
	def funcF(self, x):
		n1 = self.__f_size
		n1_n0 = self.__t_size # the number of all genes
		# If x > n1, then print error and exit calculation.
		# (MASL may not follow the monotonic decrease.)
		all_x = n1_n0 - n1
		
		# x < n1 <= n0
		if x <= n1:
			ans = self.__probability(x, x)
			return ans
		# n1 <= x <= n0
		elif x <= all_x:
			ans = self.__probability(x, n1)
			return ans
		# x > n0:
		else:
			sys.stderr.write("Error: x > n1, n0. This code cannot consider this case.\n")
			sys.exit()
	
			
	##
	# Calculate p-value by using fisher's exact test.
	# transaction_list: List of transactions
	# flag_itemset_id: Transactions which have items
	##
	def calPValue(self, transaction_list, flag_transactions_id):
		ovalues = self.contingencyTable( transaction_list, flag_transactions_id, self.__t_size, self.__f_size )
		total_col1 = self.__f_size
		total_row1 = sum( ovalues[0] )
		p = self.__pvalTable.getValue( total_row1, ovalues[0][0] )
		if p < 0: # calculate P-value and save to the table
			p0 = self.__probability(total_row1, ovalues[0][0])
#			sys.stdout.write("p0: %s\n" % p0)
			p = p0; pos_max = min( total_row1, total_col1 )
			# when the alternative hypothesis is "two.sided",
			# the lower case probability is cumulated. 
			if self.alternative < 1:
				a = 0
				# cumulate the lower case probability.
				while a < ovalues[0][0]:
					pa = self.__probability( total_row1, a )
#					sys.stdout.write( "x: %d, a:%d, pa: %s\n" % (total_row1, a, pa) )
					if (pa - p0 > 1.E-16): # pa > p0
						break
					p = p + pa
					a = a + 1
				# cumulate the upper case probability.
				a = pos_max
				while ( a > ovalues[0][0] ):
					pa = self.__probability( total_row1, a )
#					sys.stdout.write( "x: %d, a:%d, pa: %s\n" % (total_row1, a, pa) )
					if (pa - p0 > 1.E-16): # pa > p0
						break
					p = p + pa
					a = a - 1
			# when the alternative hypothesis is "greater" or "less",
			# the higher/less case probability is cumulated.  
			else:
				a = ovalues[0][0] + 1
				while ( a <= pos_max ):
					pa = self.__probability( total_row1, a )
					p = p + pa
					a = a + 1
			self.__pvalTable.putValue( total_row1, ovalues[0][0], p )
			self.calTime = self.calTime + 1
#		sys.stdout.write( "x: %d, a:%d, p: %s\n" % (total_row1, ovalues[0][0], p) )
		return p, ovalues[0][0]
	
	##
	# Calculate probability of occurrence probability about table.
	# x: total of the first row (the number of targetting gene)
	# a: the top-left of the table
	# Return C(n1, a)*C(n0, b)/C(n1+n0, a+b)
	##
	def __probability(self, x, a):
		p = self.__occrTable.getValue( x, a )
		if p < 0:
			n = self.__t_size
			n1 = self.__f_size
			n0 = n - n1
			b = x - a
			i = 0
			p = 1.0
			while i < a:
				p = p*(n1 - i)/(a - i) # c(n1, a)
				p = p*(x - i)/(n - i) # c(n1+n0, x)
				i = i + 1
			i = 0
			while i < b:
				p = p*(n0 - i)/(b - i) # c(n0, b)
				minus_denominator = a + i
				p = p*(x - minus_denominator)/(n - minus_denominator) # c(n1+n0, x)
				i = i + 1
			self.__occrTable.putValue( x, a, p )
		return p

def run(xls_file, value_file, itemset_str_lst, delimiter, alternative):
	global readFile
	import readFile
	transaction_list, columnid2name = readFile.readFiles(xls_file, value_file, delimiter)
	max_lambda = maxLambda(transaction_list)
	if alternative < 0:
		global lamp
		from lamp import reverseValue
		transaction_list = reverseValue( transaction_list, "fisher" )
	func = FunctionOfX(transaction_list, max_lambda, abs( alternative ) )
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
	sys.stdout.write("p-value: %s (N: %s, n1: %s, x: %s, a: %s)\n"
					 % (p_value, n, n1, len(flag_transactions_id), stat_value))
	return (p_value, len(flag_transactions_id))

##
# Return max lambda. That is, max size itemset.
##
def maxLambda(transaction_list):
	# Count each item size
	item_sizes = {}
	for t in transaction_list:
		for item in t.itemset:
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


if __name__ == "__main__":
	if (len(sys.argv) < 4):
		sys.stderr.write("Error: functions4fisher.py [item-file] [value-file] [itemset]\n")
		sys.exit()

	xls_file = sys.argv[1]
	value_file = sys.argv[2]
	itemset_str_lst = sys.argv[3].split(',')
	delimiter = ','
	p_value, down_size = run(xls_file, value_file, itemset_str_lst, delimiter)
