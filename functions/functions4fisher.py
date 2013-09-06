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
# This source include calculate combination
# and function f which is defined in paper (C(n1, x)/C(n0+n1, x)).
# @author Terada, 26, June, 2011
# @editor Terada, 24, Feb, 2012,
#     If x > n1 in caluclation MASL, then raise error and exit (This case does not treat the case).
# @editor Terada, 16, Apr, 2013,
#     Acceletate of the calculation P-value by storing the calculated P-value.

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
##
class FunctionOfX(fs.FunctionsSuper):
	def __init__(self, transaction_list, row_size):
		fs.FunctionsSuper.__init__(self)
		self.__t_size = len(transaction_list) # all transaction size
		self.__f_size = self.sumValue(transaction_list) # transaction size which have flag = 1 (n1)
		self.__pvalTable = pvalTable.PvalTable( row_size ) # P-value table
		self.__occrTable = pvalTable.PvalTable( row_size ) # occurence table for calculate P-value
		self.calTime = 0 # Total number of calculate P-value
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
		
		# check the support size.
		# If support size larger than half of all data size, raise error.
		# Because this version does not treat x > (n1+n0)/2.
		"""
		if self.__f_size > (self.__t_size/2):
			sys.stderr.write("Reduce the number of rows whose second column is 1 to be less than the half of the total rows.\n")
			sys.stderr.write("This code cannot analyze this case.\n")
			sys.exit()
		"""

			
	def getN1(self):
		return self.__f_size
	
	def funcF(self, x):
		n1 = self.__f_size
		n1_n0 = self.__t_size # the number of all genes
#		print "x: " + str(x)
		# If x > n1, then print error and exit calculation.
		# (MASL may not follow the monotonic decrease.)
		all_x = n1_n0 - n1
#		print "all transaction size: " + str(n1_n0)
#		print "n1+n0-x: " + str(all_x)
#		print "n1: " + str(n1)
		# if x > n1+n0-x, calculate x=n1+n0-x
		
		# x < n1 <= n0
		if x <= n1:
			ans = self.__probability(x, x)
#			print "f(" + str(x) + ") = " + str(ans)
			return ans
		# n1 <= x <= n0
		elif x <= all_x:
			ans = self.__probability(x, n1)
			return ans
		# x > n0:
		else:
#			print x
#			print n1
#			print all_x
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
#		print ovalues
		total_row1 = sum( ovalues[0] )
		p = self.__pvalTable.getValue( total_row1, ovalues[0][0] )
		if p < 0: # calculate P-value and save to the table
			p0 = self.__probability(total_row1, ovalues[0][0])
#			print "p0: " + str(p0)
			p = 0
			a = ovalues[0][0]
			while ( a <= total_row1 ) and ( a <= total_col1 ):
				pa = self.__probability(total_row1, a)
				p = p + pa
				a = a + 1
			self.__pvalTable.putValue( total_row1, ovalues[0][0], p )
			self.calTime = self.calTime + 1
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
			time = 0
			p = 1.0
			while time < a:
#				print "c(n1, a): " + str(n1-time) + "/" + str(a-time)
				p = p*(n1-time)/(a-time) # c(n1, a)
#				print "c(n1+n0, x)" + str(n1_n0-time) + "/" + str(a_b-time)
				p = p*(x-time)/(n-time) # c(n1+n0, x)
				time = time + 1
			time = 0
			while time < b:
#				print "c(n0, b): " + str(n0-time) + "/" + str(b-time)
				p = p*(n0-time)/(b-time) # c(n0, b)
				minus_denominator = a+time
#				print "c(n1+n0, x): " + str(n1_n0-minus_denominator) + "/" + str(a_b-minus_denominator)
				p = p*(x-minus_denominator)/(n-minus_denominator) # c(n1+n0, x)
				time = time + 1
#				print p
			self.__occrTable.putValue( x, a, p )
		return p

def run(xls_file, value_file, itemset_str_lst, delimiter):
	transaction_list, columnid2name, lcm2transaction_id = readFile.readFiles(xls_file, value_file, delimiter)
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


if __name__ == "__main__":
	if (len(sys.argv) < 4):
		sys.stderr.write("Error: functions4fisher.py [item-file] [value-file] [itemset]\n")
		sys.exit()

	xls_file = sys.argv[1]
	value_file = sys.argv[2]
	itemset_str_lst = sys.argv[3].split(',')
	delimiter = ','
	p_value, down_size = run(xls_file, value_file, itemset_str_lst, delimiter)
