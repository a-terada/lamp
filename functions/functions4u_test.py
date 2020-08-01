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
# This source include method p-value by Mann-Whitne's U-test.
# Function f which represents the minimum p-value. (MASL)
# @author Terada, 8, Nov, 2011
# @editor Terada, 28, Jan, 2011
#         Add main method. When user run this code alone, calculate p-value by U-test.
# @editor Terada, 10, Feb, 2012
#         For speeding up, the binary search list (t_group_y in uValue) cut off the smaller part
#         that find the less value in previous fase.

from __future__ import division
import sys, math, os
from . import functionsSuper as fs

pardir = os.path.dirname(os.path.dirname(os.path.abspath( __file__ )))
sys.path.append(pardir)

##
# Define class
# This class calculate function f that means minimum p-value (MASL).
# transaction_list: list of transactions
##
class FunctionOfX(fs.FunctionsSuper):
	def __init__(self, transaction_list, alternative):
		fs.FunctionsSuper.__init__(self)
		self.__t_size = len(transaction_list) # all transaction size
		self.__transaction_list = transaction_list[:]
		self.alternative = alternative # alternative hypothesis. greater -> 1, less -> -1, two.sided -> 0.
		self.calTime = 0 # Total number of calculate P-value
#		self.__range20_1 = range(1, 21)
#		self.__range20_1.reverse() # the integer list from 20 to 1. this is used by standard normal probability

	##
	# Search the group which t.value less than threshold.
	# threshold: search value
	# tgroup: list to search threshold
	# left_index: start index of tgroup to search
	# right_index: end index of tgroup to search
	##
	def __binarySearch(self, threshold, tgroup, left_index, right_index):
		# print tgroup and threshold
#		print "--- binary search ---"
#		print "threshold: " + str(threshold)
#		for t in tgroup:
#			print t.value,
#		print ""
#
#		print "left:" + str(left_index)
#		print "right:" + str(right_index)

		if (left_index >= len(tgroup)):
			return len(tgroup), len(tgroup)
		
		# compare threshold to min and max value
		if (tgroup[left_index].value > threshold):
			return left_index, left_index
		if (tgroup[right_index].value < threshold):
			return right_index, right_index
		
		# serach the index which all larger indexes are more than threshold.
		mid_index = -1
		while (left_index <= right_index):
			mid_index = (left_index + right_index) // 2
			mid_transaction = tgroup[mid_index]
			
#			print "left: " + str(left_index) + ", right: " + str(right_index)
#			print "mid: " + str(mid_transaction.value)

			# When the mid.value = threshold, finish the search.
			if (mid_transaction.value == threshold):
				break
			# When the check value less than threshod, go to serach right.
			elif (mid_transaction.value < threshold):
				left_index = mid_index + 1
			# When the check value >= threshold, go to search left.
			else:
				right_index = mid_index - 1

#		print "left: " + str(left_index) + ", right: " + str(right_index)
#		print "mid: " + str(mid_index) + " " + str(mid_transaction.value)
		# search the same range of the threshold
		mid_transaction = tgroup[mid_index]
		if (mid_transaction.value == threshold):
			min_index = mid_index
			max_index = mid_index
			min_transaction = tgroup[min_index]
			max_transaction = tgroup[max_index]
			while (min_transaction.value >= threshold):
				min_index = min_index - 1
				if (min_index < 0):
					break
				min_transaction = tgroup[min_index]
			while (max_transaction.value <= threshold):				
				max_index = max_index + 1
				if (max_index >= len(tgroup)):
					break
				max_transaction = tgroup[max_index]
			return min_index + 1, max_index
		# not found the threshold in tgroup.
		# in this case, min_index > max_index
		elif (mid_transaction.value < threshold):
			while (mid_transaction.value < threshold):
				mid_index = mid_index + 1
				mid_transaction = tgroup[mid_index]
			return mid_index, mid_index
		# not found the threshod in tgroup and the case of mid_transaction value > threshold
		# In this case, min_index > max_index
		else:
			while (mid_transaction.value > threshold):
				mid_index = mid_index - 1
				mid_transaction = tgroup[mid_index]
			return mid_index + 1, mid_index + 1

	
	##
	# Calculate u value which measurs difference rank sum of two groups.
	# tgroup_x: test group 1. This is consisted of transactions.
	# tgroup_y: test group 2. This is consisted of transactions.
	# tgroup_x and t_group_y already sorted by transaction value.
	##
	def __uValue(self, tgroup_x, tgroup_y):
#		for t in tgroup_y:
#			print t.value,
#		print ""
	
		u_value = 0.0
		previous_u_x_min = 0 # The rank of transaction in previous search
		previous_u_x_max = 0 # The rank of transaction in previous search
		previous_value = None # The previous expression value
		left_index = 0 # The start point of searching value.
		right_index = len(tgroup_y) - 1 # The end point of searching value.
		for t_x in tgroup_x:
#			print t_x.value
			# u_x_min: rank sum of transaction which the value < t_x in tgroup_y
			# u_x_max: rank sum of transaction which the value <= t_x in tgroup_y
			u_x_min = None
			u_y_max = None
			# If t_x.value is equal to previous one, u_x_min and u_x_max are also equals.
			if (t_x.value == previous_value):
				u_x_min = previous_u_x_min
				u_x_max = previous_u_x_max
			# Caluclate u_value because tgroup_x value exists between tgroup_y range
			else:
				u_x_min, u_x_max = self.__binarySearch(t_x.value, tgroup_y, left_index, right_index)
#				for t_y in tgroup_y:
#					print t_y.value,
#				print ""
				left_index = u_x_max
			
#			print "min: " + str(u_x_min) + ", max: " + str(u_x_max)
			# Add rank of t_x to u_value
			u_value = u_value + (u_x_min+u_x_max)/2
#			print "u_value: %f" % u_value
			previous_u_x_min = u_x_min
			previous_u_x_max = u_x_max
			previous_value = t_x.value
#		print "final u_value: %f" % u_value
		return u_value
	
	##
	# This function returns mean and variance which is used in U test.
	##
	def __calStatValue(self, tgroup_x, tgroup_y):
		size_x = len(tgroup_x)
		size_y = len(tgroup_y)
		mean_u =  (size_x*size_y)/2
		var_u = size_x*size_y*(size_x + size_y + 1)/12
		return (mean_u, var_u)
	
	##
	# Calculate p-value by using Mann-Whitney U test
	# tgroup_x: test group 1. This is consisted of transactions.
	# tgroup_y: test group 2. This is consisted of transactions.
	# tgroup_x and t_group_y already sorted by transaction value.
	##
	def __uTest(self, tgroup_x, tgroup_y):
		u_value = self.__uValue(tgroup_x, tgroup_y) # u-value of two groups.
		# z value of u-value
		mean_u, var_u = self.__calStatValue(tgroup_x, tgroup_y)
		if var_u == 0:
			return 1.0, 0
		z_value = (u_value - mean_u)/math.sqrt(var_u)
#		sys.stderr.write("u: %s, mean: %s, var: %s, z_value: %s, " % (u_value, mean_u, var_u, z_value))
		
		# calculate p-value from z_value
		# this value approximation of standard normal distribution
		return self.stdNorDistribution(z_value), z_value
		
	##
	# divide transaction_list to two groups.
	# One group is transactions which included itemset, the other is not.
	# itemset: itemset. transaction_list is divided which includes this itemset or not.
	##
	def __divideGroup(self, frequent_itemset):
		in_t_list = [] # transactions which have itemset.
		out_t_list = [] # transactions does not have itemset.
		# If itemset of t contains test itemset, t puts in_t_list.
		# Else, t puts out_t_list
		for i in range(0, len(self.__transaction_list)):
			t = self.__transaction_list[i]
			if i in frequent_itemset:
				in_t_list.append(t)
			else:
				out_t_list.append(t)
#		sys.stderr.write("in_t_list: \n")
#		for t in in_t_list:
#			sys.stderr.write("%s %s\n" % (t.name, t.value) )
#		sys.stderr.write("out_t_list: \n")
#		for t in out_t_list:
#			sys.stderr.write("%s %s\n" % (t.name, t.value) )
		return in_t_list, out_t_list

	##
	# This function calculates the minimum p-value which support size is x.
	# That is, calculates MASL.
	# The z-value that minimum p-value is mean/var
	##
	def funcF(self, x):
		# calculate p-value if the group_x is consisted from max x transactions.
		min_t_list = self.__transaction_list[0:x]
		max_t_list = self.__transaction_list[x:]
#		for t in min_t_list:
#			print t.name + " ",
#			print t.itemset,
#			print " " + str(t.value)
#		print "---"
#		for t in max_t_list:
#			print t.name + " ",
#			print t.itemset,
#			print " " + str(t.value)
		mean_u, var_u = self.__calStatValue(min_t_list, max_t_list)
		min_z = mean_u/math.sqrt(var_u) # minimum z-value limited x.
#		print "  mean: " + str(mean_u) + " var: " + str(var_u) + " z-value: " + str(min_z)
		p = self.stdNorDistribution(min_z) # p-value if transaction divided into max x and other.
		if (self.alternative == 0):
			p_value = min( p * 2., 1.0 )
		return p
	
	##
	# calculate p-value by t test.
	# transaction_list: list of transactions
	# frequent_itemset: frequent transactions
	##	
	def calPValue(self, transaction_list, frequent_itemset):
#		print "--- calPValue ---"
#		print frequent_itemset
		all_size = self.__t_size # the number of all transaction (n1 + n0)
		in_t_list, out_t_list = self.__divideGroup(frequent_itemset) # dvide transactions to itemset or not.
		p_value, z_value = self.__uTest(in_t_list, out_t_list)
		if (self.alternative == 0):
			p_value = min( p_value * 2., 1.0 )
		else:
			if (z_value < 0):
				p_value = 1. - p_value
			if (self.alternative < 0):
				z_value = -z_value
		self.calTime = self.calTime + 1
		return p_value, z_value
	
	##
	# test code about stdNorDistribution
	##
#	def testStdNorDistribution(self):
#		for i in range(1, 11):
#			i = pow(10, i)
#			p = self.stdNorDistribution(i)
#			print str(i) + " " + str(p)

##
# Make mapping from column_name to column ID
##
def columnName2ID(columnid2name):
	columname2id = {}
	i = 0
	for name in columnid2name:
		columname2id[name] = i
		i = i + 1
	return columname2id

def comma2List(item_list_str, columnid2name):
	columnname2id = columnName2ID(columnid2name) # Make mapping from column_name to column ID
	item_list = []
	item_list_s = item_list_str.split(",")
	for i in item_list_s:
		columnid = columnname2id[i]
		item_list.append(columnid)
	return item_list

def run(xls_file, value_file, itemset_str_lst, delimiter, alternative):
	global readFile
	import readFile
	transaction_list, columnid2name = readFile.readFiles(xls_file, value_file, delimiter)
	
	if alternative < 0:
		global lamp
		from lamp import reverseValue
		transaction_list = reverseValue( transaction_list, "u_test" )
		
	func = FunctionOfX( transaction_list, alternative )
	colname2id_dict = readFile.colname2id(columnid2name)
	itemset = set()
	for i in itemset_str_lst:
		item_id = colname2id_dict[i]
		itemset.add(item_id + 1)

	flag_transactions_id = []
	for i in range( len(transaction_list) ):
		t = transaction_list[i]
		if len( itemset & t.itemset ) == len(itemset):
			flag_transactions_id.append( i )
	p_value, stat_score = func.calPValue(transaction_list, flag_transactions_id)
#		print i
#		print item_id
#		print columnid2name[item_id]
#	p, stat_score = func.calPValue(transaction_list, itemset)
	n = len(transaction_list)
	
	sys.stdout.write("p-value: %g (N: %s, x: %s, z-score: %f)\n" \
					 % (p_value, n, len(flag_transactions_id), stat_score))
	return (p_value, len(flag_transactions_id))

if __name__ == "__main__":
	# chech the arguments
	if (len(sys.argv) < 4):
		sys.stderr.write("usage: functions4u_test.py [item-file] [value-file] [itemset]\n")
		sys.exit()

	transaction_file = sys.argv[1]
	value_file = sys.argv[2]
	itemset_str_lst = sys.argv[3].split(',')

	# check existance of files
	if (not os.path.isfile(transaction_file)):
		print("IOError: No such file: \'" + transaction_file + "\'")
		sys.exit()
	if (not os.path.isfile(value_file)):
		print("IOError: No such file: \'" + value_file + "\'")
		sys.exit()
	
	delimiter = ','
	p_value, down_size = run(transaction_file, value_file, itemset_str_lst, delimiter)
