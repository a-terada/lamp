#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Define fuctions that is used in multiple_test.py
# This source include calculate combination
# and function f which is defined in paper (C(n1, x)/C(n0+n1, x)).
# @author Terada, 26, June, 2011
# @editor Terada, 24, Feb, 2012,
#     If x > n1 in caluclation MASL, then raise error and exit (This case does not treat the case).

from __future__ import division
import sys
import functionsSuper as fs
import readFile

##
# Define class
# This class calculate function f that define in paper. (C(n1, x)/C(n0+n1, x))
# transaction_list: list of transactions
##
class FunctionOfX(fs.FunctionsSuper):
	def __init__(self, transaction_list):
		fs.FunctionsSuper.__init__(self)
		self.__t_size = len(transaction_list) # all transaction size
		self.__f_size = self.sumValue(transaction_list) # transaction size which have flag = 1 (n1)
		
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
		if self.__f_size > (self.__t_size/2):
			e_out = "The support size larger than half of all transaction size.\n"
			e_out = e_out + "                 This version does not treat this case."
			sys.exit()

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
			ans = self.__probabilityTable(x, 0, n1, all_x)
#			print "f(" + str(x) + ") = " + str(ans)
			return ans
		# n1 <= x <= n0
		elif x <= all_x:
			ans = self.__probabilityTable(n1, x - n1, n1, all_x)
		# x > n0:
		else:
#			print x
#			print n1
#			print all_x
			sys.stderr.write("Error: x > n1, n0. This code cannot consider this case.\n")
			sys.exit()
	
	def funcF2(self, x):
		n1 = self.__f_size
		n1_n0 = self.__t_size
#		print "x: " + str(x)
		all_x = n1_n0 - x
#		print "all transaction size: " + str(n1_n0)
#		print "n1+n0-x: " + str(all_x)
#		print "n1: " + str(n1)
		# if x > n1+n0-x, calculate x=n1+n0-x
		if x > all_x:
			ans = self.funcF(all_x)
#			print "f(" + str(x) + ") = " + str(ans)
			return ans
		
		# if x <= n1, calculate C(n1, x)/C(n1+n0, x)
		if x <= n1:
			sum_f = 1.0
			t = 0
			while t < x:
#				print str(n1-t) + "/" + str(n1_n0-t)
				sum_f = sum_f*(n1-t)/(n1_n0-t)
				t = t+1
#				print "f(" + str(x) + ") = " + str(sum_f)
			return sum_f
		# if x >= n1(else), calculate C(n0, x-n1)/C(n1+n0, x)
		else:
			n0 = n1_n0 - n1
			t = 0
			sum_f = 1.0
			# from 0 to x-n1, sum_f*(n0-t)*(x-t)/((n0+n1-t)*(x-n1-t))
			while t < (x-n1):
#				print str(n0-t) + "*" + str(x-t) + "/(" + str(n1_n0-t) + "*" + str(x-n1-t) + ")"
				sum_f = sum_f*(n0-t)*(x-t)/((n1_n0-t)*(x-n1-t))
				t = t+1
			# from x-n1 to x, sumf*(n0-t)/(n0+n1-t)
			while t < (x):
				sum_f = sum_f*(n0-t)/(n1_n0-t)
#				print str(n0-t) + "/" + str(n1_n0-t)
				t = t+1
			return sum_f

	##
	# Calculate p-value by using fisher's exact test.
	# transaction_list: List of transactions
	# flag_itemset_id: Transactions which have items
	##
	def calPValue(self, transaction_list, flag_transactions_id):
		all_size = self.__t_size # the number of all transaction (n1 + n0)
		all_flag_size = self.__f_size # the number of all flag 1 transaction (n1)
		all_not_flag_size = all_size - all_flag_size # the number of all flag 0 transactio (n0)
		# count trahsaction which contains itemset and flag is 1. (This is indicate a of paper.)
		item_flag = 0 # count size that flag = 1 and having itemset (a of paper)
		item_all = len(flag_transactions_id) # count all size that flag = 1 (x of paper)
		for i in flag_transactions_id:
			t = transaction_list[i]
			# If t flag = 1, then sum_has_flag ++.
#			print t.itemset,
#			print " flag: " + str(t.value)
#			print t.name
			item_flag = item_flag + t.value
		item_not_flag = item_all - item_flag # the number of transaction which contains itemset and flag is 0 (This is indicate b of paper)
		
#		print "a: " + str(item_flag) + ", b: " + str(item_not_flag) + ", x: " + str(item_all)
		p0 = self.__probabilityTable(item_flag, item_not_flag, all_flag_size, all_not_flag_size)
#		print "p0: " + str(p0)

		# Calculate hypergeometric distribution p
		# Change a from 0 to x, and if probability of table less than p0, then add p.
		p = 0
		a = item_flag
		while ( a <= item_all ) and ( a <= all_flag_size ):
			b = item_all - a
			pa = self.__probabilityTable(a, b, all_flag_size, all_not_flag_size)
#			print "a: p = " + str(pa)
			p = p + pa
			a = a + 1
#		print p
		return p, item_flag
	
	##
	# Calculate p-value by using fisher's exact test
	# transaction_list: list of transactions
	# itemset: testing item set
	##
	def calPValue_pre(self, transaction_list, itemset):
		print itemset
		all_size = self.__t_size # the number of all transaction (n1 + n0)
		all_flag_size = self.__f_size # the number of all flag 1 transaction (n1)
		all_not_flag_size = all_size - all_flag_size # the number of all flag 0 transactio (n0)
		# count trahsaction which contains itemset and flag is 1. (This is indicate a of paper.)
		item_flag = 0 # count size that flag = 1 and having itemset (a of paper)
		item_all = 0 # count all size that flag = 1 (x of paper)
		i = 0
		for t in transaction_list:
			# If itemset of t contains test itemset, sum_has_flag ++.
			# Moreover, t flag = 1, then sum_has_flag ++.
			if len(t.itemset & itemset) == len(itemset):
				print i
				print t.itemset,
				print " flag: " + str(t.value)
				print t.name
				item_flag = item_flag + t.value
				item_all = item_all + 1
			i = i + 1
		item_not_flag = item_all - item_flag # the number of transaction which contains itemset and flag is 0 (This is indicate b of paper)
		
		print "a: " + str(item_flag) + ", b: " + str(item_not_flag) + ", x: " + str(item_all)
		p0 = self.__probabilityTable(item_flag, item_not_flag, all_flag_size, all_not_flag_size)
#		print "p0: " + str(p0)

		# Calculate hypergeometric distribution p
		# Change a from 0 to x, and if probability of table less than p0, then add p.
		p = 0
		a = 0
		while a <= item_all:
			b = item_all - a
			pa = self.__probabilityTable(a, b, all_flag_size, all_not_flag_size)
#			print "a: p = " + str(pa)
			if pa <= p0:
				p = p + pa
			a = a + 1
		return p, item_flag
	
	##
	# Calculate probability of occurrence probability about table.
	# Return C(n1, a)*C(n0, b)/C(n1+n0, a+b)
	# a: Top left of table
	# b: Top right of table
	# n1: Sum of top and bottom left (a + c)
	# n0: Sum of top and bottom right (b + d)
	##
	def __probabilityTable(self, a, b, n1, n0):
		n1_n0 = n1 + n0 # the number of all transaction (n1 + n0)
		a_b = a + b # the number of support size (x = a + b)
		time = 0
		p = 1
		while time < a:
#			print "c(n1, a): " + str(n1-time) + "/" + str(a-time)
			p = p*(n1-time)/(a-time) # c(n1, a)
#			print "c(n1+n0, x)" + str(n1_n0-time) + "/" + str(a_b-time)
			p = p*(a_b-time)/(n1_n0-time) # c(n1+n0, x)
			time = time + 1
		time = 0
		while time < b:
#			print "c(n0, b): " + str(n0-time) + "/" + str(b-time)
			p = p*(n0-time)/(b-time) # c(n0, b)
			minus_denominator = a+time
#			print "c(n1+n0, x): " + str(n1_n0-minus_denominator) + "/" + str(a_b-minus_denominator)
			p = p*(a_b-minus_denominator)/(n1_n0-minus_denominator) # c(n1+n0, x)
			time = time + 1
#			print p
		return p	


def run(xls_file, value_file, itemset_str_lst):
	transaction_list, columnid2name, lcm2transaction_id = readFile.readFiles(xls_file, value_file)
	func = FunctionOfX(transaction_list)
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


if __name__ == "__main__":
	if (len(sys.argv) < 4):
		sys.stderr.write("Error: functions4fisher.py [xls_file] [value_file] [itemset]\n")
		sys.exit()

	xls_file = sys.argv[1]
	value_file = sys.argv[2]
	itemset_str_lst = sys.argv[3].split(',')
	run(xls_file, value_file, itemset_str_lst)
