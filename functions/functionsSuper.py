#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Define functions which is used the each test method (fisher, t_test, u_test and so on.).
# Definition about Errors and combinations.
# @author Terada 10, Nov., 2011
import sys, math

class TestMethodError:
	def __init__(self, e):
		sys.stderr.write("TestMethodError: " + e + "\n")

class FunctionsSuper:
	def __init__(self):
		self.__range20_1 = range(1, 21)
		self.__range20_1.reverse() # the integer list from 20 to 1. this is used by standard normal probability
	
	##
	# Calculate combination C(n, m)
	# return the number of pattern that select m things from n things.
	# This function use for correcting by Bonferroni method and fisher's exact test
	##
	def combination(self, n, m):
#		print "calcluate combination ..."
#		print str(n) + 'C' + str(m)
		# If (n-m) < m, calculate C(n, n-m) that same as C(n, m)
		if (n-m) < m:
			return self.combination(n, n-m)
		t = m
		c_sum = 1.0
		while t > 0:
			c_sum = c_sum*(n-m+t)/t
			t = t-1
		return c_sum

	##
	# Divide transaction_list to two groups.
	# One group is transactions which included itemset, the other is not.
	# itemset: itemset. transaction_list is divided which includes this itemset or not.
	##
	def __divideGroup(self, itemset):
		in_t_list = [] # transactions which have itemset.
		out_t_list = [] # transactions does not have itemset.
		# If itemset of t contains test itemset, t puts in_t_list.
		# Else, t puts out_t_list
		for t in self.transaction_list:
			if len(t.itemset & itemset) == len(itemset):
				print t.itemset,
				print "value: " + str(t.value)
				in_t_list.append(t)
			else:
				out_t_list.append(t)
		return in_t_list, out_t_list

	def sumValue(self, transaction_list):
#		print "count flag size ..."
		f_size = 0
		for t in transaction_list:
			f_size = f_size + t.value
#		print "finish: " + str(f_size)
		return f_size

	##
	# Calculate probability of standard normal distribution.
	# this function returns the probability of one-sided test.
	# x: 
	##
	def stdNorDistribution(self, x):
		pi2 = 0.398942280401432677940
		is_value = -1
		y = abs(x)
		c = y*y
		p = 0.0
		z = math.exp(-c*0.5)*pi2
		if (y < 2.5):
			for i in self.__range20_1:
				p = i*c/(i*2+1+is_value*p)
				is_value = -is_value
			p = 0.5-z*y/(1.0-p)
		else:
			for i in self.__range20_1:
				p = i/(y+p)
			p = z/(y+p)
#		p = 2 * p # double p-value because returens about two-sided test.
#		print str(x) + " " + str(p)
		return p

	##
	# Make the contingency table.
	# This function is used by fisher and chi-square test.
	##
	def contingencyTable( self, transaction_list, flag_transactions_id, total, total_col1 ):
		ovalues = [ [0, 0], [0, 0] ]
		total_col2 = total - total_col1 # the number of all flag 0 transactio (n0)
		# count trahsaction which contains itemset and flag is 1. (This is indicate a of paper.)
		total_row1 = len(flag_transactions_id) # count all size that flag = 1 (x of paper)
		for i in flag_transactions_id:
			t = transaction_list[i]
			# If t flag = 1, then sum_has_flag ++.
			ovalues[0][0] = ovalues[0][0] + t.value
		ovalues[0][0] = int(ovalues[0][0])
		ovalues[0][1] = total_row1 - ovalues[0][0] # the number of transaction which contains itemset and flag is 0 (This is indicate b of paper)
		ovalues[1][0] = total_col1 - ovalues[0][0]
		ovalues[1][1] = total_col2 - ovalues[0][1]
		return ovalues
