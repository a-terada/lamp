#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Define functions which is used the each test method (fisher, t_test, u_test and so on.).
# Definition about Errors and combinations.
# @author Terada 10, Nov., 2011
import sys

class TestMethodError:
	def __init__(self, e):
		sys.stderr.write("TestMethodError: " + e + "\n")

class FunctionsSuper:
	def __init__(self): pass
	
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
