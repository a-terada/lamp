#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Node class
# This clas contains:
#   1. minimum p-value
#   2. the total of itemset that support >= this node's support
#   3. List of (itemset, transactions)

import sys

class Node():
	def __init__(self):
		self.bound = 2.0 # the value to use upper/lower bound
		self.total = -1
		self.itemset_list = []

	def setBound(self, bound):
		self.bound = bound

	def setTotal(self, total):
		self.total = total

	def addItemSet(self, item_tuple):
		self.itemset_list.append(item_tuple)

	def getItemSet(self, i):
		return self.itemset_list[i][0]
	
	def getTransactionSet(self, i):
		return self.itemset_list[i][1]
	
	def output(self):
		sys.stderr.write("lower_bound: %s, m: %s\n" % (self.bound, self.total ) )
		for i in xrange( 0, len(self.itemset_list) ):
			sys.stderr.write("   [%s] (" % i)
			itemset = self.itemset_list[i][0]; transaction_set = self.itemset_list[i][1]
			# Print itemset
			for item in itemset:
				sys.stderr.write(" %s" % item)
			sys.stderr.write(" ): ")
			# Print transaction_set
			for t in transaction_set:
				sys.stderr.write(" %s" % t)
			sys.stderr.write("\n")
