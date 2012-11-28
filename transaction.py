#!/usr/local/bin/python
# -*- coding: utf-8 -*-

# Define Class that indicate a transaction.
# The transaction has name and item set
# @author aika 26, June, 2011
# @edotor aika 1, Nov, 2011
#         advance for t-test. before: flag, after: value

##
# Define Class
# name: transaction name
# itemset: items that is belonged to transaction
# value: indicate this transaction related to feature.
#        if fisher's exact test, the value is 1 or 0.
#        t-test, the example is gene expression.
##
class Transaction:
	def __init__(self, name):
		self.name = name # gene name
		self.itemset = set() # item set
		self.value = None
	
	##
	# This function is used for sort transaction list in t_test
	##
	def __cmp__(self, other):
#		print str(self.value) + " " + str(other.value)
#		print cmp(self.value, other.value)
		return cmp(self.value, other.value)

	##
	# Add item to this instance
	# item: add item ID
	##
	def addItem(self, item):
		self.itemset.add(item)
