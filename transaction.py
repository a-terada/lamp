#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Define Class that indicate a transaction.
# The transaction has name and item set
# @author aika 26, June, 2011

##
# Define Class
# name: Transaction name (Gene name)
# itemset: Items that is belonged to transaction (associated TFs set)
# value: Indicate this transaction related to feature.
#        If fisher's exact test, the value is 1 or 0.
#        If Mann-Whitney's u-test, the value takes any value.
##
class Transaction:
	def __init__(self, name):
		self.name = name # gene name
		self.itemset = set() # item set
		self.value = None
	
	##
	# This function is used for sort transaction list in u_test
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
