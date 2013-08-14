#!/usr/bin/env python

"""
Copyright (c) 2013, LAMP Team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the LAMP Team nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL LAMP TEAM BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

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
