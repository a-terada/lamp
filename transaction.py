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
