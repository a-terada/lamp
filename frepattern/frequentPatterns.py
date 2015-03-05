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

# Get frequent pattern in transaction lists.
# @ author Teradas 26, June, 2011
# @ editor Terada 28, Nov, 2011
#     Change the default of lcm program to lcm25, and closed frequent pattern.
#     (The old version uses lcm_basic and frequent pattern.)
# @ editor Terada 8, Dec, 2011
#     Limit the max size of the item set.
#     In this advance, multiple_test can only calculate from 1 to max_size set.
# @ editor Terada 26, June, 2012
#     Fix a bug to have max limit size of the item set.
#     From this fix, this program use lcm52.
# @ editor Terada 10, Jan. 2013
#     

import subprocess, os, time, sys
import nodeClass

class LCMError(Exception):
	def __init__(self, e):
		sys.stderr.write("LCMError: " + e + "\n")


class LCM():
	def __init__(self, lcm_path, max_support, outlog):
		self.frequent_list = [] # index: max-support - support, In the list: Instance of Node class
		self.max_support = max_support # maximum value of the minimum support.
		self.constructed_index = -1 # frequent_list is constructed that its index is less than the value.
		self.outlog = outlog
	
		# Initialize the frequent_list.
		for i in xrange(0, self.max_support):
			self.frequent_list.append( nodeClass.Node() )
		
		# set LCM code path
		if lcm_path == None:
			current_dir = os.getcwd() # curent directory
			self.__LCMPATH = current_dir + "/lcm53/lcm"
			self.__LCMNAME = "lcm"
		else:
			self.__LCMPATH = lcm_path
			lcm_path_s = lcm_path.split("/")
			self.__LCMNAME = lcm_path_s[len(lcm_path_s) - 1]
		# check the LCM code exists
		if not (os.path.isfile(self.__LCMPATH)):
			sys.stderr.write("Please set lcm binary path with --lcm option\n")
			sys.exit()
			
	def getIndex(self, support):
		return self.max_support - support
	
	##
	# Return the itemset size that support >= min_sup
	##
	def getTotal(self, min_sup):
		return self.frequent_list[ self.getIndex(min_sup) ].total
	
	##
	# Return the bound
	##
	def getBound(self, min_sup):
		return self.frequent_list[ self.getIndex(min_sup) ].bound

	def setBound(self, min_sup, bound):
		node = self.frequent_list[ self.getIndex(min_sup) ]
		node.setBound( bound )

	def getFrequentList(self, support):
		return self.frequent_list[ self.getIndex(support) ].itemset_list
	
	##
	# Make file for running LCM. 
	# LCM was downloaded from http://research.nii.ac.jp/~uno/codes-j.html
	# transaction_list: list of transactions
	##
	def makeFile4Lem(self, transaction_list, output_file):
		fw = open(output_file, 'w')
		for t in transaction_list:
			for item in t.itemset:
				fw.write(str(item) + " ")
			fw.write('\n')
		fw.close()
	
	
	##
	# Read LCM result file and return itemset list.
	# result_lcm_file: the result filename of running LCM
	# low_sup: the minimum support when LCM has run
	# upper_sup: the maximum support when LCM has run (the frequent_list re-constructed between min_sup to max_sup)
	##
	def readResultLCMFile(self, result_lcm_file, low_sup, upper_sup):
		# Initialize re-constructed nodes
		for i in xrange( low_sup, upper_sup + 1 ):
			node = nodeClass.Node()
			self.frequent_list[ self.getIndex( i ) ] = node
		
		# convert output of lcm_basic to item set list
		try:
			f = open(result_lcm_file, 'r')
			itemset_line = f.readline()
			transactions_line = ""
			while itemset_line:
				transactions_line = f.readline()
				# if line startswith space, this line is ignored
				if (itemset_line.startswith(" ")):
					itemset_line = f.readline()
					continue
				# convert output of lcm to item set list
				s = itemset_line[:-1].split(' ')
				transactions_line = transactions_line[1:]
				transactions = transactions_line[:-1].split(' ')
				# if itemset is not empty, add itemset to itemset_list
				if len(s) > 1:
					itemset = set()
					for i in range(0, len(s)-1):
						itemset.add(int(s[i]))
					if transactions_line[:-1] == "":
						transactions = []
					else:
						for i in range(0, len(transactions)):
							transactions[i] = int(transactions[i])
					support = len( transactions )
					node = None; node_index = self.getIndex( support )
					if ( node_index < 0 ):
						node = self.frequent_list[ 0 ]
					else:
						node = self.frequent_list[ self.getIndex( support ) ]
					node.addItemSet( tuple([itemset, transactions]) )
				itemset_line = f.readline()
			f.close()
		except IOError, e:
			sys.stderr.write("%s" % e)
			sys.exit()
	
		
	##
	# Construct frequent patterns list.	This method use lcm53 program.
	# input_file:
	# low_sup: The number of minimum support. get item set that appeare abobe min_sup.
	# arity_limit: The limit to appriori depth.
	##
	def frequentPatterns(self, input_file, low_sup, arity_limit):
		# If frequent pattern has already serched, then return.
		if self.getIndex( low_sup ) <= self.constructed_index:
			return
		
		out_dir = input_file + ".results." + self.__LCMNAME
		if not os.path.exists(out_dir):
			os.mkdir(out_dir)
		out_file_s = input_file.split("/")
		out_file_name = out_file_s[len(out_file_s)-1]
		out_file_pre = out_dir + "/" + out_file_name

		upper_sup = self.max_support - self.constructed_index - 1
		# Run LCM
		try:
			# If the arity size is not limited, run LCM to get closed frequent pattern.
			if ( arity_limit < 0 and self.constructed_index > -1 ):
				out_file = out_file_pre + ".lowsup" + str( low_sup ) + ".upsup" + str( upper_sup ) + ".closed"
				subprocess.check_call([self.__LCMPATH, "CIf", "-U", str(upper_sup), \
									   input_file, str(low_sup), out_file], \
									  stdout=self.outlog, stderr = self.outlog)
			elif ( arity_limit < 0 and self.constructed_index == -1 ):
				out_file = out_file_pre + ".lowsup" + str( low_sup ) + ".closed"
				subprocess.check_call([self.__LCMPATH, "CIf", input_file, str(low_sup), out_file], \
									  stdout=self.outlog, stderr = self.outlog)
			elif ( arity_limit >= 0 and self.constructed_index > -1):
				out_file = out_file_pre + ".lowsup" + str( low_sup ) + ".upsup" + str( upper_sup ) \
						   + ".aritylim" + str( arity_limit )
				subprocess.check_call([self.__LCMPATH, "FIf", "-U", str(upper_sup), "-u", \
									   str( arity_limit ), input_file, str(low_sup), out_file], \
									  stdout=self.outlog, stderr = self.outlog)
			else:
				out_file = out_file_pre + ".lowsup" + str( low_sup ) +".aritylim" + str( arity_limit )
				subprocess.check_call([self.__LCMPATH, "FIf", "-u", str( arity_limit ), \
									   input_file, str(low_sup), out_file], \
									  stdout=self.outlog, stderr = self.outlog)
		except subprocess.CalledProcessError, (p):
			sys.stderr.write('subprocess.CalledProcessError: cmd:%s returncode:%s\n' % (p.cmd, p.returncode) )
			sys.exit()
		# Read the file of LCM result
		self.readResultLCMFile( out_file, low_sup, upper_sup )

		# Update the total number of transactions
		total = 0
		if (low_sup < self.max_support) and (self.constructed_index > -1):
			total = self.frequent_list[ self.getIndex( low_sup ) - 1 ].total
		for i in xrange( upper_sup, low_sup - 1, -1 ):
			node = self.frequent_list[ self.getIndex( i ) ]
			total = total + len( node.itemset_list )
			node.total = total
		
		self.constructed_index = self.getIndex( low_sup )
	
		
	##
	# Run LCM-LAMP and return the optimal minimum support. 
	# input_file: filename for LCM.
	# arity_limit: limit for appriori depth.
	# n1: the number of positive samples.
	# sig_level: the significance level.
	# p_mode: the integer that indicates the kind of statistical test.
	#         1 -> Fisher's exact test,  2 -> chi-square test
	##
	def runLCMLAMP( self, input_file, arity_limit, n1, sig_level, p_mode ):
		out_dir = input_file + ".results." + self.__LCMNAME
		if not os.path.exists(out_dir):
			os.mkdir(out_dir)
		out_file_s = input_file.split("/")
		out_file_name = out_file_s[len(out_file_s)-1]
		out_file_pre = out_dir + "/" + out_file_name
		out_file = ""
		
		if ( arity_limit < 0 ):
			out_file = out_file_pre + ".lcmlamp.closed"
			outlog_lcmlamp = open( out_file, 'w' )
			subprocess.check_call( [self.__LCMPATH, "C", "-LAMP", str(n1), "-LAMP_P", str(p_mode), \
									input_file, str(sig_level)], \
								   stdout=outlog_lcmlamp, stderr=outlog_lcmlamp )
			outlog_lcmlamp.close()
		else:
			out_file = out_file_pre + ".lcmlamp.aritylim" + str( arity_limit )
			outlog_lcmlamp = open( out_file, 'w' )
			subprocess.check_call( [self.__LCMPATH, "F", "-LAMP", str(n1), "-u", str(arity_limit), \
									input_file, str(sig_level)], \
								   stdout=outlog_lcmlamp, stderr=outlog_lcmlamp )
			outlog_lcmlamp.close()
		
		fr = open( out_file, 'r' ); line = ""; lam = -1
		for line in fr:
			if line.startswith( "frq= " ):
				s = line.split( ' ' )
				lam = int( s[1] )
		fr.close()
		lam = lam - 1
		return lam
			
