#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

import subprocess, os, time, sys
import readFile, transaction

class LCMError(Exception):
	def __init__(self, e):
		sys.stderr.write("LCMError: " + e + "\n")

class LCM():
	def __init__(self, lcm_pass):
		if lcm_pass == None:
			current_dir = os.getcwd() # curent directory
			self.__LCMPASS = current_dir + "/lcm53/lcm"
			self.__LCMNAME = "lcm53"
		else:
			self.__LCMPASS = lcm_pass
			lcm_pass_s = lcm_pass.split("/")
			self.__LCMNAME = lcm_pass_s[len(lcm_pass_s) - 1]
		# check the LCM code exists
		if not (os.path.isfile(self.__LCMPASS)):
			sys.stderr.write("Error: %s does not exist.\n" % self.__LCMPASS)
			sys.exit()
		
	
	##
	# Make file for executing lcm25 code.
	# lcm25 can be download from http://research.nii.ac.jp/~uno/codes-j.html
	# transaction_list: list of transactions
	##
	def makeFile4Lem(self, transaction_list, output_file):
		fw = open(output_file, 'w')
		for i in range(0, len(transaction_list)):
			t = transaction_list[i]
			for item in t.itemset:
				fw.write(str(item) + " ")
				#fw.write(str(item) + " ")
			if len(t.itemset) > 0:
				fw.write('\n')
		fw.close()

	##
	# Read LCM result file and return itemset list.
	##
	def readResultLCMFile(self, result_lcm_file):
		# convert output of lcm_basic to item set list
		itemset_list = []
		last_line = ""
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
				# convert output of lcm_basic to item set list
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
					itemset_list.append([itemset, s[len(s)-1], transactions])
				last_line = itemset_line[:-1]
				itemset_line = f.readline()
			# check the lcm result file finishing successfuly.
			if not last_line.endswith(")"):
				input_file = result_lcm_file.split("/")
				print input_file
				e_out = "LCM result file may be broken. Please command\n"
				e_out = e_out + "          $ rm -rf "
				for i in range(0, len(input_file) - 2):
					e_out = e_out + input_file[i] + "/"
				e_out = e_out + input_file[-2] + "*"
				e_out = e_out + "\n       Then rerun this script.\n"
				sys.stdout.write("Error: %s" % e_out)
				sys.exit()
			f.close()
		except IOError, e:
			sys.stderr.write("%s" % e)
			sys.exit()
		return itemset_list

		
	##
	# Return frequent patterns list.	
	# This method use lcm25 program.
	# input_file:
	# mim_sup: The number of minimum support. get item set that appeare abobe min_sup.
	# max_size: The number of the maximum item set size in each set.
	##
	def frequentPatterns(self, input_file, min_sup, max_size):
		out_dir = input_file + ".results." + self.__LCMNAME
		if not os.path.exists(out_dir):
			os.mkdir(out_dir)
		out_file_s = input_file.split("/")
		out_file_name = out_file_s[len(out_file_s)-1]
		out_file = out_dir + "/" + out_file_name + "." + str(min_sup)
		# If user set upper bound to item set size,
		# file name for LCM result is add "usize" + max_size
		if (max_size > 0):
			out_file = out_file + ".usize" + str(max_size) + ".all"
		itemset_list = []
		try:
			# If user does not set limit to itemse size,
			# run LCM to get closed frequent pattern.
			if (max_size < 0):
				subprocess.check_call([self.__LCMPASS, "CIf", input_file, str(min_sup), out_file], stdout=subprocess.PIPE)
				itemset_list = self.readResultLCMFile(out_file)
			# If user set limit to itemset size,
			# Run normal LCM and reduce itemset to obtain closed itemset.
			else:
				subprocess.check_call([self.__LCMPASS, "FIf", "-u", str(max_size), input_file, str(min_sup), out_file], stdout=subprocess.PIPE)
		except subprocess.CalledProcessError, (p):
			print 'subprocess.CalledProcessError: cmd:%s returncode:%s' % (p.cmd, p.returncode)
			sys.exit()
		
		# convert output of lcm_basic to item set list
		itemset_list = self.readResultLCMFile(out_file)
		return itemset_list
