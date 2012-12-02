#!/usr/local/bin/python
# -*- coding: utf-8 -*-

# get frequent pattern in transaction lists.
# @ author aika 26, June, 2011
# @ editor aika 28, Nov, 2011
#    change the default of lcm program to lcm25, and closed frequent pattern.
#    (the old version uses lcm_basic and frequent pattern.)
# @ editor aika 8, Dec, 2011
#    limit the max size of the item set.
#    In this advance, multiple_test can only calculate from 1 to max_size set.
# @ editor aika 26, June, 2012
#    fix a bug to have max limit size of the item set.
#    From this fix, this program use lcm52.

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
#			self.__LCMPASS = current_dir.rstrip("/src")
#			self.__LCMPASS = self.__LCMPASS + "/lcm53/lcm"
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
	# make file for executing lcm25 code.
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

#	def getClosedPattern(self, result_lcm_file):
#		 # all frequtne pattern which is satisfied with support size is larther than lambda
#		itemset_list_all = self.readResultLCMFile(result_lcm_file) # all frequtne pattern which is satisfied with support size is larther than lambda
#		itemset_list_closed = []
#
#		for i in itemset_list_all:
#			print i
#		
#		return itemset_list_closed

	##
	# read LCM result file and return itemset list.
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
#				e_out = result_lcm_file + " is broken\n"
#				e_out = "          Please remake the file."
				e_out = "LCM result file may be broken. Please command\n"
				e_out = e_out + "          $ rm -rf "
				for i in range(0, len(input_file) - 2):
					e_out = e_out + input_file[i] + "/"
				e_out = e_out + input_file[-2] + "*"
				e_out = e_out + "\n       Then rerun this script.\n"
				sys.stdout.write("Error: %s" % e_out)
				sys.exit()
#				raise LCMError, e_out
			f.close()
		except IOError, e:
			sys.stderr.write("%s" % e)
			sys.exit()
		return itemset_list

		
	##
	# return frequent patterns list.	
	# this method use lcm25 program.
	# input_file:
	# mim_sup: the number of minimum support. get item set that appeare abobe min_sup.
	# max_size: the number of the maximum item set size in each set.
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
		# If out_file does not exist, execute lcm program and output the file.
		if not os.path.exists(out_file):
			itemset_list = []
			try:
				# if user does not set limit to itemse size,
				# execute LCM to get closed frequent pattern.
				if (max_size < 0):
					subprocess.check_call([self.__LCMPASS, "CIf", input_file, str(min_sup), out_file], stdout=subprocess.PIPE)
					itemset_list = self.readResultLCMFile(out_file)
				# if user set limit to itemset size,
				# execute normal LCM and reduce itemset to obtain closed itemset.
				else:
					subprocess.check_call([self.__LCMPASS, "FIf", "-u", str(max_size), input_file, str(min_sup), out_file], stdout=subprocess.PIPE)
					# obtain closed itemset from all
#					itemset_list = self.getClosedPattern(out_file)
			except subprocess.CalledProcessError, (p):
				print 'subprocess.CalledProcessError: cmd:%s returncode:%s' % (p.cmd, p.returncode)
				sys.exit()
#			while not os.path.exists(out_file):
#				print out_file
#				print "wait ..."
#				time.sleep(1)
		
		# convert output of lcm_basic to item set list
		itemset_list = self.readResultLCMFile(out_file)
#		itemset_list = []
#		last_line = ""
#		for line in open(out_file, 'r'):
			# if line startswith space, this line is ignored
#			if (line.startswith(" ")):
#				continue
			# convert output of lcm_basic to item set list
#			s = line[:-1].split(' ')
			# if itemset is not empty, add itemset to itemset_list
#			if len(s) > 1:
#				itemset = set()
#				for i in range(0, len(s)-1):
#					itemset.add(int(s[i]))
#				itemset_list.append([itemset, s[len(s)-1]])
#			last_line = line[:-1]
		# check the lcm result file finishing successfuly.
#		if not last_line.endswith(")"):
#			e_out = out_file + " is broken\n"
#			e_out = e_out + "          Please remake the file."
#			raise LCMError, e_out
		return itemset_list

	##
	# return frequent patterns list.	
	# this method use lcm30 program.
	# input_file:
	# mim_sup: the number of minimum support. get item set that appeare abobe min_sup.
	##
#	def frequentPatterns2(self, input_file, min_sup):
#		p = subprocess.Popen([self.__LCMPASS, input_file, str(min_sup)], stdout=subprocess.PIPE)
#		# convert output of lcm_basic to item set list
#		itemset_list = []
#		while 1:
#			line = p.stdout.readline()
#			if not line:
#				break
#			# convert output of lcm_basic to item set list
#			s = line[:-1].split(' ')
#			print s
#			# if itemset is not empty, add itemset to itemset_list
#			if len(s) > 1:
#				itemset = set()
#				for i in range(0, len(s)-1):
#					itemset.add(int(s[i]))
#				itemset_list.append([itemset, s[len(s)-1]])
#		return itemset_list
