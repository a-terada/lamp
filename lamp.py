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

# Run multiple testing correction.
# This script need transaction file, expression-file and significance-level.
# transaction-file: The file includes associations between TFs and genes.
#     Each line indicates a gene.
#     If gene is targeted by the TF, then value is 1, otherwise 0.
# expression-file: Each line indicates a gene. The column1 is gene name.
#     If gene has the feature, the column2 is 1. The other case 0.
# significance-level: The statistical significance threshold.
# @author Terada 26, June, 2011

import sys, os.path, time, datetime
import transaction
import readFile
import frepattern.frequentPatterns as frequentPatterns
from optparse import OptionParser

import functions.functionsSuper as fs
import functions.functions4fisher as functions4fisher
import functions.functions4u_test as functions4u_test
import functions.functions4chi as functions4chi

set_opts = ("fisher", "u_test", "chi") # methods which used each test

__version__ = "1.0"

class MASLError(Exception):
	def __init__(self, e):
		sys.stderr.write("MASLError: " + e + "\n")
	

##
# Return the bound of given minimum support.
##
def calBound( func_f, min_sup, fre_pattern ):
	if min_sup == 0:
		return 1.0
	# If lower bound is not calculated, calculate the value and save to fre_pattern.
	if fre_pattern.getBound( min_sup ) > 1:
		bound = func_f.funcF( min_sup ) # minimum support value
		fre_pattern.setBound( min_sup, bound ) # save
	return fre_pattern.getBound( min_sup ) 
	
##
# Run multiple test.
# transaction_list: List of itemset and expression value.
# trans4lcm: File name for argument of LCM program. This file is made in this method.
# threshold: The statistical significance threshold.
# columnid2name: Mapping between TS id to TF name.
# lcm2transaction_id: Mapping between LCM ID to transaction id.
# set_method: The procedure name for calibration p-value (fisher/u_test).
##
def runMultTest(transaction_list, trans4lcm, threshold, set_method, lcm_pass, max_comb, outlog):
	max_lambda = maxLambda(transaction_list)
	lam_star = 1; func_f = None;
	try:
		if set_method == "fisher":
			func_f = functions4fisher.FunctionOfX(transaction_list, max_lambda)
		elif set_method == "u_test":
			func_f = functions4u_test.FunctionOfX(transaction_list)
		elif set_method == "chi":
			func_f = functions4chi.FunctionOfX(transaction_list, max_lambda)
		else:
			sys.stderr.write("Error: choose \"fisher\", \"chi\" or \"u_test\".\n")
			sys.exit
				
	except fs.TestMethodError, e:
		sys.exit()
		
	try:
		lam = max_lambda
		
		# check a MASL of max_lambda
		if (set_method == 'fisher') or (set_method == 'chi'):
			n1 = func_f.sumValue(transaction_list)
			if (n1 < max_lambda):
				max_lambda = int( n1 )
				lam = int( n1 )
		
		fre_pattern = frequentPatterns.LCM(lcm_pass, max_lambda, outlog)
		fre_pattern.makeFile4Lem(transaction_list, trans4lcm) # make itemset file for lcm
		# If file for Lcm exist, comfiem that overwrite the file.
		# solve K and lambda
		while lam > 1:
			outlog.write("--- lambda: " + str(lam) + " ---\n")
			# if lambda == 1, all tests which support >= 1 are tested.
			if lam == 1:
				lam_star = lam
				fre_pattern.frequentPatterns( trans4lcm, lam, max_comb ) # line 3 of Algorithm
				k = fre_pattern.getTotal( lam )
				break

			fre_pattern.frequentPatterns( trans4lcm, lam, max_comb ) # line 3 of Algorithm
			m_lambda = fre_pattern.getTotal( lam ) # line 4 of Algorithm
			outlog.write("  m_lambda: " + str(m_lambda) + "\n")
			
			f_lam_1 = calBound( func_f, lam-1, fre_pattern ) # f(lam-1)
			outlog.write("  f(" + str(lam-1) + ") = " + str(f_lam_1) + "\n")
			if (f_lam_1 == 0):
				bottom = sys.maxint
			else:
				bottom = (threshold//f_lam_1) + 1 # bottom of line 5 of Algorithm
			f_lam = calBound( func_f, lam, fre_pattern ) # f(lam)
			outlog.write("  f(" + str(lam) + ") = " + str(f_lam) + "\n")
			# If f(lambda) > f(lambda-1), raise error.
			# Because MASL f(x) is smaller if x is larger.
			if f_lam > f_lam_1:
				e_out = "f(" + str(lam) + ") is larger than f(" + str(lam-1) + ")"
				outlog.write("MASLError: " + e_out + "\n")
				sys.exit()
			if (f_lam == 0):
				top = sys.maxint
			else:
				top = threshold//f_lam # top of line 5 of Algorithm
			outlog.write("  " + str(bottom) + " <= m_lam:" + str(m_lambda) + " <= " + str(top) + "?\n")
			if bottom <= m_lambda and m_lambda <= top: # branch on condition of line 5
				k = m_lambda
				lam_star = lam
				break
			outlog.write("  " + str(m_lambda) + " > " + str(top) + "?\n")
			if m_lambda > top: # branch on condition of line 8
				lam_star = lam
				break
			lam = lam -1
	except fs.TestMethodError, e:
		sys.exit()
	except frequentPatterns.LCMError, e:
		sys.exit()
	
	try:
		fre_pattern.frequentPatterns( trans4lcm, lam_star, max_comb ) # P_lambda* at line 13
		k = fre_pattern.getTotal( lam_star )
	except frequentPatterns.LCMError, e:
		sys.exit()

	# multiple test by using k and lambda_star
	outlog.write("finish calculation of K: %s\n" % k)
	# If lam_star > max_lambda, m_lambda set to max_lambda.
	# This case cause when optimal solution is found at first step.
	outlog.write("%s\n" % lam_star)
	if (lam_star > max_lambda):
		lam_star = max_lambda

	correction_term_time = time.time()
	return (fre_pattern, lam_star, max_lambda, correction_term_time, func_f)

def outputResult( transaction_file, flag_file, threshold, set_method, max_comb, columnid2name, lam_star, k, \
				  enrich_lst, transaction_list, func_f ):
	flag_size = -1
	if not set_method == "u_test":
		flag_size = func_f.getN1()
	# output setting
	sys.stdout.write("# LAMP ver. %s\n" % __version__)
	sys.stdout.write("# item-file: %s\n" % (transaction_file))
	sys.stdout.write("# value-file: %s\n" % (flag_file))
	sys.stdout.write("# significance-level: %s\n" % threshold)
	sys.stdout.write("# P-value computing procedure: %s\n" % set_method)
	sys.stdout.write("# # of tested elements: %d, # of samples: %d" % ( len(columnid2name), len(transaction_list) ))
	if flag_size > 0:
		sys.stdout.write(", # of positive samples: %d" % flag_size)
	sys.stdout.write("\n")
	sys.stdout.write("# Adjusted significance level: %.5g, " % (threshold/k) )
	sys.stdout.write("Correction factor: " + str(k) + " (# of target rows >= " + str(lam_star) + ")\n" )
	sys.stdout.write("# # of significant combinations: " + str(len(enrich_lst)) + "\n")
	# output header
	if len(enrich_lst) > 0:
		sys.stdout.write("Rank\tRaw p-value\tAdjusted p-value\tCombination\tArity\t# of target rows\t")
		if set_method == "u_test":
			sys.stdout.write("z-score\n")
		else:
			sys.stdout.write("# of positives in the targets\n")
		enrich_lst.sort(lambda x,y:cmp(x[1], y[1]))
		rank = 0
		for l in enrich_lst:
			rank = rank + 1
			sys.stdout.write("%d\t%.5g\t%.5g\t" % (rank, l[1], k*l[1]))
#			sys.stdout.write(str(rank) + "\t" + str(l[1]) + "\t" + str(k*l[1]) + "\t")
			out_column = ""
			for i in l[0]:
				out_column = out_column + columnid2name[i-1] + ","
			sys.stdout.write("%s\t%d\t%d\t" % (out_column[:-1], len(l[0]), l[2]) )
			if set_method == "u_test":
				sys.stdout.write("%.5g\n" % l[3])
			else:
				sys.stdout.write("%d\n" % l[3])
#			sys.stdout.write(out_column[:-1] + "\t" + str(l[2]) + "\t" + str(l[3]) + "\n")

# list up the combinations p_i <= alpha/k
def fwerControll(transaction_list, fre_pattern, lam_star, max_lambda, threshold, lcm2transaction_id, func_f, columnid2name, outlog):
	k = fre_pattern.getTotal( lam_star )
	enrich_lst = []
	i = 0
	max_itemset_size = 0 # the maximum itemset size in detection of our method. This value is used for Bonferroni correction.
 	for l in reversed( xrange( lam_star, max_lambda + 1 )):
		item_trans_list = fre_pattern.getFrequentList( l )
		for item_set_and_size in item_trans_list:
			i = i + 1
			item_set = item_set_and_size[0]
			outlog.write("--- testing " + str(i) + " : ")
			outlog.write("%s" % item_set)
#			print item_set,
			flag_transaction_list = [] # transaction list which has all items in itemset.
			for t in item_set_and_size[1]:
#				print t,
#				print " " + str(lcm2transaction_id[t]) + ", ",
				flag_transaction_list.append(lcm2transaction_id[t])
#			print " ---"
			p, stat_score = func_f.calPValue(transaction_list, flag_transaction_list)
			outlog.write("p: " + str(p) + "\n")
			if p < (threshold/k):
				enrich_lst.append([item_set, p, len( item_set_and_size[1] ), stat_score])
				item_set_size = len(item_set)
				if ( item_set_size > max_itemset_size ):
					max_itemset_size = item_set_size
#			print "p: " + str(p)+ "  ",
#			print item_set

	finish_test_time = time.time()
	return ( enrich_lst, finish_test_time ) # return the number of enrich set for permutation
	#return (len(enrich_lst), finish_test_time) # return the number of enrich set for permutation
		
##
# Return max lambda. That is, max size itemset.
##
def maxLambda(transaction_list):
	# Count each item size
	item_sizes = {}
	for t in transaction_list:
		for item in t.itemset:
			# If item does not exist in item_size, then make mapping to 0
			if not item_sizes.has_key(item):
				item_sizes[item] = 0
			item_sizes[item] = item_sizes[item] + 1
	
	# Get max value in item_sizes
	max_value = 0
	for i in item_sizes.itervalues():
		if i > max_value:
			max_value = i
			
	return max_value

##
# Run multiple test.
# itemset_file: The file includes associations between TFs and genes.
#     Each line indicates a gene.
#     If gene is targeted by the TF, then value is 1, otherwise 0.
# flag_file: Each line indicates a gene. The column1 is gene name.
#     If gene has the feature, the column2 is 1. The other case 0.
# threshold: The statistical significance threshold.
# set_method: The procedure name for calibration p-value (fisher/u_test).
# max_comb: the maxmal size which the largest combination size in tests set.
# min_p_times: the integer whether permutation test (minP) is executed.
#     When the value is over than 0, the minP is run.
# fdr_flag: A flag to determine FWER or FDR control.
# delm: delimiter of transaction_file and flag_file
##
def run(transaction_file, flag_file, threshold, set_method, lcm_pass, max_comb, log_file, delm):
	# read 2 files and get transaction list
	transaction_list = set()
	try:
		transaction_list, columnid2name, lcm2transaction_id = readFile.readFiles(transaction_file, flag_file, delm)
		if (max_comb == None):
			max_comb = -1
	except ValueError, e:
		return
	except KeyError, e:
		return
	
	# run multiple test
	transaction4lcm53 = transaction_file + ".4lcm53"
	# run
	try:
		outlog = open( log_file, 'w' )

		starttime = time.time()
		fre_pattern, lam_star, max_lambda, correction_term_time, func_f \
					 = runMultTest(transaction_list, transaction4lcm53, threshold, set_method, \
								   lcm_pass, max_comb, outlog)
		enrich_lst, finish_test_time \
					= fwerControll(transaction_list, fre_pattern, lam_star, max_lambda, \
								   threshold, lcm2transaction_id, func_f, columnid2name, outlog)
		
	except IOError, e:
		outlog.close()
	
	# output result
	k = fre_pattern.getTotal( lam_star )
	outputResult( transaction_file, flag_file, threshold, set_method, max_comb, \
				  columnid2name, lam_star, k, enrich_lst, transaction_list, func_f )
	# output time cost
	sys.stdout.write("Time (sec.): Computing correction factor %.3f, P-value %.3f, Total %.3f\n" \
					 % (correction_term_time-starttime, finish_test_time - correction_term_time, finish_test_time - starttime))

	return enrich_lst, k, columnid2name

if __name__ == "__main__":
	usage = "usage: %prog [options] transaction_file value_file significance_probability"
	p = OptionParser(usage = usage, version = "%s" % __version__)
	p.add_option('-p', '--pvalue', dest = "pvalue_procedure", help = "Choose the p-value calculation procedure from 'fiehser' (Fisher's exact test), 'chi' (Chi-square test) or 'u_test' (Mann-Whitney's U-test)")

	p.add_option('--lcm', dest = "lcm_path", default = "./lcm53/lcm", \
				 help = "Set LCM program path if you do not use ./lcm53/lcm")

	p.add_option('--max_comb', dest = "max_comb", help = "Set the maximum size of combination to be tested.")
	
	p.add_option('-e', dest = "log_filename", default = "", help = "The file name to output log.\n")

#	p.add_option('-d', dest = "delimiter", default = ",", help = "The delimiter for two input files.\n")

	opts, args = p.parse_args()
	
	# check argsuments
	if len(args) != 3:
		sys.stderr.write("Error: input [target-file], [expression-file] and [significance-level].\n")
		sys.exit()
	max_comb = None
	if (not (opts.max_comb == None)):
		if (opts.max_comb.isdigit()):
			max_comb = int(opts.max_comb)
		else:
			sys.stderr.write("Error: max_comb must be an integer value.\n")
			sys.exit()
	
	# check p-vlaue procedure
	if not opts.pvalue_procedure in set_opts:
		sys.stderr.write("Error: Choose \"fisher\" or \"u_test\" by using -p option\n")
		sys.exit()
	
	# check the file exist.
	if not os.path.isfile(args[0]):
		sys.stderr.write("IOError: No such file: \'" + args[0] + "\'\n")
		sys.exit()
	if not os.path.isfile(args[1]):
		sys.stderr.write("IOError: No such file: \'" + args[1] + "\'\n")
		sys.exit()
	try:
		sig_pro = float(args[2])
	except ValueError:
		sys.stderr.write("ArgumentsError: significance probabiliy must be an float value from 0.0 to 1.0.\n")
		sys.exit()

	if (sig_pro < 0) or (sig_pro > 1):
		sys.stderr.write("ArgumentsError: significance probabiliy must be an float value from 0.0 to 1.0.\n")
		sys.exit()

	
	# change log file
	d = datetime.datetime.today()
	log_file = "lamp_log_" + d.strftime("%Y%m%d") + "_" + d.strftime("%H%M%S") + ".txt"
	if len(opts.log_filename) > 0:
		log_file = opts.log_filename

	opts.delimiter = ','
	
	transaction_file = args[0]; flag_file = args[1]; threshold = float(args[2])
	enrich_lst, k, columnid2name \
				= run(transaction_file, flag_file, threshold, opts.pvalue_procedure, \
					  opts.lcm_path, max_comb, log_file, opts.delimiter)
