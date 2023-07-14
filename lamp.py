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

import sys, os.path, time, datetime, math
import transaction
import readFile
import frepattern.frequentPatterns as frequentPatterns
from optparse import OptionParser

import functions.functionsSuper as fs
import functions.functions4fisher as functions4fisher
import functions.functions4u_test as functions4u_test
import functions.functions4chi as functions4chi

#set_opts = ("fisher", "u_test", "chi") # methods which used each test

__version__ = "2.0.3"

BINARY_METHODS = tuple( [ "fisher", "chi" ] )

class MASLError(Exception):
	def __init__(self, e):
		sys.stderr.write("MASLError: " + e + "\n")

def version():
	return __version__

##
# Convert a string for the limit to a combination size to an integer.
# Without the arity limit, return -1.
# max_comb: a string represents the arity limit.
# item_size: the number of items in the dataset. 
##
def convertMaxComb( max_comb, item_size ):
	if max_comb == "all":
		return -1
	elif max_comb >= item_size: 
		return -1
	else:
		return int( max_comb )

##
# Reverse the observed values for alternative = 'less'.
# transaction_list: list of itemset and expression value.
# set_method: statistical test name. 
##
def reverseValue( transaction_list, set_method ):
	if set_method in BINARY_METHODS:
		map( lambda t: t.setValue( 1 - t.value ), transaction_list )
	else:
		map( lambda t: t.setValue( 0 - t.value ), transaction_list )
#		map( lambda t: t.setValue( math.fabs( t.value ) ), transaction_list )
	transaction_list.reverse()
	return transaction_list

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
# alternative: hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
##
def runMultTest(transaction_list, trans4lcm, threshold, set_method, lcm_path, max_comb, outlog, alternative):
	max_lambda = maxLambda(transaction_list)
	lam_star = 1; func_f = None;
	try:
		if set_method == "fisher":
			func_f = functions4fisher.FunctionOfX(transaction_list, max_lambda, abs(alternative))
		elif set_method == "u_test":
			func_f = functions4u_test.FunctionOfX(transaction_list, alternative)
		elif set_method == "chi":
			func_f = functions4chi.FunctionOfX(transaction_list, max_lambda, abs( alternative))
		else:
			sys.stderr.write("Error: choose \"fisher\", \"chi\" or \"u_test\" by using -p option.\n")
			outlog.close()
			sys.exit()
		
		lam = max_lambda
		
		# check a MASL of max_lambda
		if set_method in BINARY_METHODS: 
			n1 = func_f.sumValue(transaction_list)
			if (n1 < max_lambda):
				max_lambda = int( n1 )
				lam = int( n1 )
		
		fre_pattern = frequentPatterns.LCM(lcm_path, max_lambda, outlog)
		fre_pattern.makeFile4Lem(transaction_list, trans4lcm) # make itemset file for lcm
		
		# If Fisher's exact test or chi-square test is used for computing P-value, 
		# LCM-LAMP is run to find optimal lambda.
		if set_method == "fisher":
			neg_size = func_f.getAllSize() - func_f.getN1()
			n1 = min( n1, neg_size )
			# # of positives == # of negatives, and two.sided hypothesis test.
			if ( func_f.getN1() == neg_size ) and (alternative == 0):
				fre_pattern, lam_star = depthFirst( trans4lcm, fre_pattern, max_comb, n1, 0.5*threshold, 1 )
			else:
				fre_pattern, lam_star = depthFirst( trans4lcm, fre_pattern, max_comb, n1, threshold, 1 )
		elif set_method == "chi":
			neg_size = func_f.getAllSize() - func_f.getN1()
			n1 = min( n1, neg_size )
			# two-sided hypothesis test
			if alternative == 0:
				fre_pattern, lam_star = depthFirst( trans4lcm, fre_pattern, max_comb, n1, 0.5*threshold, 2 )
			# one-sided
			else:
				fre_pattern, lam_star = depthFirst( trans4lcm, fre_pattern, max_comb, n1, threshold, 2 )
		# If Mann-Whitney U test of Chi-square test is used,
		# LAMP ver 1. is run for computing the optimal lambda. 
		else:
			# two-sided hypothesis test
			if alternative == 0:
				fre_pattern, lam_star = breadthFirst( trans4lcm, fre_pattern, func_f, max_comb, 0.5*threshold, lam, outlog )
			# one-sided hypothesis test
			else:
				fre_pattern, lam_star = breadthFirst( trans4lcm, fre_pattern, func_f, max_comb, threshold, lam, outlog )
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

##
# Find the optimal lambda by breadth first algorithm.
# This function is called when Mann-Whitney U- test of Chi-square test is selected as the statistical test. 
# trans4lcm: File name to run LCM. 
# fre_pattern: Instance to run LCM.
# func_f: Instance to perform the statistical test. 
# max_comb: The maximum arity limit.
# threshold: Significance level.
# lam: The initializing value of lambda. 
##
def breadthFirst( trans4lcm, fre_pattern, func_f, max_comb, threshold, lam, outlog ):
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
			sys.stderr.write("MASLError: f(%s) = %.3g is larger than f(%s) = %.3g\n" \
							 % (lam, f_lam, lam-1, f_lam_1))
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
	return fre_pattern, lam

##
# Find the optimal lambda by depth first algorithm.
# This function is called when Fisher's exact test or Chi-square is selected. 
# trans4lcm: File name to run LCM. 
# fre_pattern: Instance to run LCM. 
# max_comb: The maximum arity limit. 
# n1: The number of positive samples.
# threshold: Significance level.
# p_mode: the integer. 1 -> Fisher's exact test, 2 -> chi-square test
##
def depthFirst( trans4lcm, fre_pattern, max_comb, n1, threshold, p_mode ):
	lam = fre_pattern.runLCMLAMP( trans4lcm, max_comb, n1, threshold, p_mode )
	fre_pattern.frequentPatterns( trans4lcm, lam, max_comb )
	return fre_pattern, lam


def outputResult( transaction_file, flag_file, threshold, set_method, max_comb, columnid2name, lam_star, k, \
				  enrich_lst, transaction_list, func_f, alternative ):
	flag_size = -1
	if set_method in BINARY_METHODS: 
		flag_size = func_f.getN1()
	# output setting
	sys.stdout.write("# LAMP ver. %s\n" % __version__)
	sys.stdout.write("# item-file: %s\n" % (transaction_file))
	sys.stdout.write("# value-file: %s\n" % (flag_file))
	sys.stdout.write("# significance-level: %s\n" % threshold)
	sys.stdout.write("# P-value computing procedure: %s" % set_method)
	if alternative > 0:
		sys.stdout.write(" (greater)\n")
	elif alternative < 0:
		sys.stdout.write(" (less)\n")
	else:
		sys.stdout.write(" (two.sided)\n")
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
			out_column = ""
			for i in l[0]:
				out_column = out_column + columnid2name[i-1] + ","
			sys.stdout.write("%s\t%d\t%d\t" % (out_column[:-1], len(l[0]), l[2]) )
			if set_method == "u_test":
				sys.stdout.write("%.5g\n" % l[3])
			else:
				sys.stdout.write("%d\n" % l[3])


# list up the combinations p_i <= alpha/k
def fwerControl(transaction_list, fre_pattern, lam_star, max_lambda, threshold, func_f, columnid2name, outlog):
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
			flag_transaction_list = [] # transaction list which has all items in itemset.
			for t in item_set_and_size[1]:
				flag_transaction_list.append( t )
			p, stat_score = func_f.calPValue(transaction_list, flag_transaction_list)
			outlog.write("p: " + str(p) + "\n")
			if p < (threshold/k):
				enrich_lst.append([item_set, p, len( item_set_and_size[1] ), stat_score])
				item_set_size = len(item_set)
				if ( item_set_size > max_itemset_size ):
					max_itemset_size = item_set_size
	finish_test_time = time.time()
	return ( enrich_lst, finish_test_time ) # return the number of enrich set for permutation


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

	# check the max lambda to the nuber of transactions
	if max_value > ( len(transaction_list)):
		max_value = len(transaction_list)
	
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
# max_comb: the maximal size which the largest combination size in tests set.
# delm: delimiter of transaction_file and flag_file
##
def run(transaction_file, flag_file, threshold, set_method, lcm_path, max_comb, log_file, alternative):
	# read 2 files and get transaction list
	sys.stderr.write( "Read input files ...\n" )
	transaction_list = set()
	try:
		transaction_list, columnid2name = readFile.readFiles(transaction_file, flag_file, ',')
		# If the alternative hypothesis is 'less',
		# the positive and negative of observe values are reversed, 
		# and conduct the identical procedure to 'greater'.
		if alternative < 0:
			transaction_list = reverseValue( transaction_list, set_method )
		max_comb = convertMaxComb( max_comb, len(columnid2name) )
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
		sys.stderr.write( "Compute the optimal correction factor ..." )
		fre_pattern, lam_star, max_lambda, correction_term_time, func_f \
					 = runMultTest(transaction_list, transaction4lcm53, threshold, set_method, \
								   lcm_path, max_comb, outlog, alternative)
		k = fre_pattern.getTotal( lam_star )
		sys.stderr.write( " %s\n" % k )
		sys.stderr.write( "Compute P-values of testable combinations ...\n" )
		enrich_lst, finish_test_time \
					= fwerControl(transaction_list, fre_pattern, lam_star, max_lambda, \
								   threshold, func_f, columnid2name, outlog)
		
		outlog.close()
	except IOError, e:
		outlog.close()

	sys.stderr.write( "Output results ...\n" )
	# If the positives and negatives are reversed, the number of positives is calculated. 
	if ( alternative < 0 ) and ( set_method in BINARY_METHODS ):
		for l in enrich_lst:
			l[3] = l[2] - l[3]
			
	# output result
	outputResult( transaction_file, flag_file, threshold, set_method, max_comb, \
				  columnid2name, lam_star, k, enrich_lst, transaction_list, func_f, alternative )
	# output time cost
	sys.stdout.write("Time (sec.): Computing correction factor %.3f, Enumerating significant combinations %.3f, Total %.3f\n" \
					 % (correction_term_time-starttime, finish_test_time - correction_term_time, finish_test_time - starttime))

	return enrich_lst, k, lam_star, columnid2name

if __name__ == "__main__":
	usage = "usage: %prog [options] transaction_file value_file significance_probability"
	p = OptionParser(usage = usage, version = "%s" % __version__)
	p.add_option('-p', '--pvalue', dest = "pvalue_procedure", help = "Choose the p-value calculation procedure from 'fiehser' (Fisher's exact test), 'chi' (Chi-square test) or 'u_test' (Mann-Whitney's U-test)")

	p.add_option('--lcm', dest = "lcm_path", \
				 default = os.path.dirname(os.path.abspath( __file__ )) + "/lcm53/lcm", \
				 help = "Set LCM program path if you do not use ./lcm53/lcm")

	p.add_option('--max_comb', dest = "max_comb", default = "all", \
				 help = "Set the maximum size of combination to be tested.")
	
	p.add_option('-e', dest = "log_filename", default = "", help = "The file name to output log.\n")

	p.add_option('--alternative', dest = "alternative", default = "greater", help = "Indicate which alternative hypothesis is used. Select \"greater\", \"less\" or \"two.sided\"\n, and the default is \"greater\".")

#	p.add_option('-d', dest = "delimiter", default = ",", help = "The delimiter for two input files.\n")

	opts, args = p.parse_args()
	
	# check arguments
	if len(args) != 3:
		sys.stderr.write("Error: input [target-file], [expression-file] and [significance-level].\n")
		sys.exit()

	opts.max_comb = opts.max_comb.lower()
	if not opts.max_comb == "all":
		if (opts.max_comb.isdigit()):
			opts.max_comb = int(opts.max_comb)
		else:
			sys.stderr.write("Error: max_comb must be an integer value or all.\n")
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
		if (sig_pro < 0) or (sig_pro > 1):
			raise ValueError
	except ValueError:
		sys.stderr.write("Error: significance probability must be a float value from 0.0 to 1.0.\n")
		sys.exit()

	# check the value of alternative hypothesis
	if opts.alternative == "greater":
		opts.alternative = 1
	elif opts.alternative == "less":
		opts.alternative = -1
	elif opts.alternative == "two.sided":
		opts.alternative = 0
	else:
		sys.stderr.write( "Error: \"alternative\" should be one of {\"greater\", \"less\", \"two.sided\"}\n" )
		sys.exit()
	
	# change log file
	d = datetime.datetime.today()
	log_file = "lamp_log_" + d.strftime("%Y%m%d") + "_" + d.strftime("%H%M%S") + ".txt"
	if len(opts.log_filename) > 0:
		log_file = opts.log_filename
	
	opts.delimiter = ','
	
	transaction_file = args[0]; flag_file = args[1]; threshold = float(args[2])
	enrich_lst, k, lam_star, columnid2name \
				= run(transaction_file, flag_file, threshold, opts.pvalue_procedure, \
					  opts.lcm_path, opts.max_comb, log_file, opts.alternative)
