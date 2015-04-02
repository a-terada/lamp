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
# @author Terada 10, March, 2014

import sys, os.path, time, datetime, random, math
import transaction, readFile, lamp
import frepattern.frequentPatterns as frequentPatterns
from optparse import OptionParser

import functions.functionsSuper as fs
import functions.functions4fisher as functions4fisher
import functions.functions4u_test as functions4u_test
import functions.functions4chi as functions4chi

set_opts = ("fisher", "u_test", "chi") # methods which used each test

__version__ =  "1.0.1" + " (LAMP ver." + lamp.__version__ + ")"


def version():
	return __version__


def getValuesList( transaction_list ):
	values_list = []
	for t in transaction_list:
		values_list.append( t.value )
	return values_list

##
# Generate the permuted values dataset.
# transaction_list: the original transaction list.
# org_values_list: the values list contained the transaction_list
##
def permuteData( transaction_list, org_values_list ):
	random_index_list = range( 0, len(org_values_list) )
	random.shuffle( random_index_list )
	permute_transaction_list = []; org2shuffled_list = [-1]*len( transaction_list )
	for i in xrange( 0, len( random_index_list ) ):
		random_index = random_index_list[ i ]
		t = transaction_list[ random_index ]
		new_t = t.copy()
		new_t.setValue( org_values_list[ i ] )
		permute_transaction_list.append( new_t )
		org2shuffled_list[ t.getID() ] = i
	return permute_transaction_list, org2shuffled_list

##
# Calculate the minimum p-value in the permute_transaction_list
# transaction_list: the permuted transaction list
# trans4lcm: the filename to run LCM
# fre_pattern: the instance to obtain the frequent pattern from the min_sup
# func_f: the instance to calculate MASL and the p-value
# max_comb: the limit to maximum combination size
# org2shuffled_list: the mapping from transaction ID from the row dataset to shuffled dataset
##
def calculateMinimumPValue( permute_transaction_list, trans4lcm, fre_pattern, func_f, \
							max_comb, org2shuffled_list ):
	min_p = 1.0; min_p_pattern = None # the minimum p-value of the permuted set
	flag = True; low_sup = fre_pattern.max_support; i = 0
	freq_time = 0
	while flag:
		starttime = time.time() # time to construct apriori
		fre_pattern.frequentPatterns( trans4lcm, low_sup, max_comb ) # construct frequent patterns
		bound = fre_pattern.getBound( low_sup )
 		if bound > 1:
			bound = func_f.funcF( low_sup ) # minimum support value
			fre_pattern.setBound( low_sup, bound )
		cal_list = fre_pattern.getFrequentList( low_sup ) # Itemset calculated its P-value
		endtime = time.time()
		freq_time += endtime - starttime # time to construct apriori
		for cal_item_set, cal_transaction_list in cal_list:
			i = i + 1
#			sys.stderr.write("--- testing %s [ " % i)
#			for j in cal_item_set:
#				sys.stderr.write("%s " % j)
#			sys.stderr.write("]: ")
			flag_transaction_list = [] # transaction list which has all items in itemset.
			for t in cal_transaction_list:
				shuffled_id = org2shuffled_list[ t ]
				flag_transaction_list.append( shuffled_id )
#			sys.stderr.write("%s" % flag_transaction_list)
			p, stat_score = func_f.calPValue( permute_transaction_list, flag_transaction_list )
#			sys.stderr.write("p " + str(p) + ", stat_score %s\n" % stat_score)
			if p < min_p:
				min_p = p; min_p_pattern = cal_item_set
		
#		sys.stderr.write( "min_p: %s, low_bound: %s, min_sup: %s\n" % (min_p, bound, low_sup) )
		# If the minimum p-value is less than the lower bound of P-value, finish the calculation.
		if (min_p < bound) or (low_sup <= 1):
			flag = False
		# If the minimum p-value is over than MASL, the minimum support is small and repeat the calculation.
		else:
			low_sup = low_sup - 1
	return min_p, low_sup, freq_time

##
# Generate a probability distribution of the minimum P-value using permuted datasets 
# transaction_list: List of itemset and expression value.
# trans4lcm: File name for argument of LCM program. This file is made in this method.
# threshold: The statistical significance threshold.
# set_method: The procedure name for calibration p-value (fisher/u_test).
# lcm_path: LCM path
# max_comb: the limit to maximum combination size
# permute_num: the number of permuted dataset used in FastWY
# outlog: file object to output logs
##
def generateMinPDist(transaction_list, trans4lcm, threshold, set_method, lcm_path, \
					 max_comb, permute_num, outlog, alternative):
#	sys.stderr.write("--- original dataset ---\n")
#	for j in transaction_list:
#		j.output()
#	sys.stderr.write("------\n")

	starttime = time.time()

	# Initialize the apriori and functinos using LAMP. 
	fre_pattern, lam_star, max_lambda, correction_term_time, func_f \
				 = lamp.runMultTest( transaction_list, trans4lcm, threshold, set_method, \
									 lcm_path, max_comb, outlog, alternative )
	
	# calculate the set of minimum p-values using permuted data
	min_p_list = [] # the list stores the minimum p-values
	org_values_list = getValuesList( transaction_list ) # Raw (non-permuted) dataset
	
	# estimate the probability distribution of the minimum p-value using permuted datasets.
	for i in xrange( 0, permute_num ):
		per_start = time.time()
		permute_transaction_list, org2shuffled_list = permuteData( transaction_list, org_values_list ) # generate the permuted dataset.
#		for j in permute_transaction_list:
#			sys.stderr.write("%s %s " % (j.id, j.name))
#			sys.stderr.write("%s" % j.itemset)
#			sys.stderr.write(" %s\n" % j.value)
		func_f.calTime = 0
		min_p, low_sup, freq_time = calculateMinimumPValue( permute_transaction_list, trans4lcm, fre_pattern,
															func_f, max_comb, org2shuffled_list )
		per_time = time.time() - per_start
		if ( i == 0 ):
			per_time = time.time() - starttime
		min_p_list.append( tuple( [ min_p, low_sup, fre_pattern.getTotal( low_sup ), freq_time, per_time, func_f.calTime ] ) )
		
		outlog.write( "[permute %s] minP %s, minSupport %s, totalTest %s, freqTime %s, totalTime %s, #ofPvalue %s\n" \
						  % (i, min_p_list[i][0], min_p_list[i][1], min_p_list[i][2], \
							 min_p_list[i][3], min_p_list[i][4], min_p_list[i][5]))
		
	return min_p_list, fre_pattern, func_f

##
# Calculate the adjusted significance level
# min_p_list: the list of minimum P-values used by FastWY
# threshold: the statistical significance threshold.
# permute_num: the number of permuted dataset used in FastWY
##
def adjustedThreshold( min_p_list, threshold, permute_num ):
	# calculate the adjusted significance level
	min_p_index = max( int( math.floor(permute_num * threshold) ) - 1, 0 )
	sorted_min_p_list =  sorted( min_p_list, cmp = lambda x,y:cmp(x[0], y[0]))
	adjusted_threshold = sorted_min_p_list[ min_p_index ][0] # the adjusted significance level.
	
	if min_p_index + 1 >= len(min_p_list):
		return adjusted_threshold, sorted_min_p_list
	while adjusted_threshold == sorted_min_p_list[ min_p_index + 1][0]:
		min_p_index = min_p_index - 1
		adjusted_threshold = sorted_min_p_list[ min_p_index ][0] # the adjusted significance level.
		if min_p_index < 0:
			min_p_index = 0
			break
	
	return adjusted_threshold, sorted_min_p_list

##
# Enumerate significant combinations (P-value <= adjusted threshold)
# transaction_list: the original transaction list.
# trans4lcm: the filename to run LCM
# fre_pattern: the instance to obtain the frequent pattern from the min_sup
# func_f: the instance to calculate MASL and the p-value
# max_comb: the limit to maximum combination size
# adjusted_threshold: adjusted threshold for P-value
# outlog: file object to output logs
##
def enumerateSigComb(transaction_list, trans4lcm, fre_pattern, func_f, \
					 max_comb, adjusted_threshold, outlog):
	# test the raw (non-permuted) dataset
	i = 0; enrich_lst = []; flag = True; low_sup = fre_pattern.max_support; freq_time = 0
	start_time = time.time()
	while flag:
		freq_start_time = time.time()
		fre_pattern.frequentPatterns( trans4lcm, low_sup, max_comb ) # construct frequent patterns
		bound = fre_pattern.getBound( low_sup )
		if bound > 1:
			bound = func_f.funcF( low_sup ) # minimum support value
			fre_pattern.setBound( low_sup, bound )
		cal_list = fre_pattern.getFrequentList( low_sup ) # Itemset calculated its P-value
		freq_time = freq_time + time.time() - freq_start_time
		func_f.calTime = 0
		for item_set, item_transid_list in cal_list:
			i = i + 1
			outlog.write("--- testing %s: " % i)
			flag_transaction_list = [] # transaction list which has all items in itemset.
			for t in item_transid_list:
				flag_transaction_list.append( t )
			p, stat_score = func_f.calPValue( transaction_list, flag_transaction_list )
			outlog.write( "p " + str(p) + ", stat_score %s\n" % stat_score )
			if ( p <= adjusted_threshold ):
				enrich_lst.append([item_set, p, low_sup, stat_score])
		# If the minimum p-value is less than MASL, finish the calculation.
		if (adjusted_threshold < bound) or (low_sup <= 1):
			flag = False
		# If the minimum p-value is over than MASL, the minimum support is small and repeat the calculation.
		else:
			low_sup = low_sup - 1
		time_enumerate_total = time.time() - start_time
	return enrich_lst, freq_time, time_enumerate_total

##
# Adjusted threshold.
# Return the adjusted P-value 
# and the number of P-values up to the row P-value in the minimum P-value distribution. 
# pvalue: a raw P-value
# sorted_min_p_list: the sorted minimum P-values
# start_index: pvalue is compared with values from the start_index-th and the subsequent P-values of the sorted_min_p_list. 
##
def adjustPval( pvalue, sorted_min_p_list, start_index ):
	for i in xrange( start_index, len(sorted_min_p_list) ):
		min_p = sorted_min_p_list[i][0]
		if min_p > pvalue:
			break
	if pvalue < sorted_min_p_list[0][0]:
		adj_pvalue = -1
	elif pvalue == sorted_min_p_list[i - 1][0]:
		adj_pvalue = i / float(len( sorted_min_p_list ))
	else:
		adj_pvalue = (i + 1) / float(len( sorted_min_p_list ))
	return adj_pvalue, i

##
# Output results.
# transaction_file: The filename includes associations between TFs and genes.
# flag_file: the filename represents gene expressions.
# permute_num: the number of permuted dataset used in FastWY
# threshold: the statistical significance threshold.
# set_method: the procedure name for calibration p-value (fisher/u_test).
# max_comb: the limit to maximum combination size
# column2name: the map between column ID and column name (TF name). 
# enrich_lst: list contains significant combinations 
# adjusted_threshold: adjusted threshold for P-value
# transaction_list: the original transaction list.
# func_f: the instance to calculate MASL and the p-value
# sorted_min_p_list: the list of minimum P-values sorted used by FastWY. This list is sorted by P-values. 
##
def outputResult( transaction_file, flag_file, threshold, permute_num, set_method, max_comb, \
				  columnid2name, enrich_lst, adjusted_threshold, transaction_list, func_f, sorted_min_p_list, alternative ):
	flag_size = -1
	if not set_method == "u_test":
		flag_size = func_f.getN1()

	# output setting
	sys.stdout.write("# FastWY ver. %s\n" % __version__)
	sys.stdout.write("# item-file: %s\n" % (transaction_file))
	sys.stdout.write("# value-file: %s\n" % (flag_file))
	sys.stdout.write("# significance-level: %s, # of permuted datasets: %s\n" % (threshold, permute_num))
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
	sys.stdout.write("# Adjusted significance level: %.5g\n" % (adjusted_threshold) )
	sys.stdout.write("# # of significant combinations: " + str(len(enrich_lst)) + "\n")

	# output header
	if len(enrich_lst) > 0:
		sys.stdout.write("Rank\tRaw p-value\tAdjusted p-value\tCombination\tArity\t# of target rows\t")
		if set_method == "u_test":
			sys.stdout.write("z-score\n")
		else:
			sys.stdout.write("# of positives in the targets\n")
		enrich_lst.sort(lambda x,y:cmp(x[1], y[1]))
		rank = 0; minp_index = 0; smallest_p = 1 / float( permute_num )
		
		for l in enrich_lst:
			rank = rank + 1
			sys.stdout.write("%d\t%.5g\t" % (rank, l[1]))
			adj_pval, minp_index = adjustPval( l[1], sorted_min_p_list, minp_index )
			if adj_pval < 0:
				sys.stdout.write("< %.5g\t" % smallest_p )
			else:
				sys.stdout.write("%.5g\t" % adj_pval)
			out_column = ""
			for i in l[0]:
				out_column = out_column + columnid2name[i-1] + ","
			sys.stdout.write("%s\t%d\t%d\t" % (out_column[:-1], len(l[0]), l[2]) )
			if set_method == "u_test":
				sys.stdout.write("%.5g\n" % l[3])
			else:
				sys.stdout.write("%d\n" % l[3])

##
# Output minimum P-value distribution.
# min_p_list: the list of minimum P-values used by FastWY
##
def outputMinP( min_p_list ):
	# output minP distribution
	sys.stdout.write("--- minimum P-values ---\n")
	sys.stdout.write("[id]\tminP\tminSupport\t#ofCandComb\tfreqTime\ttotalTime\t#ofPvalueCalculation\n")
	for i in range( 0, len( min_p_list ) ):
		sys.stdout.write("[%s]\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
						 (i, min_p_list[i][0], min_p_list[i][1], min_p_list[i][2], min_p_list[i][3], min_p_list[i][4], min_p_list[i][5]))

##
# Run FastWY.
# itemset_file: The file includes associations between TFs and genes.
#     Each line indicates a gene.
#     If gene is targeted by the TF, then value is 1, otherwise 0.
# flag_file: Each line indicates a gene. The column1 is gene name.
#     If gene has the feature, the column2 is 1. The other case 0.
# threshold: The statistical significance threshold.
# min_p_times: the integer whether permutation test (minP) is executed.
#     When the value is over than 0, the minP is run.
# set_method: The procedure name for calibration p-value (fisher/u_test).
# max_comb: the maximal size which the largest combination size in tests set.
# alternative: alternative hypothesis. 1 -> greater, -1 -> less, 0 -> two.sided.
##
def run(transaction_file, flag_file, threshold, k, set_method, lcm_path, max_comb, log_file, alternative):
	# read 2 files and get transaction list
	sys.stderr.write( "Read input files ...\n" )
	transaction_list = set()
	try:
		transaction_list, columnid2name = readFile.readFiles(transaction_file, flag_file, ",")
		# If the alternative hypothesis is 'less',
		# the positive and negative of observe values are reversed, 
		# and conduct the identical procedure to 'greater'.
		if alternative < 0:
			transaction_list = lamp.reverseValue( transaction_list, set_method )
		max_comb = lamp.convertMaxComb( max_comb, len(columnid2name) )
	except ValueError, e:
		return
	except KeyError, e:
		return

	trans4lcm = transaction_file + ".4lcm53" # the filename for outputting logs 

	# run multiple test
	try:
		outlog = open( log_file, 'w' )
	except IOError, e:
		outlog.close()

	start_time = time.time()
	# generate null distribution
	sys.stderr.write( "Calculate the minimum p-value distribution using the permutation test ...\n" )
	outlog.write("Calculate the minimum p-value distribution using the permutation test ...\n")
	min_p_list, fre_pattern, func_f = \
				generateMinPDist(transaction_list, trans4lcm, threshold, set_method, \
								 lcm_path, max_comb, k, outlog, alternative)
	# adjusted significance level
	outlog.write("Adjust significance level ...\n")
	adjusted_threshold, sorted_min_p_list = adjustedThreshold( min_p_list, threshold, k )
	outlog.write("Adjusted significance level: %s\n" % adjusted_threshold)
	correction_term_time = time.time()
	# enumerate combination whose P-value up to adjusted threshold
	outlog.write("Calculate the p-values in the given data set ...\n")	
	enrich_lst, time_enumerate_freq, time_enumerate_total = \
				enumerateSigComb( transaction_list, trans4lcm, fre_pattern, func_f, \
								  max_comb, adjusted_threshold, outlog )
	
	finish_test_time = time.time()

	# output the significant combinations
	outputResult( transaction_file, flag_file, threshold, k, set_method, max_comb, columnid2name, \
				  enrich_lst, adjusted_threshold, transaction_list, func_f, sorted_min_p_list, alternative )
	
	# output time cost
	sys.stdout.write("Time (sec.): Computing correction factor %.3f, Enumerating significant combinations %.3f, Total %.3f\n" \
					 % (correction_term_time-start_time, time_enumerate_total, finish_test_time - start_time))

	# output the minimum P-values
	outputMinP( min_p_list )
	
	outlog.close()
	
	return enrich_lst, adjusted_threshold, columnid2name


if __name__ == "__main__":
	usage = "usage: %prog [options] transaction_file value_file significance_probability k"
	p = OptionParser(usage = usage, version = "%s" % __version__)
	p.add_option('-p', '--pvalue', dest = "pvalue_procedure", help = "Choose the p-value calculation procedure from 'fiehser' (Fisher's exact test), 'chi' (Chi-square test) or 'u_test' (Mann-Whitney's U-test)")

	p.add_option('--lcm', dest = "lcm_path", default = "./lcm53/lcm", \
				 help = "Set LCM program path if you do not use ./lcm53/lcm")

	p.add_option('--max_comb', dest = "max_comb", default = "all", \
				 help = "Set the maximum size of combination to be tested.")
	
	p.add_option('-e', dest = "log_filename", default = "", help = "The file name to output log.\n")

	p.add_option('--alternative', dest = "alternative", default = "greater", help = "Indicate which alternative hypothesis is used. Select \"greater\", \"less\" or \"two.sided\"\n, and the default is \"greater\".")
	
	opts, args = p.parse_args()
	
	# check arguments
	if len(args) != 4:
		sys.stderr.write("Error: input [target-file], [expression-file], [significance-level], and [k].\n")
		sys.exit()

	opts.max_comb = opts.max_comb.lower()
	if not opts.max_comb == "all":
		try:
			opts.max_comb = int( opts.max_comb )
		except ValueError:
			sys.stderr.write("Error: max_comb must be an integer value or all.\n")
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

	# check the significance level
	try:
		sig_pro = float(args[2])
		if (sig_pro < 0) or (sig_pro > 1):
			raise ValueError
	except ValueError:
		sys.stderr.write("Error: significance probability must be an float value from 0.0 to 1.0.\n")
		sys.exit()
		
	# check the number of permuted datasets used in FastWY
	k = -1
	try:
		k = int(args[3])
		if ( k < 1 ):
			raise ValueError
	except ValueError:
		sys.stderr.write("Error: k be an integer value.\n")
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
	log_file = "fastwy_log_" + d.strftime("%Y%m%d") + "_" + d.strftime("%H%M%S") + ".txt"
	if len(opts.log_filename) > 0:
		log_file = opts.log_filename
	
	opts.delimiter = ','
	
	transaction_file = args[0]; flag_file = args[1]; threshold = float(args[2])
	enrich_lst, adjusted_threshold, columnid2name \
				= run(transaction_file, flag_file, threshold, k, opts.pvalue_procedure, \
					  opts.lcm_path, opts.max_comb, log_file, opts.alternative)
