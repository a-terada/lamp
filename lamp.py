#!/usr/bin/env python

# Run multiple testing correction.
# This script need transaction file, expression-file and significance-level.
# transaction-file: The file includes associations between TFs and genes.
#     Each line indicates a gene.
#     If gene is targeted by the TF, then value is 1, otherwise 0.
# expression-file: Each line indicates a gene. The column1 is gene name.
#     If gene has the feature, the column2 is 1. The other case 0.
# significance-level: The statistical significance threshold.
# @author Terada 26, June, 2011
# @editor Terada 8, Nov, 2011
#     Add option -t u_test testing method. (fisher or t-test)

import sys, os.path, time
import readFile, transaction
import frequentPatterns
from optparse import OptionParser

import functionsSuper as fs
import functions4fisher, functions4u_test
set_opts = ("fisher", "u_test") # methods which used each test

class MASLError(Exception):
	def __init__(self, e):
		sys.stderr.write("MASLError: " + e + "\n")

##
# Run multiple test.
# transaction_list: List of itemset and expression value.
# trans4lcm: File name for argument of LCM program. This file is made in this method.
# threshold: The statistical significance threshold.
# columnid2name: Mapping between TS id to TF name.
# lcm2transaction_id: Mapping between LCM ID to transaction id.
# set_method: The procedure name for calibration p-value (fisher/u_test).
##
def executeMultTest(transaction_list, trans4lcm, threshold, columnid2name, lcm2transaction_id, set_method, lcm_pass, max_comb):
	starttime = time.clock()
	max_lambda = maxLambda(transaction_list)
	lam_star = 1
	func_f = None
	try:
		if set_method == "fisher":
			func_f = functions4fisher.FunctionOfX(transaction_list)
		else:
			func_f = functions4u_test.FunctionOfX(transaction_list)
				
	except fs.TestMethodError, e:
		sys.exit()
	try:
#		print max_lambda
		lam = max_lambda
		
		# check a MASL of max_lambda
		if (set_method == 'fisher'):
			n1 = func_f.sumValue(transaction_list)
			if (n1 < max_lambda):
				lam = int(n1)
			
		fre_pattern = frequentPatterns.LCM(lcm_pass)
		# If the file for Lcm does not exist, comfiem that overwrite the file.
		if not os.path.isfile(trans4lcm):
			fre_pattern.makeFile4Lem(transaction_list, trans4lcm) # make itemset file for lcm
		# If file for Lcm30 exist, comfiem that overwrite the file.
		# solve K and lambda
		while lam > 1:
			sys.stderr.write("--- lambda: " + str(lam) + " ---\n")
			# if lambda == 1, all tests which support >= 1 are tested.
			if lam == 1:
				lam_star = lam
				frequent_list = fre_pattern.frequentPatterns(trans4lcm, lam, max_comb) # line 3 of Algorithm
				k = len(frequent_list)
				break
			
			frequent_list = fre_pattern.frequentPatterns(trans4lcm, lam, max_comb) # line 3 of Algorithm
			m_lambda = len(frequent_list) # line 4 of Algorithm
#			print "  frequent pattern:",
#			print frequent_list
			sys.stderr.write("  m_lambda: " + str(m_lambda) + "\n")
			
			f_lam_1 = func_f.funcF(lam-1) # f(lam-1)
			sys.stderr.write("  f(" + str(lam-1) + ") = " + str(f_lam_1) + "\n")
			sys.stderr.write(str(threshold) + "//" + str(f_lam_1) + "\n")
			if (f_lam_1 == 0):
				bottom = sys.maxint
			else:
				bottom = (threshold//f_lam_1) + 1 # bottom of line 5 of Algorithm
#			print bottom
			f_lam = func_f.funcF(lam) # f(lam)
			sys.stderr.write("  f(" + str(lam) + ") = " + str(f_lam) + "\n")
			# If f(lambda) > f(lambda-1), raise error.
			# Because MASL f(x) is smaller if x is larger.
			if f_lam > f_lam_1:
				e_out = "f(" + str(lam) + ") is larger than f(" + str(lam-1) + ")"
				sys.stderr.write("MASLError: " + e_out + "\n")
				sys.exit()
			if (f_lam == 0):
				top = sys.maxint
			else:
				top = threshold//f_lam # top of line 5 of Algorithm
			sys.stderr.write("  " + str(bottom) + " <= m_lam:" + str(m_lambda) + " <= " + str(top) + "?\n")
			if bottom <= m_lambda and m_lambda <= top: # branch on condition of line 5
				k = m_lambda
				lam_star = lam
				break
			sys.stderr.write("  " + str(m_lambda) + " > " + str(top) + "?\n")
			if m_lambda > top: # branch on condition of line 8
				lam_star = lam
				break
			lam = lam -1
	except fs.TestMethodError, e:
		sys.exit()
	except frequentPatterns.LCMError, e:
		sys.exit()
	
	try:
		frequent_list = fre_pattern.frequentPatterns(trans4lcm, lam_star, max_comb) # P_lambda* at line 13
		k = len(frequent_list)
	except frequentPatterns.LCMError, e:
		sys.exit()

	# multiple test by using k and lambda_star
	sys.stderr.write("finish calculation of K: %s\n" % k)
	# If lam_star > max_lambda, m_lambda set to max_lambda.
	# This case cause when optimal solution is found at first step.
	sys.stderr.write("%s\n" % lam_star)
	if (lam_star > max_lambda):
		lam_star = max_lambda
		
	correction_term_time = time.clock()
	
	enrich_lst = []
	i = 0
	max_itemset_size = 0 # the maximum itemset size in detection of our method. This value is used for Bonferroni correction.
	for item_set_and_size in frequent_list:
		i = i + 1
		item_set = item_set_and_size[0]
		sys.stderr.write("--- testing " + str(i) + " : ")
		flag_transaction_list = [] # transaction list which has all items in itemset.
#		print item_set_and_size[2]
		for t in item_set_and_size[2]:
#			print t,
#			print " " + str(lcm2transaction_id[t]) + ", ",
			flag_transaction_list.append(lcm2transaction_id[t])
#		print " ---"
		p, stat_score = func_f.calPValue(transaction_list, flag_transaction_list)
		sys.stderr.write("p: " + str(p) + "\n")
		if p < (threshold/k):
			enrich_lst.append([item_set, p, item_set_and_size[1], stat_score])
			item_set_size = len(item_set)
			if ( item_set_size > max_itemset_size ):
				max_itemset_size = item_set_size
#			print "p: " + str(p)+ "  ",
#			print item_set

	finish_test_time = time.clock()
	
	# calculate correction term of Bonferroni
	# the correction term is defined as the sum of the combination to max itemset_size.
	bonferroni_k = 0
	for i in range(1, max_itemset_size+1):
		bonferroni_k = bonferroni_k + func_f.combination(len(columnid2name), i)
			
	sys.stdout.write("--- results ---\n")
#	print enrich_lst
	if lam == 1:
		sys.stdout.write("Warning: lamda = 1, this means all tests which # target genes >= 1 are tested.\n")
	if (len(frequent_list) < 1):
		sys.stdout.write("Warning: there is no test which satisfying # target genes >= " + str(lam_star) + ".\n")
	sys.stdout.write("Threshold: " + str(threshold/k) + ", ")
	sys.stdout.write("Correction factor: " + str(k) + " (# of target genes >= " + str(lam_star) + ")\n" )
	sys.stdout.write("# of significance: " + str(len(enrich_lst)) + "\n")
	if len(enrich_lst) > 0:
		sys.stdout.write("Raw p-value\tAdjusted p-value\tCombination\t# of target genes\tStatistic score\n")
		enrich_lst.sort(lambda x,y:cmp(x[1], y[1]))
	for l in enrich_lst:
		sys.stdout.write(str(l[1]) + "\t" + str(k*l[1]) + "\t")
		out_column = ""
		for i in l[0]:
			out_column = out_column + columnid2name[i-1] + ","
		sys.stdout.write(out_column[:-1] + "\t" + l[2][1:-1] + "\t" + str(l[3]) + "\n")
		
	# output time cost
	sys.stdout.write("Time (sec.): Correction factor %s, P-value %s, Total %s\n" % (correction_term_time-starttime, finish_test_time - correction_term_time, finish_test_time-starttime))
	return len(enrich_lst) # return the number of enrich set for permutation
		
##
# Return max lambda. That is, max size itemset.
##
def maxLambda(transaction_list):
	# Count each item size
	item_sizes = {}
	for t in transaction_list:
#		print t.itemset
		for item in t.itemset:
#			print item
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

# transaction-file: The file includes associations between TFs and genes.
#     Each line indicates a gene.
#     If gene is targeted by the TF, then value is 1, otherwise 0.
# expression-file: Each line indicates a gene. The column1 is gene name.
#     If gene has the feature, the column2 is 1. The other case 0.
# significance-level: The statistical significance threshold.
	
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
##
def run(transaction_file, flag_file, threshold, set_method, lcm_pass, max_comb):
	# read 2 files and get transaction list
	transaction_list = set()
	try:
		transaction_list, columnid2name, lcm2transaction_id = readFile.readFiles(transaction_file, flag_file)
		if (max_comb == None):
			max_comb = -1
#			max_comb = len(transaction_list)
	except ValueError, e:
		return
	except KeyError, e:
		return
	except readFile.ReadFileError, e:
		return
	
	# run multiple test
	transaction4lcm30 = transaction_file + ".4lcm53"
	executeMultTest(transaction_list, transaction4lcm30, threshold, columnid2name, lcm2transaction_id, set_method, lcm_pass, max_comb)

if __name__ == "__main__":
	usage = "usage: %prog [options] transaction_file value_file significance_probability"
	p = OptionParser(usage = usage)
	p.add_option('-p', '--pvalue', dest = "pvalue_procedure", help = "Chose the p-value calculation procedure from 'fiehser' (Fisher's exact test) or 'u_test' (Mann-Whitney's U-test)")

#	p.add_option('--lcm', dest = "lcm_pass", help = "Set LCM program pass if you do not have it the directory in multiple_test/lcm25/fim_closed")

#	p.add_option('--max_comb', dest = "max_comb", help = "Set the maxmal item set size which consider the tests set.")
	
	opts, args = p.parse_args()

	# check argsuments
	if len(args) != 3:
		sys.stderr.write("Error: input [target-file], [expression-file] and [significance-level].\n")
		sys.exit()
	max_comb = None
#	if (not (opts.max_comb == None)):
#		if (opts.max_comb.isdigit()):
#			max_comb = int(opts.max_comb)
#		else:
#			sys.stderr.write("Error: max_comb must be an integer value.\n")
#			sys.exit()
#	opts.test_method = "fisher"
	opts.lcm_pass = None

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
		if (sig_pro < 0) or (sig_pro > 1):
			sys.stderr.write("ArgumentsError: significance probabiliy must be an float value from 0.0 to 1.0.\n")
			sys.exit()
		run(args[0], args[1], float(args[2]), opts.pvalue_procedure, opts.lcm_pass, max_comb)
	except ValueError:
		sys.stderr.write("ArgumentsError: significance probabiliy must be an float value from 0.0 to 1.0.\n")
		sys.exit()
