/*
 * Copyright (c) 2016, LAMP development team
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the LAMP development team nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL LAMP DEVELOPMENT TEAM BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "LampCore.h"

#include <boost/version.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

/**
 * Constructor
 */
LampCore::LampCore() {
	func_f = NULL;
	fre_pattern = NULL;
}

/**
 * Destructor
 */
LampCore::~LampCore() {
	clean();
}

/**
 * cleanup function
 */
void LampCore::clean() {
	if (func_f != NULL)
		delete func_f;
	func_f = NULL;
	for (enrich_t* enrich : enrich_lst) {
		delete enrich;
	}
	enrich_lst.clear();
	if (fre_pattern != NULL)
		delete fre_pattern;
	fre_pattern = NULL;
}

/**
 * Run multiple test.
 *
 * @param transaction_file The file includes associations between TFs and genes.
 *                          Each line indicates a gene.
 *                          If gene is targeted by the TF, then value is 1, otherwise 0.
 * @param flag_file Each line indicates a gene. The column1 is gene name.
 *                   If gene has the feature, the column2 is 1. The other case 0.
 * @param threshold The statistical significance threshold.
 * @param set_method The procedure name for calibration p-value (fisher/u_test/chi).
 * @param max_comb the maximal size which the largest combination size in tests set.
 * @param log_file File name for logging.
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
void LampCore::run(std::string& transaction_file, std::string& flag_file, double threshold,
        std::string& set_method, int max_comb,
        std::string& log_file, int alternative) {
	clean();
	// read 2 files and get transaction list
	std::cerr << "Read input files ..."<< std::endl;
	try {
		readFile.readFiles(transaction_file, flag_file, ',');
		// If the alternative hypothesis is 'less',
		// the positive and negative of observe values are reversed, 
		// and conduct the identical procedure to 'greater'.
		if (alternative < 0)
			reverseValue( readFile.getTransaction_list(), set_method );
		if ((int)readFile.getColumnid2name().size() <= max_comb)
			max_comb = -1;
	} catch (std::string &msg) {
		throw msg;
	} catch (...) {
		throw std::string("Error: An unexpected error occurred while trying to read input files.");
	}
	
	// run multiple test
	std::string transaction4lcm53 = transaction_file + ".4lcm53";
	// run
	try {
		FILE* fp_log = fopen(log_file.c_str(), "w");
		if (fp_log == NULL)
			throw std::string("Can't open file : " + log_file);
		boost::iostreams::stream<boost::iostreams::file_descriptor_sink> outlog;
#if (BOOST_VERSION >= 104900) 
		outlog.open(fileno(fp_log), boost::iostreams::close_handle);
#else
		outlog.open(fileno(fp_log), true);
#endif
		starttime = std::chrono::high_resolution_clock::now();
		std::cerr << "Compute the optimal correction factor ...";
		double max_lambda = maxLambda(readFile.getTransaction_list());
		fre_pattern = new LCM(max_lambda, fileno(fp_log));
		lam_star = runMultTest(readFile.getTransaction_list(), transaction4lcm53,
				threshold, set_method, max_comb, outlog, alternative,
				max_lambda);
		correction_term_time = std::chrono::high_resolution_clock::now();

		k = fre_pattern->getTotal( lam_star );
		std::cerr << " " << k << std::endl;
		std::cerr << "Compute P-values of testable combinations ..." << std::endl;
		fwerControll(max_lambda, threshold, outlog);
		finish_test_time = std::chrono::high_resolution_clock::now();
		
		outlog.close();
		
		std::cerr << "Output results ..." << std::endl;
		// If the positives and negatives are reversed, the number of positives is calculated. 
		if (( alternative < 0 ) && (std::find(BINARY_METHODS.begin(), BINARY_METHODS.end(), set_method) != BINARY_METHODS.end())) {
			for (enrich_t* l : enrich_lst) {
				l->stat_score = l->len - l->stat_score;
			}
		}
	} catch (std::string &msg) {
		throw msg;
	} catch (...) {
		throw std::string("Error: An unexpected error occurred while trying to test.");
	}
}

/**
 * Print results of multiple test.
 * @param transaction_file The file includes associations between TFs and genes.
 *                          Each line indicates a gene.
 *                          If gene is targeted by the TF, then value is 1, otherwise 0.
 * @param flag_file Each line indicates a gene. The column1 is gene name.
 *                   If gene has the feature, the column2 is 1. The other case 0.
 * @param threshold The statistical significance threshold.
 * @param set_method The procedure name for calibration p-value (fisher/u_test/chi).
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
void LampCore::print(std::string& transaction_file, std::string& flag_file, double threshold,
        std::string& set_method, int alternative) {
	try {
		// output result
		outputResult( transaction_file, flag_file, threshold, set_method,
				readFile.getColumnid2name(), readFile.getTransaction_list(), alternative );

		// output time cost
		long long ccf = std::chrono::duration_cast<std::chrono::microseconds>(correction_term_time - starttime).count();
		long long esc = std::chrono::duration_cast<std::chrono::microseconds>(finish_test_time - correction_term_time).count();
		long long tt  = std::chrono::duration_cast<std::chrono::microseconds>(finish_test_time - starttime).count();
		std::cout << "Time (sec.): Computing correction factor "
				<< boost::format("%.3f") % ((double)ccf / 1000000)
				<< ", Enumerating significant combinations "
				<< boost::format("%.3f") % ((double)esc / 1000000)
				<< ", Total "
				<< boost::format("%.3f") % ((double)tt / 1000000)
				<< std::endl;
	} catch (std::string &msg) {
		throw msg;
	} catch (...) {
		throw std::string("Error: An unexpected error occurred while trying to print.");
	}
}

/**
 * Reverse the observed values for alternative = 'less'.
 * 
 * @param transaction_list list of itemset and expression value.
 * @param set_method statistical test name. 
 */
void LampCore::reverseValue(std::vector<Transaction*>& transaction_list, const std::string& set_method ) {
	if (std::find(BINARY_METHODS.begin(), BINARY_METHODS.end(), set_method) != BINARY_METHODS.end()) {
		for (Transaction* t : transaction_list) {
			t->setValue(1.0f - t->getValue());
		}
	} else {
		for (Transaction* t : transaction_list) {
			t->setValue(0.0f - t->getValue());
		}
	}
	std::reverse(transaction_list.begin(), transaction_list.end()); 
}

/**
 * Run multiple test.
 * @param transaction_list List of itemset and expression value.
 * @param trans4lcm File name for argument of LCM program. This file is made in this method.
 * @param threshold The statistical significance threshold.
 * @param set_method The procedure name for calibration p-value (fisher/u_test).
 * @param max_comb The maximum arity limit.
 * @param outlog File name for logging.
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 * @param max_lambda Return value of maximum lambda.
 * @return 
 */
int LampCore::runMultTest(const std::vector<Transaction*>& transaction_list, std::string& trans4lcm,
		double threshold, std::string& set_method, int max_comb,
		std::ostream& outlog, int alternative, double max_lambda) {
	int lam_star = 1;
	double k = -1;
	try {
		if (set_method.compare("fisher") == 0) {
			func_f = new Functions4fisher(transaction_list, std::abs(alternative));
		}
		else if (set_method.compare("u_test") == 0) {
			func_f = new Functions4u_test(transaction_list, alternative);
		}
		else if (set_method.compare("chi") == 0) {
			func_f = new Functions4chi(transaction_list, std::abs( alternative));
		}
		else {
			throw std::string("Error: choose \"fisher\", \"chi\" or \"u_test\" by using -p option.");
		}
		double lam = max_lambda;
		
		// check a MASL of max_lambda
		double n1 = 0.0;
		if (std::find(BINARY_METHODS.begin(), BINARY_METHODS.end(), set_method) !=  BINARY_METHODS.end()) {
			n1 = func_f->sumValue();
			if (n1 < max_lambda) {
				max_lambda = int( n1 );
				lam = int( n1 );
			}
		}
		
		fre_pattern->makeFile4Lem(transaction_list, trans4lcm); // make itemset file for lcm
		
		// If Fisher's exact test or chi-square test is used for computing P-value, 
		// LCM-LAMP is run to find optimal lambda.
		if (set_method == "fisher") {
			int neg_size = func_f->getAllSize() - func_f->getN1();
			n1 = std::min( n1, (double)neg_size );
			// # of positives == # of negatives, and two.sided hypothesis test.
			if (( func_f->getN1() == neg_size ) && (alternative == 0)) {
				lam_star = depthFirst( trans4lcm, max_comb, n1, 0.5*threshold, 1 );
			} else {
				lam_star = depthFirst( trans4lcm, max_comb, n1, threshold, 1 );
			}
		}
		else if (set_method == "chi") {
			int neg_size = func_f->getAllSize() - func_f->getN1();
			n1 = std::min( n1, (double)neg_size );
			// two-sided hypothesis test
			if (alternative == 0) {
				lam_star = depthFirst( trans4lcm, max_comb, n1, 0.5*threshold, 2 );
			}
			// one-sided
			else {
				lam_star = depthFirst( trans4lcm, max_comb, n1, threshold, 2 );
			}
		}
		// If Mann-Whitney U test of Chi-square test is used,
		// LAMP ver 1. is run for computing the optimal lambda. 
		else {
			// two-sided hypothesis test
			if (alternative == 0) {
				lam_star = breadthFirst( trans4lcm, max_comb, 0.5*threshold, lam, outlog );
			}
			// one-sided hypothesis test
			else {
				lam_star = breadthFirst( trans4lcm, max_comb, threshold, lam, outlog );
			}
		}

		fre_pattern->frequentPatterns( trans4lcm, lam_star, max_comb ); // P_lambda* at line 13
		k = fre_pattern->getTotal( lam_star );
	} catch (std::string &msg) {
		throw msg;
	} catch (...) {
		throw std::string("Error: An unexpected error occurred while trying multiple test.");
	}

	// multiple test by using k and lambda_star
	outlog << "finish calculation of K: " << k << std::endl;
	// If lam_star > max_lambda, m_lambda set to max_lambda.
	// This case cause when optimal solution is found at first step.
	outlog << lam_star << std::endl;
	if (max_lambda < lam_star)
		lam_star = max_lambda;

	return lam_star;
}

/**
 * Return max lambda. That is, max size itemset.
 * @param transaction_list
 * @return 
 */
double LampCore::maxLambda(std::vector<Transaction*>& transaction_list) {
	// Count each item size
	std::map<int, int> item_sizes;
	for (Transaction* t : transaction_list) {
		for (int item : t->getItemset()) {
			// If item does not exist in item_size, then make mapping to 0
			if (item_sizes.find(item) == item_sizes.end())
				item_sizes[item] = 0;
			item_sizes[item] = item_sizes[item] + 1;
		}
	}
	
	// Get max value in item_sizes
	double max_value = 0.0f;
	for (auto itr = item_sizes.begin(); itr != item_sizes.end(); ++itr) {
		if (max_value < itr->second)
			max_value = itr->second;
	}

	// check the max lambda to the nuber of transactions
	if (( transaction_list.size() / 2.0f ) < max_value )
		max_value = transaction_list.size() / 2.0f;
	
	return max_value;
}

/**
 * Find the optimal lambda by depth first algorithm.
 * This function is called when Fisher's exact test or Chi-square is selected.
 * @param trans4lcm File name to run LCM. 
 * @param max_comb The maximum arity limit. 
 * @param n1 The number of positive samples.
 * @param threshold Significance level.
 * @param p_mode the integer. 1 -> Fisher's exact test, 2 -> chi-square test
 * @return 
 */
int LampCore::depthFirst( const std::string& trans4lcm, int max_comb,
		int n1, double threshold, int p_mode ) {
	int lam = fre_pattern->runLCMLAMP( trans4lcm, max_comb, n1, threshold, p_mode );
	fre_pattern->frequentPatterns( trans4lcm, lam, max_comb );
	return lam;
}

/**
 * Find the optimal lambda by breadth first algorithm.
 * This function is called when Mann-Whitney U- test of Chi-square test is selected as the statistical test. 
 * @param trans4lcm File name to run LCM. 
 * @param max_comb The maximum arity limit.
 * @param threshold Significance level.
 * @param lam The initializing value of lambda. 
 * @param outlog Output stream for logging.
 * @return 
 */
int LampCore::breadthFirst( const std::string& trans4lcm, int max_comb,
		double threshold, int lam, std::ostream& outlog ) {
	// solve K and lambda
	while (1 < lam) {
		outlog << "--- lambda: " << lam << " ---" << std::endl;
		// if lambda == 1, all tests which support >= 1 are tested.
		if (lam == 1) {
			fre_pattern->frequentPatterns( trans4lcm, lam, max_comb ); // line 3 of Algorithm
			break;
		}
	
		fre_pattern->frequentPatterns( trans4lcm, lam, max_comb ); // line 3 of Algorithm
		int m_lambda = fre_pattern->getTotal( lam ); // line 4 of Algorithm
		outlog << "  m_lambda: " << m_lambda << std::endl;
		
		double f_lam_1 = calBound( lam-1 ); // f(lam-1)
		outlog << "  f(" << (lam-1) << ") = " << (f_lam_1) << std::endl;
		int bottom;
		if (f_lam_1 == 0) {
			bottom = std::numeric_limits<int>::max();//.maxint;
		} else {
			bottom = (int)(threshold / f_lam_1) + 1; // bottom of line 5 of Algorithm
		}
		double f_lam = calBound( lam ); // f(lam)
		outlog << "  f(" << lam << ") = " << f_lam << std::endl;
		// If f(lambda) > f(lambda-1), raise error.
		// Because MASL f(x) is smaller if x is larger.
		if (f_lam > f_lam_1) {
			throw std::string("MASLError: f(" + std::to_string(lam) + ") = " +
					(boost::format("%.3g") % f_lam).str() +
					" is larger than f(" + std::to_string(lam - 1) + ") = " +
					(boost::format("%.3g") % f_lam_1).str());
		}
		int top;
		if (f_lam == 0) {
			top = std::numeric_limits<int>::max();
		} else {
			top = (int)(threshold / f_lam); // top of line 5 of Algorithm
		}
		outlog << "  " << bottom << " <= m_lam:" << m_lambda << " <= " << top << "?" << std::endl;
		if (bottom <= m_lambda && m_lambda <= top) { // branch on condition of line 5
			break;
		}
		outlog << "  " << m_lambda << " > " << top << "?" << std::endl;
		if (top < m_lambda) { // branch on condition of line 8
			break;
		}
		lam = lam -1;
	}
	return lam;
}

/**
 * Return the bound of given minimum support.
 * @param min_sup
 * @return 
 */
double LampCore::calBound(int min_sup ) {
	if (min_sup == 0) {
		return 1.0;
	}
	// If lower bound is not calculated, calculate the value and save to fre_pattern.
	if (1 < fre_pattern->getBound( min_sup )) {
		double bound = func_f->funcF( min_sup ); // minimum support value
		fre_pattern->setBound( min_sup, bound ); // save
	}
	return fre_pattern->getBound( min_sup );
}

/**
 * list up the combinations p_i <= alpha/k
 * @param max_lambda
 * @param threshold The statistical significance threshold.
 * @param outlog
 */
void LampCore::fwerControll(double max_lambda, double threshold, std::ostream& outlog) {
	int k = fre_pattern->getTotal( lam_star );
	int i = 0;
	int max_itemset_size = 0; // the maximum itemset size in detection of our method. This value is used for Bonferroni correction.
	for (int l = max_lambda; lam_star <= l; l--) {
		std::vector<Node::itemset_t*> item_trans_list = fre_pattern->getItemsetList( l );
		for (int j = 0; j < (int)item_trans_list.size(); j++) {
			i = i + 1;
			std::vector<int>* item_set = item_trans_list[j]->item_list;//item_set_and_size[0];
			outlog << "--- testing " << std::to_string(i) << " : set([";
			std::sort(item_set->begin(), item_set->end());
			for (int sval : *item_set) {
				outlog << sval;
				if (sval != *item_set->rbegin()) outlog << ", ";
			}
			outlog << "])";
			std::vector<int> flag_transaction_list; // transaction list which has all items in itemset.
			for (int t : *item_trans_list[j]->tran_list) {
				flag_transaction_list.push_back( t );
			}
			double stat_score;
			double p = func_f->calPValue(flag_transaction_list, stat_score);
			outlog << "p: " << std::to_string(p) << std::endl;
			if (p < (threshold / k)) {
				enrich_t* enrich = new enrich_t{ item_set, p,
						(int)item_trans_list[j]->tran_list->size(), stat_score};
				enrich_lst.push_back(enrich);
				int item_set_size = item_set->size();
				if (max_itemset_size < item_set_size)
					max_itemset_size = item_set_size;
			}
		}
	}
}

/**
 * Print result to standart output.
 * @param transaction_file The file includes associations between TFs and genes.
 * @param flag_file Each line indicates a gene. The column1 is gene name.
 * @param threshold The statistical significance threshold.
 * @param set_method The procedure name for calibration p-value (fisher/u_test/chi).
 * @param columnid2name
 * @param transaction_list
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
void LampCore::outputResult(std::string& transaction_file, std::string& flag_file,
		double threshold, std::string& set_method, std::vector<std::string*>& columnid2name,
		std::vector<Transaction*>& transaction_list,
		int alternative ) {
	int flag_size = -1;
	if (std::find(BINARY_METHODS.begin(), BINARY_METHODS.end(), set_method) != BINARY_METHODS.end()) 
		flag_size = func_f->getN1();
	// output setting
	std::cout << "# LAMP ver. " << __LAMP_VER__ << std::endl;
	std::cout << "# item-file: " << transaction_file << std::endl;
	std::cout << "# value-file: " << flag_file << std::endl;
	std::cout << "# significance-level: " << threshold << std::endl;
	std::cout << "# P-value computing procedure: " << set_method;
	if (0 < alternative) {
		std::cout << " (greater)" << std::endl;
	}
	else if (alternative < 0) {
		std::cout << " (less)" << std::endl;
	}
	else {
		std::cout << " (two.sided)" << std::endl;
	}
	std::cout << "# # of tested elements: "<< columnid2name.size() << ", # of samples: " << transaction_list.size();
	if (0 < flag_size)
		std::cout << ", # of positive samples: " << flag_size;
	std::cout << std::endl;
	std::cout << "# Adjusted significance level: " << boost::format("%.5g") % (threshold/k) << ", ";
	std::cout << "Correction factor: " << k << " (# of target rows >= " << lam_star << ")" << std::endl;
	std::cout << "# # of significant combinations: " << enrich_lst.size() << std::endl;
	// output header
	if (0 < enrich_lst.size()) {
		std::cout << "Rank\tRaw p-value\tAdjusted p-value\tCombination\tArity\t# of target rows\t";
		if (set_method == "u_test")
			std::cout << "z-score" << std::endl;
		else
			std::cout << "# of positives in the targets" << std::endl;
		std::sort(enrich_lst.begin(), enrich_lst.end(), LampCore::cmpEnrich);
		int rank = 0;
		for (enrich_t* l : enrich_lst) {
			rank = rank + 1;
			std::cout << rank << "\t" << boost::format("%.5g") % l->p << "\t" << boost::format("%.5g") % (k*l->p);
			std::string out_column = "\t";
			for (int i : *l->item_set) {
				out_column += *columnid2name[i-1] + ",";
			}
			out_column.erase( --out_column.end() );
			std::cout << out_column << "\t" << l->item_set->size() << "\t" << l->len << "\t";
			if (set_method == "u_test")
				std::cout << boost::format("%.5g") % l->stat_score << std::endl;
			else
				std::cout << boost::format("%d") % l->stat_score << std::endl;
		}
	}
}
