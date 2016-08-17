/*
 * Copyright (c) 2013, LAMP development team
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

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <algorithm>
#include <boost/foreach.hpp>
#include <gtest/gtest.h>
#include <boost/date_time.hpp>
#include "true_comb_list.h"
#include "LampCore.h"

using namespace std;
using namespace boost;

/** @file */
/**
 * Test case of LAMP function
 */
class Test_lamp : public testing::Test {
public:
	string RESULT_FILE; /**< output file lamp result */
	string LOG_FILE; /**<  File name for logging */

	string csv_file; /**< provides the associations between TFs and their target genes. */
	string flag_file; /**< provides the gene expressions levels. */
	string flag_less_file; /**< provides the gene expressions levels. */
	string value_file; /**< provides the gene expressions levels. */
	string value_less_file; /**< provides the gene expressions levels. */

	//# files to test # of positives = # of negatives.
	string csv_file2; /**< provides the associations between TFs and their target genes. */
	string flag_file2; /**< provides the gene expressions levels. */
	string flag_less_file2; /**< provides the gene expressions levels. */

	double sig_level; /**< threshold Significance level */

public:

	/**
	 * Function to prepare the objects for each test.
	 */
	virtual void SetUp() {
		string file_date = file_name();
		RESULT_FILE = "lamp_test_" + file_date + "_result.txt";
		LOG_FILE = "lamp_test_" + file_date + "_log.txt";
		csv_file = "sample/sample_item.csv"; 
		flag_file = "sample/sample_expression_over1.csv";
		flag_less_file = "sample/sample_expression_less1.csv";
		value_file = "sample/sample_expression_value.csv"; 
		value_less_file = "sample/sample_expression_value_rev.csv";

		csv_file2 = "sample/sample_item2.csv";
		flag_file2 = "sample/sample_expression2_over1.csv";
		flag_less_file2 = "sample/sample_expression2_less1.csv";
		sig_level = 0.05;
	}

	/**
	 * Return time stamp in the format "%Y%m%d_%H%M%S"
	 * @return string
	 */
	string file_name() {
		namespace pt = boost::posix_time;
		namespace gg = boost::gregorian;

		// date format
		auto facet = new pt::time_facet("%Y%m%d_%H%M%S");
		stringstream ss;
		ss.imbue(locale(cout.getloc(), facet));

		// current time
		auto now_time = pt::second_clock::local_time();
		ss << now_time;
		return ss.str();
	}

	/**
	 * Whether every element of array A is contained in array B.
	 * @param arr
	 * @param sub
	 * @return bool
	 */
	bool isSubset(vector<string> arr, vector<string> sub) {
		int m = arr.size();
		int n = sub.size();
		for (int i = 0; i < n; i++) {
			int j;
			for (j = 0; j < m; j++) {
				if ((!arr[j].empty() && !sub[i].empty()) && sub[i] == arr[j]) {
					break;
				}
			}
			if (j == m) {
				return false;
			}
		}
		return true;
	}

	/**
	 * compares calculated results with correct value.
	 * @param csv_file provides the associations between TFs and their target genes.
	 * @param value_file provides the gene expressions levels.
	 * @param method The procedure name for calibration p-value (fisher/u_test/chi).
	 * @param arity_lim  maximal size which the largest combination size in tests set.
	 * @param log_file log_file File name for logging.
	 * @param true_k true lamp result.
	 * @param true_lam true lamp result.
	 * @param true_comb_list right data set.
	 * @param alternative 1 -> greater, 0 -> two sided, -1 -> less
	 * @param sig_level The statistical significance threshold.
	 */
	void checkResults(string csv_file, string value_file, string method, int arity_lim, string log_file, int true_k, int true_lam, vector<True_comb_list>true_comb_list, int alternative, double sig_level) {
		LampCore lampCore;
		lampCore.run(csv_file, value_file, (double) sig_level, method, arity_lim, log_file, alternative);
		vector<LampCore::enrich_t*> enrich_lst = lampCore.getEnrich_lst();
		int k = lampCore.getK();
		int lam = lampCore.getLam_star();
		vector<string*> columnid2name = lampCore.getColumnid2name();

		cerr << "check minimum support..." << endl;
		ASSERT_EQ(lam, true_lam);
		cerr << "check correction factor..." << endl;
		ASSERT_EQ(k, true_k);
		cerr << "\n" << endl;
		cerr << "check the significance combinations...\n" << endl;

		for (LampCore::enrich_t* comb : enrich_lst) {
			const vector<int>* item_set = comb->item_set;
			vector<string> detect_set_vec;
			for (int i : *item_set) {
				detect_set_vec.push_back(*columnid2name[i - 1]);
			}
			bool flag = false;
			double true_p = -1;
			int true_support = 0;
			double true_score = -1;
			for (True_comb_list true_list : true_comb_list) {
				vector<string>true_comb_vec = true_list.true_comb;
				true_p = true_list.true_p;
				true_support = true_list.true_support;
				true_score = true_list.true_score;

				if (isSubset(detect_set_vec, true_comb_vec) && isSubset(true_comb_vec, detect_set_vec)) {
					flag = true;
					break;
				}
			}
			ASSERT_TRUE(flag);
			ASSERT_FLOAT_EQ(comb->p, true_p);
			ASSERT_EQ(comb->len, true_support);
			ASSERT_FLOAT_EQ(comb->stat_score, true_score);
		}
		ASSERT_EQ(enrich_lst.size(), true_comb_list.size());
	}
};

/**
 * test fixture for Fisher's exact test
 */
TEST_F(Test_lamp, TestFisher) {

	int true_k;
	int true_lam;

	cerr << "\n\n#######################################\n" << endl;
	cerr << "  Test LAMP using Fisher's exact test\n" << endl;
	cerr << "#######################################\n" << endl;
	cerr << "--- without arity limit (default) ---\n" << endl;
	true_k = 5;
	true_lam = 5;
	sig_level = 0.05;

	vector<True_comb_list> true_comb_list;
	True_comb_list true_list;
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.00699300699301;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, flag_file, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	cerr << "\n--- arity limit = 2 ---\n" << endl;
	true_k = 7;
	true_lam = 5;

	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2"};
	true_list.true_p = 0.00699300699301;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF1", "TF3"};
	true_list.true_p = 0.00699300699301;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2", "TF3"};
	true_list.true_p = 0.00699300699301;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, flag_file, "fisher", 2, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	cerr << "\n--- alternative=\"greater\" ---\n" << endl;
	//# # of positives != # of negatives
	true_k = 5;
	true_lam = 3;
	sig_level = 0.5;

	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.00699300699301;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2"};
	true_list.true_p = 0.034965034965;
	true_list.true_support = 6;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF3"};
	true_list.true_p = 0.034965034965;
	true_list.true_support = 6;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, flag_file, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	true_k = 5;
	true_lam = 3;
	sig_level = 0.5;
	true_comb_list.clear();
	checkResults(csv_file, flag_less_file, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	cerr << "\n--- alternative=\"two.sided\" ---\n" << endl;
	true_k = 5;
	true_lam = 3;
	sig_level = 0.5;

	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.00699300699301;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2"};
	true_list.true_p = 0.0405594405594;
	true_list.true_support = 6;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF3"};
	true_list.true_p = 0.0405594405594;
	true_list.true_support = 6;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, flag_file, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, 0, sig_level);

	true_k = 5;
	true_lam = 3;
	sig_level = 0.5;

	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.00699300699301;
	true_list.true_support = 5;
	true_list.true_score = 0;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2"};
	true_list.true_p = 0.0405594405594;
	true_list.true_support = 6;
	true_list.true_score = 1;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF3"};
	true_list.true_p = 0.0405594405594;
	true_list.true_support = 6;
	true_list.true_score = 1;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, flag_less_file, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, 0, sig_level);


	cerr << "\n--- alternative=\"less\" ---\n" << endl;
	true_k = 5;
	true_lam = 3;
	sig_level = 0.5;
	true_comb_list.clear();
	checkResults(csv_file, flag_file, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, -1, sig_level);

	true_k = 5;
	true_lam = 3;
	sig_level = 0.5;
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.00699300699301;
	true_list.true_support = 5;
	true_list.true_score = 0;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2"};
	true_list.true_p = 0.034965034965;
	true_list.true_support = 6;
	true_list.true_score = 1;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF3"};
	true_list.true_p = 0.034965034965;
	true_list.true_support = 6;
	true_list.true_score = 1;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, flag_less_file, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, -1, sig_level);

	//# # of positives == # of negatives (greater)
	cerr << "\n--- # of positives == # of negatives ---\n" << endl;
	true_k = 5;
	true_lam = 4;
	sig_level = 0.3;

	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.0104895104895;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2"};
	true_list.true_p = 0.0512820512821;
	true_list.true_support = 6;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF3"};
	true_list.true_p = 0.0512820512821;
	true_list.true_support = 6;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);
	checkResults(csv_file2, flag_file2, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	true_comb_list.clear();
	checkResults(csv_file2, flag_less_file2, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	//# # of positives == # of negatives (two.sided)
	true_k = 5;
	true_lam = 5;
	sig_level = 0.3;
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.020979020979;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);
	checkResults(csv_file2, flag_file2, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, 0, sig_level);

	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.020979020979;
	true_list.true_support = 5;
	true_list.true_score = 0;
	true_comb_list.push_back(true_list);
	checkResults(csv_file2, flag_less_file2, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, 0, sig_level);

	//# # of positives == # of negatives (less)
	true_k = 5;
	true_lam = 4;
	sig_level = 0.3;
	true_comb_list.clear();
	checkResults(csv_file2, flag_file2, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, -1, sig_level);

	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.0104895104895;
	true_list.true_support = 5;
	true_list.true_score = 0;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2"};
	true_list.true_p = 0.0512820512821;
	true_list.true_support = 6;
	true_list.true_score = 1;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF3"};
	true_list.true_p = 0.0512820512821;
	true_list.true_support = 6;
	true_list.true_score = 1;
	true_comb_list.push_back(true_list);
	checkResults(csv_file2, flag_less_file2, "fisher", -1, LOG_FILE, true_k, true_lam, true_comb_list, -1, sig_level);
}

/**
 * test fixture for Mann-Whitney U-test
 */
TEST_F(Test_lamp, testUTest) {
	int true_k;
	int true_lam;

	cerr << "\n\n#######################################\n" << endl;
	cerr << "  Test LAMP using Mann-Whitney U-test\n" << endl;
	cerr << "#######################################\n" << endl;
	cerr << "--- without arity limit ---\n" << endl;
	true_k = 5;
	true_lam = 3;

	vector<True_comb_list> true_comb_list;
	True_comb_list true_list;
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.00602414187918;
	true_list.true_support = 5;
	true_list.true_score = 2.510727;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, value_file, "u_test", -1, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	cerr << "\n--- arity limit = 2 ---\n" << endl;
	true_k = 7;
	true_lam = 3;
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2"};
	true_list.true_p = 0.00602414187918;
	true_list.true_support = 5;
	true_list.true_score = 2.510727;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF1", "TF3"};
	true_list.true_p = 0.00602414187918;
	true_list.true_support = 5;
	true_list.true_score = 2.510727;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2", "TF3"};
	true_list.true_p = 0.00602414187918;
	true_list.true_support = 5;
	true_list.true_score = 2.510727;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, value_file, "u_test", 2, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	cerr << "\n--- alternative=\"greater\" ---\n" << endl;
	//# # of positives != # of negatives
	true_k = 5;
	true_lam = 3;
	sig_level = 0.05;
	true_comb_list.clear();
	checkResults(csv_file, value_less_file, "u_test", -1, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	cerr << "\n--- alternative=\"two.sided\" ---\n" << endl;
	true_k = 5;
	true_lam = 3;
	sig_level = 0.1;
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.01204828;
	true_list.true_support = 5;
	true_list.true_score = 2.510727;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, value_file, "u_test", -1, LOG_FILE, true_k, true_lam, true_comb_list, 0, sig_level);

	true_k = 5;
	true_lam = 3;
	sig_level = 0.1;
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.01204828;
	true_list.true_support = 5;
	true_list.true_score = -2.510727;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, value_less_file, "u_test", -1, LOG_FILE, true_k, true_lam, true_comb_list, 0, sig_level);

	cerr << "\n--- alternative=\"less\" ---\n" << endl;
	true_k = 5;
	true_lam = 3;
	sig_level = 0.05;
	true_comb_list.clear();
	checkResults(csv_file, value_file, "u_test", -1, LOG_FILE, true_k, true_lam, true_comb_list, -1, sig_level);

	true_k = 5;
	true_lam = 3;
	sig_level = 0.05;
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.00602414187918;
	true_list.true_support = 5;
	true_list.true_score = -2.510727;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, value_less_file, "u_test", -1, LOG_FILE, true_k, true_lam, true_comb_list, -1, sig_level);
}

/**
 * test fixture for Chi-square test
 */
TEST_F(Test_lamp, testChiSquareTest) {
	int true_k;
	int true_lam;

	cerr << "\n\n#######################################\n" << endl;
	cerr << "  Test LAMP using the Chi-square test\n" << endl;
	cerr << "#######################################\n" << endl;

	cerr << "--- without arity limit (default) ---\n" << endl;
	true_k = 5;
	true_lam = 5;
	sig_level = 0.05;

	vector<True_comb_list> true_comb_list;
	True_comb_list true_list;
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.0086855750272;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, flag_file, "chi", -1, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	cerr << "\n--- arity limit = 2 ---\n" << endl;
	true_k = 7;
	true_lam = 5;
	sig_level = 0.1;
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2"};
	true_list.true_p = 0.0086855750272;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF1", "TF3"};
	true_list.true_p = 0.0086855750272;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2", "TF3"};
	true_list.true_p = 0.0086855750272;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, flag_file, "chi", 2, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	cerr << "\n--- alternative=\"greater\" ---\n" << endl;
	//# # of positives != # of negatives
	true_k = 5;
	true_lam = 3;
	sig_level = 0.5;
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.0086855750272;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2"};
	true_list.true_p = 0.036251012711;
	true_list.true_support = 6;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF3"};
	true_list.true_p = 0.036251012711;
	true_list.true_support = 6;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, flag_file, "chi", -1, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	true_k = 5;
	true_lam = 3;
	sig_level = 0.5;
	true_comb_list.clear();
	checkResults(csv_file, flag_less_file, "chi", -1, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	cerr << "\n--- alternative=\"two.sided\" ---\n" << endl;
	true_k = 5;
	true_lam = 4;
	sig_level = 0.5;
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.0173711500544;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2"};
	true_list.true_p = 0.0725020254219;
	true_list.true_support = 6;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF3"};
	true_list.true_p = 0.072502025419;
	true_list.true_support = 6;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, flag_file, "chi", -1, LOG_FILE, true_k, true_lam, true_comb_list, 0, sig_level);

	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.0173711500544;
	true_list.true_support = 5;
	true_list.true_score = 0;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2"};
	true_list.true_p = 0.0725020254219;
	true_list.true_support = 6;
	true_list.true_score = 1;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF3"};
	true_list.true_p = 0.072502025419;
	true_list.true_support = 6;
	true_list.true_score = 1;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, flag_less_file, "chi", -1, LOG_FILE, true_k, true_lam, true_comb_list, 0, sig_level);

	cerr << "\n--- alternative=\"less\" ---\n" << endl;
	true_k = 5;
	true_lam = 3;
	sig_level = 0.5;
	true_comb_list.clear();
	checkResults(csv_file, flag_file, "chi", -1, LOG_FILE, true_k, true_lam, true_comb_list, -1, sig_level);

	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.0086855750272;
	true_list.true_support = 5;
	true_list.true_score = 0;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2"};
	true_list.true_p = 0.036251012711;
	true_list.true_support = 6;
	true_list.true_score = 1;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF3"};
	true_list.true_p = 0.036251012711;
	true_list.true_support = 6;
	true_list.true_score = 1;
	true_comb_list.push_back(true_list);
	checkResults(csv_file, flag_less_file, "chi", -1, LOG_FILE, true_k, true_lam, true_comb_list, -1, sig_level);

	//# # of positives == # of negatives (greater)
	cerr << "\n--- # of positives == # of negatives ---\n" << endl;
	true_k = 5;
	true_lam = 4;
	sig_level = 0.3;
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.0128374712894;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2"};
	true_list.true_p = 0.05259625256;
	true_list.true_support = 6;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF3"};
	true_list.true_p = 0.05259625256;
	true_list.true_support = 6;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);
	checkResults(csv_file2, flag_file2, "chi", -1, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	true_comb_list.clear();
	checkResults(csv_file2, flag_less_file2, "chi", -1, LOG_FILE, true_k, true_lam, true_comb_list, 1, sig_level);

	//# # of positives == # of negatives (two.sided)
	true_k = 5;
	true_lam = 5;
	sig_level = 0.3;
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.0256749425788;
	true_list.true_support = 5;
	true_list.true_score = 5;
	true_comb_list.push_back(true_list);
	checkResults(csv_file2, flag_file2, "chi", -1, LOG_FILE, true_k, true_lam, true_comb_list, 0, sig_level);
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.0256749425788;
	true_list.true_support = 5;
	true_list.true_score = 0;
	true_comb_list.push_back(true_list);
	checkResults(csv_file2, flag_less_file2, "chi", -1, LOG_FILE, true_k, true_lam, true_comb_list, 0, sig_level);

	//# # of positives == # of negatives (less)
	true_k = 5;
	true_lam = 4;
	sig_level = 0.3;
	true_comb_list.clear();
	true_list.true_comb = {"TF1", "TF2", "TF3"};
	true_list.true_p = 0.0128374712894;
	true_list.true_support = 5;
	true_list.true_score = 0;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF2"};
	true_list.true_p = 0.05259625256;
	true_list.true_support = 6;
	true_list.true_score = 1;
	true_comb_list.push_back(true_list);

	true_list.true_comb = {"TF3"};
	true_list.true_p = 0.05259625256;
	true_list.true_support = 6;
	true_list.true_score = 1;
	true_comb_list.push_back(true_list);
	checkResults(csv_file2, flag_less_file2, "chi", -1, LOG_FILE, true_k, true_lam, true_comb_list, -1, sig_level);
}

/**
 * Run unit tests
 */
int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
