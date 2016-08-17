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
/* 
 * File:   LCM.cpp
 * Author: miura
 * 
 * Get frequent pattern in transaction lists.
 * 
 * Created on 2016/01/13, 14:46
 */

#include <ostream>

#include "LCM.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

const int LCM::FILE_NAME_BUF = 1024;

/**
 * Constructor
 * @param max_support maximum value of the minimum support.
 * @param fd_outlog file descriptor of log file.
 */
LCM::LCM(int max_support, int fd_outlog) : fd_outlog(fd_outlog) {
	this->max_support = max_support; // maximum value of the minimum support.
	constructed_index = -1; // frequent_list is constructed that its index is less than the value.
	
	// Initialize the frequent_list.
	// index: max-support - support, In the list: Instance of Node class
	for (int i = 0; i < max_support; i++) {
		frequent_list.push_back( new Node() );
	}
}

/**
 * Destructor
 */
LCM::~LCM() {
	for (Node* n : frequent_list) {
		delete n;
	}
}

/**
 * Make file for running LCM.
 * @param transaction_list list of transactions
 * @param output_file File name of output.
 */
void LCM::makeFile4Lem(const std::vector<Transaction*>& transaction_list, std::string& output_file) {
	try {
		std::ofstream fw(output_file);
		if (!fw.is_open())
			throw std::string("Can't open file. " + output_file);
		for (Transaction* t : transaction_list) {
			for (int item : t->getItemset()) {
				fw << item << " ";
			}
			fw << std::endl;
		}
		fw.close();
	} catch (std::string &msg) {
		throw msg;
	} catch (...) {
		throw std::string("Error: An unexpected error occurred while trying to make File  for Lem.");
	}
}

/**
 * Run LCM-LAMP and return the optimal minimum support. 
 * @param input_file filename for LCM.
 * @param arity_limit limit for appriori depth.
 * @param n1  the number of positive samples.
 * @param sig_level the significance level.
 * @param p_mode the integer that indicates the kind of statistical test.
 *                1 -> Fisher's exact test,  2 -> chi-square test
 * @return 
 */
int LCM::runLCMLAMP( const std::string& input_file, int arity_limit, int n1, double sig_level, int p_mode ) {
	std::string out_dir(input_file + ".results.lcm");
	if (!boost::filesystem::exists(out_dir)) {
		boost::filesystem::create_directory(out_dir);
	}
	std::list<std::string> out_file_s;
	boost::split(out_file_s, input_file, boost::is_any_of("/"));
	std::string out_file_name = out_file_s.back();
	std::string out_file_pre = out_dir + "/" + out_file_name;
	char input_file_chr[FILE_NAME_BUF];
	strncpy(input_file_chr, input_file.c_str(), FILE_NAME_BUF-1);
	
	int lam = -1;
	std::cout << std::flush;
	std::cerr << std::flush;
	FILE* fp_out = NULL;
	int bk_stdout = dup(1);
	int bk_stderr = dup(2);
	try {
		if ( arity_limit < 0 ) {
			std::string out_file = out_file_pre + ".lcmlamp.closed";
			fp_out = fopen(out_file.c_str(), "w");
			if (fp_out == NULL)
				throw std::string("Can't open file : " + out_file);
			dup2(fileno(fp_out), 1);
			dup2(fileno(fp_out), 2);
			lam = LCMWrap_LAMP('C', n1, p_mode, input_file_chr, sig_level, arity_limit);
		} else {
			std::string out_file = out_file_pre + ".lcmlamp.aritylim" + std::to_string(arity_limit);
			fp_out = fopen(out_file.c_str(), "w");
			if (fp_out == NULL)
				throw std::string("Can't open file : " + out_file);
			dup2(fileno(fp_out), 1);
			dup2(fileno(fp_out), 2);
			lam = LCMWrap_LAMP('F', n1, p_mode, input_file_chr, sig_level, arity_limit);
		}
	} catch (std::string &msg) {
		throw msg;
	} catch (...) {
		throw std::string("Error: An unexpected error occurred while trying to run LCMLAMP.");
	}
	fflush(stdout);
	dup2(bk_stdout, 1);
	dup2(bk_stderr, 2);
	if (fp_out != NULL) fclose(fp_out);

	if (lam < 0)
		throw std::string("Error: An unexpected error occurred while trying to run LCMLAMP.");
	
	lam = lam - 1;
	return lam;
}

/**
 * Construct frequent patterns list.	This method use lcm53 program.
 * @param input_file
 * @param low_sup The number of minimum support. get item set that appeare abobe min_sup.
 * @param arity_limit The limit to appriori depth.
 */
void LCM::frequentPatterns(const std::string& input_file, int low_sup, int arity_limit) {
	// If frequent pattern has already serched, then return.
	if (getIndex( low_sup ) <= constructed_index)
		return;
	char input_file_chr[FILE_NAME_BUF], out_file_chr[FILE_NAME_BUF];
	strncpy(input_file_chr, input_file.c_str(), FILE_NAME_BUF-1);
	
	std::string out_dir = input_file + ".results.lcm";
	boost::filesystem::path ndir(out_dir);
	if (! boost::filesystem::exists(ndir)) {
		boost::filesystem::create_directory(ndir);
	}
	std::list<std::string> out_file_s;
	boost::split(out_file_s, input_file, boost::is_any_of("/"));
	std::string out_file_name = out_file_s.back();
	std::string out_file_pre = out_dir + "/" + out_file_name;

	int upper_sup = max_support - constructed_index - 1;
	// Run LCM
	int ret = 0;
	std::string out_file;
	int bk_stdout = dup(1);
	int bk_stderr = dup(2);
	try {
		dup2(fd_outlog, 1);
		dup2(fd_outlog, 2);
		// If the arity size is not limited, run LCM to get closed frequent pattern.
		if ( arity_limit < 0 && constructed_index > -1 ) {
			out_file = out_file_pre + ".lowsup" + std::to_string(low_sup) + ".upsup" + std::to_string(upper_sup) + ".closed";
			strncpy(out_file_chr, out_file.c_str(), FILE_NAME_BUF-1);
			ret = LCMWrap_freq('C', upper_sup, low_sup, arity_limit, input_file_chr, out_file_chr);
		}
		else if ( arity_limit < 0 && constructed_index == -1 ) {
			out_file = out_file_pre + ".lowsup" + std::to_string(low_sup) + ".closed";
			strncpy(out_file_chr, out_file.c_str(), FILE_NAME_BUF-1);
			ret = LCMWrap_freq('c', upper_sup, low_sup, arity_limit, input_file_chr, out_file_chr);
		}
		else if ( arity_limit >= 0 && constructed_index > -1) {
			out_file = out_file_pre + ".lowsup" + std::to_string(low_sup) + ".upsup" + std::to_string(upper_sup) 
					   + ".aritylim" + std::to_string(arity_limit);
			strncpy(out_file_chr, out_file.c_str(), FILE_NAME_BUF-1);
			ret = LCMWrap_freq('F', upper_sup, low_sup, arity_limit, input_file_chr, out_file_chr);
		}
		else {
			out_file = out_file_pre + ".lowsup" + std::to_string(low_sup) +".aritylim" + std::to_string(arity_limit);
			strncpy(out_file_chr, out_file.c_str(), FILE_NAME_BUF-1);
			ret = LCMWrap_freq('f', upper_sup, low_sup, arity_limit, input_file_chr, out_file_chr);
		}
	} catch (std::string &msg) {
		throw msg;
	} catch (...) {
		throw std::string("Error: An unexpected error occurred while trying to run LCM.");
	}
	fflush(stdout);
	dup2(bk_stdout, 1);
	dup2(bk_stderr, 2);
	
	if (ret != 0)
		throw std::string("Error: An unexpected error occurred while trying to run LCM.");
	
	// Read the file of LCM result
	readResultLCMFile( out_file, low_sup, upper_sup );

	// Update the total number of transactions
	int total = 0;
	if ((low_sup < max_support) && (constructed_index > -1)) 
		total = frequent_list[ getIndex( low_sup ) - 1 ]->getTotal();
	for (int i = upper_sup; low_sup - 1 < i; i--) {
		Node* node = frequent_list[ getIndex( i ) ];
		total = total + node->getItemSetSize();
		node->setTotal(total);
	}
	
	constructed_index = getIndex( low_sup );
}

/**
 * Read LCM result file and return itemset list.
 * @param result_lcm_file the result filename of running LCM
 * @param low_sup the minimum support when LCM has run
 * @param upper_sup the maximum support when LCM has run (the frequent_list re-constructed between min_sup to max_sup)
 */
void LCM::readResultLCMFile(const std::string& result_lcm_file, int low_sup, int upper_sup) {
	// Initialize re-constructed nodes
	for (int i = low_sup; i <= upper_sup; i++) {
		Node* node = new Node();
		frequent_list[ getIndex( i ) ] = node;
	}
	
	// convert output of lcm_basic to item set list
	try {
		std::ifstream f(result_lcm_file);
		if (!f.is_open())
			throw std::string("Can't open file : " + result_lcm_file);
		std::string itemset_line;
		std::getline(f, itemset_line);
		std::string transactions_line = "";
		while (f.good()) {
			std::getline(f, transactions_line);
			// if line startswith space, this line is ignored
			if (itemset_line.at(0) == ' ') {
				std::getline(f, itemset_line);
				continue;
			}
			// convert output of lcm to item set list
			boost::trim(itemset_line); 
			std::vector<std::string> s;
			boost::split(s, itemset_line, boost::is_space());
			
			boost::trim(transactions_line); 
			std::vector<std::string> transactions;
			boost::split(transactions, transactions_line, boost::is_space());
			// if itemset is not empty, add itemset to itemset_list
			if (1 < s.size()) {
				std::vector<int> itemset;
				for (int i = 0; i < (int)s.size() - 1; i++) {
					itemset.push_back(std::stoi(s[i]));
				}
				std::vector<int> transactionset;
				if (transactions_line.length() == 0) {
					transactionset.clear();
				}
				else {
					for (int i = 0; i < (int)transactions.size(); i++) {
						transactionset.push_back(std::stoi(transactions[i]));
					}
				}
				int support = transactions.size();
				Node* node = NULL;
				int node_index = getIndex( support );
				if ( node_index < 0 )
					node = frequent_list[ 0 ];
				else
					node = frequent_list[ getIndex( support ) ];
				node->addItemSet( itemset, transactionset );
			}
			std::getline(f, itemset_line);
		}
		f.close();
	} catch (std::string &msg) {
		throw msg;
	} catch (...) {
		throw std::string("Error: An unexpected error occurred while trying to read Result LCM File.");
	}
}
