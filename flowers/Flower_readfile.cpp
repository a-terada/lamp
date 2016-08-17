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

#include "Flower_readfile.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include <set>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
using namespace std;
using namespace boost;

/**
 * Constructor
 */
Flower_readfile::Flower_readfile() {
	pval_f = NULL;
}

/**
 * Destructor
 */
Flower_readfile::~Flower_readfile() {
	if (pval_f != NULL) {
		delete pval_f;
		pval_f = NULL;
	}
}

/**
 * Read result file of lamp
 * @param resfname file name of lamp results
 * @param csvfname The file includes associations between TFs and genes.
 * @param tabfname Each line indicates a gene. The column1 is gene name.
 * @return 
 */
double Flower_readfile::readResult(string resfname, string csvfname, string tabfname) {
	bool isResult = false;
	int alternative = 1;
	double significance;
	ifstream ifs(resfname, ios::in);
	if (ifs.fail()) {
		cerr << resfname << " cannot be opened." << endl;
		exit(EXIT_FAILURE);
	}
	string line;
	vector<string> parts;
	string pval_procedure;

	// read the data line-by-line
	while (getline(ifs, line)) {
		boost::trim(line);
		// read the threshold and correction factor        
		if (line.find("Correction factor:") != string::npos) {
			algorithm::split(parts, line, is_any_of(" "));
			vector<string>parts4;
			algorithm::split(parts4, parts[4], is_any_of(","));
			double thresh = lexical_cast<double>(parts4[0]);
			factor = lexical_cast<double>(parts[7]);
			threshFactor = thresh * factor;
			continue;
		}

		// read the csv filename
		if (line.find("# item-file", 0) == 0 && csvfname.length() == 0) {
			algorithm::split(parts, line, is_any_of(" "));
			csvfname = parts[2];
			continue;
		}
		// read the tab filename
		if (line.find("# value-file", 0) == 0 && tabfname.length() == 0) {
			algorithm::split(parts, line, is_any_of(" "));
			tabfname = parts[2];
			continue;
		}
		// read the significance level
		if (line.find("# significance-level:", 0) == 0) {
			algorithm::split(parts, line, is_any_of(" "));
			significance = lexical_cast<double>(parts[2]);
			continue;
		}
		// read the procedure to compute P-value
		if (line.find("# P-value computing procedure", 0) == 0) {
			algorithm::split(parts, line, is_any_of(" "));
			pval_procedure = parts[4];
			string alternative_str;
			if (parts.size() > 5) {
				string tmp_al = parts[5];
				alternative_str = tmp_al.substr(1, tmp_al.length());
				if (alternative_str == "greater") {
					alternative = 1;
				} else if (alternative_str == "less") {
					alternative = -1;
				} else if (alternative_str == "two.sided") {
					alternative = 0;
				}
			}
			continue;
		}

		//stop inputting the values
		if (line.find("Time", 0) == 0) {
			isResult = false;
		}

		//add one set of values
		if (isResult == true) {
			line = algorithm::replace_all_copy(line, " ", "\t");
			algorithm::split(parts, line, is_any_of("\t"));
			vector<string> namearray;
			string tmp = parts[3];
			algorithm::split(namearray, tmp, is_any_of(L","));
			//            string resr = namearray[1];
			int nn = namearray.size();
			if (nn == 1) {
				motifRpvalue.push_back(lexical_cast<double>(parts[1]));
				motifApvalue.push_back(lexical_cast<double>(parts[2]));
				motifName.push_back(namearray[0]);
				motifNgenes.push_back(lexical_cast<int>(parts[5]));
				motifSscore.push_back(lexical_cast<double>(parts[6]));
			} else {
				combiRank.push_back(lexical_cast<int>(parts[0]));
				combiRpvalue.push_back(lexical_cast<double>(parts[1]));
				combiApvalue.push_back(lexical_cast<double>(parts[2]));
				combiName.push_back(namearray);
				combiNgenes.push_back(lexical_cast<int>(parts[5]));
				combiSscore.push_back(lexical_cast<double>(parts[6]));
			}
		}
		// start inputting values at the next line
		if (line.find("Rank", 0) == 0) {
			isResult = true;
		}
	}
	ifs.close();
	// add motifs missed in the input file
	for (int i = 0; i < (int)combiName.size(); i++) {
		vector<string> nameset = combiName[i];
		for (int k = 0; k < (int)nameset.size(); k++) {
			// search for the motif and specify its ID
			int nameid = -1;
			for (int j = 0; j < (int)motifName.size(); j++) {
				if (motifName[j] == nameset[k]) {
					nameid = j;
					break;
				}
			}
			// append the name because it is not described in the input file
			if (nameid < 0) {
				motifRpvalue.push_back(-1.0);
				motifApvalue.push_back(-1.0);
				motifName.push_back(nameset[k]);
				motifNgenes.push_back(0);
				motifSscore.push_back(0);
			}
		}
	}
	// retrieve values of notifs not described in the input file
	ReadFile* readFile;
	readFile = new ReadFile();
	readFile->readFiles(csvfname, tabfname, ',');
	vector<Transaction*> transaction_list = readFile->getTransaction_list();
	if (alternative < 0) {
		LampCore lamp;
		lamp.reverseValue(transaction_list, pval_procedure);
	}
	if (pval_procedure == "fisher") {
		pval_f = new Functions4fisher(transaction_list, alternative);
	} else if (pval_procedure == "chi") {
		pval_f = new Functions4chi(transaction_list, alternative);
	} else if (pval_procedure == "u_test") {
		pval_f = new Functions4u_test(transaction_list, alternative);
	} else {
		throw string("Error: Unknown p-value procedure.");
	}
	map<string, int> colname2id;
	readFile->getColname2id(colname2id);
	for (int i = 0; i < (int)motifName.size(); i++) {
		if (motifApvalue[i] > 0.0) {
			continue;
		}
		// retrieve values of the current motif
		vector <string> mnamelist;
		algorithm::split(mnamelist, motifName[i], is_any_of(","));
		// cout << "... retrieving values of " << algorithm::join(mnamelist, " ");
		//        float p_value, int down_size = pval_func(csvfname, tabfname, mnamelist, ',', alternative);
		int down_size;
		double p_value = calPValue(transaction_list, colname2id, mnamelist, down_size);

		motifRpvalue[i] = p_value;
		motifApvalue[i] = p_value * factor;
		motifNgenes[i] = down_size;
	}
	// return
	return significance;
}

/**
 * Calculate p-value by using functions class
 * @param transaction_list list of transactions
 * @param colname2id dictionary about mapping column name to column id
 * @param itemset_str_lst list of item
 * @param down_size return value
 * @return p-value
 */
double Flower_readfile::calPValue(const vector<Transaction*> &transaction_list,
		map<string, int>& colname2id, const vector<string>& itemset_str_lst,
		int &down_size) {
	vector<int> itemset;
	for (string name : itemset_str_lst) {
		int item_id = colname2id[name];
		itemset.push_back(item_id + 1);
	}
	vector<int>flag_transactions_id;
	vector<int> intersect;
	for (int i = 0; i < (int)transaction_list.size(); i++) {
		Transaction* t = transaction_list[i];
		const std::vector<int>& titemset = t->getItemset();
		intersect.clear();
		set_intersection(itemset.begin(), itemset.end(), titemset.begin(), titemset.end(), back_inserter(intersect));
		if (intersect.size() == itemset.size()) {
			flag_transactions_id.push_back(i);
		}
	}
	//	double p_value, stat_value = func4.calPValue(transaction_list, flag_transactions_id);
	double stat_value;
	double p_value = pval_f->calPValue(flag_transactions_id, stat_value);
	int n = transaction_list.size();
	int n1 = pval_f->getN1();
	cout << "p-value:" << p_value << "(N:" << n << ", n1:" << n1 << ", x:" << flag_transactions_id.size() << ", a:" << stat_value << ")" << endl;
	down_size = flag_transactions_id.size();
	return (p_value);
}
