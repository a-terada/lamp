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

#include "ReadFile.h"

#include <boost/algorithm/string.hpp>

/**
 * Constructor
 */
ReadFile::ReadFile() {
}

/**
 * Destructor
 */
ReadFile::~ReadFile() {
	for (Transaction* t : transaction_list) {
		delete t;
	}
	transaction_list.clear();
	for (std::string* s : columnid2name) {
		delete s;
	}
	columnid2name.clear();
}

/**
 * Read transaction and flag file and store transaction matrix.
 * @param transaction_file The file includes associations between TFs and genes.
 *                          Each line indicates a gene.
 *                          If gene is targeted by the TF, then value is 1, otherwise 0.
 * @param value_file Each line indicates a gene. The column1 is gene name.
 *                    If gene has the feature, the column2 is 1. The other case 0.
 * @param delm Delimitor of column
 */
void ReadFile::readFiles(const std::string& transaction_file, const std::string& value_file, const char delm) {
	readTransactionFile(transaction_file, delm);
	readValueFile(value_file, delm);
	// sort transaction_list according to transaction_value
	std::sort(transaction_list.begin(), transaction_list.end(), Transaction::comparator);
	// Generate IDs to transactions
	for (int i = 0; i < (int)transaction_list.size(); i++) {
		Transaction* t = transaction_list[i];
		t->setID(i);
	}
	// check transaction names two
	checkTransName(transaction_file); // check transaction names two
}

/**
 * Read transaction file and return list of transactions.
 * @param transaction_file The file includes associations between TFs and genes.
 * @param delm Delimitor of column
 */
void ReadFile::readTransactionFile(const std::string& transaction_file, const char delm) {
	std::unordered_set<std::string> gene_set;
	int line_num = 0;
	int col_size = -1;
	try {
		std::ifstream f( transaction_file );
		std::string str;
		while ( std::getline( f, str ) ) {
			line_num = line_num + 1;
			boost::trim(str);
			std::istringstream strbuf(str);
			std::string token;
			std::vector<std::string> row_list;
			while (std::getline(strbuf, token, delm)) {
				row_list.push_back(token);
			}
			// If line is header line, read column name
			if (line_num == 1) {
				col_size = row_list.size();
				for (int i = 1; i < col_size; i++) {
					std::string colname = row_list[i];
					std::string* namePtr = new std::string(colname);
					columnid2name.push_back(namePtr);
				}
				continue;
			}
			
			std::string t_name = row_list[0];
			if (gene_set.find(t_name) != gene_set.end()) {
				throw std::string("Error: " + t_name + " is contained two or more times in " + transaction_file + ".");
			}
			gene_set.insert(t_name);
					
			// check the number of columns
			if ((int)row_list.size() != col_size) {
				throw std::string("Error in " + transaction_file + "\n" +
					"    The header line contains " + std::to_string(col_size) +
					" columns, while line " + std::to_string(line_num) +
					" contains " + std::to_string(row_list.size()) + " columns.");
			}
				
			Transaction* t = new Transaction(t_name);
			gene2id.insert(std::make_pair(t_name, (int)transaction_list.size()));
			for (int i = 1; i < (int)row_list.size(); i++) {
				int flag = std::stoi(row_list[i]);
				if (flag == 1) {
					t->addItem(i);
				}
				else if (flag == 0) {
					continue;
				}
				else {
					throw std::string("Error: line " + std::to_string(line_num) + " in \'" + transaction_file + "\' contains the value neither 0 or 1.");
				}
			}
			transaction_list.push_back(t);
		}
		f.close();
	} catch (std::string &msg) {
		throw msg;
	} catch (...) {
		throw std::string("Error: An unexpected error occurred while trying to read Transaction File.");
	}
}

/**
 * Read flag file and add information about flags to transaction list.
 * @param value_file Each line indicates a gene. The column1 is gene name.
 * @param delm  Delimitor of column
 */
void ReadFile::readValueFile(const std::string& value_file, const char delm) {
	int line_num = 0;
	std::unordered_set<std::string> gene_set;
	try {
		std::ifstream f( value_file );
		if (!f.is_open())
			throw std::string("Error: " + value_file + " cannot be found.");
		std::string str;
		while ( std::getline( f, str ) ) {
			if (str.find("#") == 0)
				continue;
			line_num = line_num + 1;
			boost::trim(str);
			std::istringstream strbuf(str);
			std::string token;
			std::vector<std::string> row_list;
			while (std::getline(strbuf, token, delm)) {
				row_list.push_back(token);
			}

			// This error raises if value file contains more than two columns.
			if (row_list.size() != 2) {
				throw std::string("Error: line " + std::to_string(line_num) + " in " + value_file + ".\n" +
					"       value-file should contain two columns.");
			}

			boost::trim(row_list[0]);
			std::string genename = row_list[0];
			boost::trim(row_list[1]);
			std::string exp_value_str = row_list[1];

			// This error raises if value cannot be converted to float.
			double exp_value;
			try {
				exp_value = std::stof(exp_value_str);
			} catch (...) {
				throw std::string("Error: line " + std::to_string(line_num) + " in " + value_file + ".\n" +
					"       \'" + exp_value_str + "\' could not be converted string to float.");
			}
			// This error raises if the identical keys are contained more than two times.
			if (gene_set.find(genename) != gene_set.end()) {
				throw std::string("Error: " + genename + " is contained two or more times in " + value_file + ".");
			}
			gene_set.insert(genename);
			// This error raise if gene does not include in itemset file.
			if (gene2id.find(genename) == gene2id.end()) {
				throw std::string("Error: line " + std::to_string(line_num) + " in " + value_file + ".\n" +
					"      \'" + genename + "\' is not contained in itemset file.");
			}
				
			int geneid = gene2id[genename];
			Transaction* t = transaction_list[geneid];
			t->setValue(exp_value);
		}
		f.close();
	} catch (std::string &msg) {
		throw msg;
	} catch (...) {
		throw std::string("Error: An unexpected error occurred while trying to read Value File.");
	}
}

/**
 * Check two files transaction name are same
 * @param transaction_file
 */
void ReadFile::checkTransName(const std::string& transaction_file) {
	for (Transaction* t : transaction_list) {
		if (std::isnan(t->getValue())) {
			throw std::string("\"" + t->getName() + "\" only appears in " + transaction_file);
		}
	}
}
