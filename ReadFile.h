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
 * File:   ReadFile.h
 * Author: miura
 *
 * Define methods to read transaction and flag files
 * 
 * Created on 2016/01/05, 13:22
 */

#ifndef READFILE_H
#define READFILE_H

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <limits>
#include <boost/algorithm/string.hpp>
#include "Transaction.h"

/** @file */
/** Read transaction and flag file and store transaction matrix 
 */
class ReadFile {
public:
	ReadFile ();
	virtual ~ReadFile ();

	void readFiles (const std::string& transaction_file, const std::string& value_file, const char delm);
	/**
	 * get list of transactions
	 * @return transactions
	 */
	std::vector<Transaction*>& getTransaction_list() { return transaction_list; }
	/**
	 * get list about mapping column id to column name
	 * @return list about mapping column id to column name
	 */
	std::vector<std::string*>& getColumnid2name() { return columnid2name; }
	/**
	 * get dictionary about mapping column name to column id
	 * @param colname2id dictionary about mapping column name to column id
	 */
	void getColname2id(std::map<std::string, int>& colname2id) {
		int index = 0;
		for (std::string* namePtr : columnid2name) {
			colname2id.insert(std::make_pair(*namePtr, (int)index));
			index++;
		}
	}

private:
	void readTransactionFile (const std::string& transaction_file, const char delm);
	void readValueFile (const std::string& value_file, const char delm);
	void checkTransName(const std::string& transaction_file);

protected:
	std::vector<Transaction*> transaction_list; /**< list of transactions */
	std::vector<std::string*> columnid2name; /**< list about mapping column id to column name */
	std::map<std::string, int> gene2id; /**< dictionary that gene name -> transaction ID */
	
};

#endif /* READFILE_H */

