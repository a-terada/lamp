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
 * File:   LCM.h
 * Author: miura
 *
 * Created on 2016/01/13, 14:46
 */

#ifndef LCM_H
#define LCM_H

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include "Node.h"
#include "../Transaction.h"
#include "../lcm53/LCMWrap.h"

/** @file */
/** Get frequent pattern in transaction lists.
 */
class LCM {
	static const int FILE_NAME_BUF; /**< size of file name buffer */
public:
	LCM(int max_support, int fd_outlog);
	virtual ~LCM();
	
	void makeFile4Lem( const std::vector<Transaction*>& transaction_list, std::string& output_file);
	int runLCMLAMP( const std::string& input_file, int arity_limit, int n1, double sig_level, int p_mode );
	void frequentPatterns(const std::string& input_file, int low_sup, int arity_limit);
	/**
	 * get number of total.
	 * @param min_sup Minimum support
	 * @return number of total
	 */
	int getTotal(int min_sup) { return frequent_list[ getIndex(min_sup) ]->getTotal(); }
	/**
	 * get value of bound.
	 * @param min_sup Minimum support
	 * @return value of bound.
	 */
	double getBound(int min_sup) { return frequent_list[ getIndex(min_sup) ]->getBound(); }
	/**
	 * set value of bound.
	 * @param min_sup Minimum support
	 * @param bound value of bound.
	 */
	void setBound(int min_sup, double bound) {
		Node* node = frequent_list[ getIndex(min_sup) ];
		node->setBound( bound );
	}
	/**
	 * get item set
	 * @param support value of support
	 * @return item set
	 */
	const std::vector<Node::itemset_t*>& getItemsetList(int support) {
		return frequent_list[ getIndex(support) ]->getItemSet();
	}
	/**
	 * get maximum value of the minimum support.
	 * @return maximum value of the minimum support.
	 */
	const int getMax_support(){return max_support;}

private:
	std::vector<Node*> frequent_list; /**< frequent patterns list */
	int max_support; /**< maximum value of the minimum support. */
	int constructed_index; /**< frequent_list is constructed that its index is less than the value. */
	int fd_outlog; /**< file descriptor of log file. */

	int getIndex(int support) { return max_support - support; }
	void readResultLCMFile(const std::string& result_lcm_file, int low_sup, int upper_sup);

};

#endif /* LCM_H */

