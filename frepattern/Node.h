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
 * File:   Node.h
 * Author: miura
 *
 * Created on 2016/01/13, 14:59
 */

#ifndef NODE_H
#define NODE_H

#include <vector>
#include "../Transaction.h"

/** @file */
/** This class manage LCM results
 */
class Node {
public:
	/** item set
	 */
	typedef struct {
		std::vector<int>* item_list; /**< item list. */
		std::vector<int>* tran_list; /**< transaction list. */
	} itemset_t;
public:
	Node();
	virtual ~Node();
	
	/**
	 * set value of bound.
	 * @param bound value of bound
	 */
	void setBound(double bound) { this->bound = bound; }
	/**
	 * get value of bound.
	 * @return value of bound
	 */
	double getBound() { return bound; }
	/**
	 * set number of total.
	 * @param total number of total
	 */
	void setTotal(int total) { this->total = total; }
	/**
	 * set number of total.
	 * @return number of total
	 */
	int getTotal() { return total; }
	/**
	 * add item set
	 * @param item_tuple item set
	 * @param tran transaction set
	 */
	void addItemSet(const std::vector<int>& item_tuple, const std::vector<int>& tran) {
		itemset_list.push_back(new itemset_t{
			new std::vector<int>(item_tuple), new std::vector<int>(tran)
		});
	}
	/**
	 * get size of item set
	 * @return size of item set
	 */
	int getItemSetSize() { return itemset_list.size(); }
	/**
	 * get item set.
	 * @return item set
	 */
	const std::vector<itemset_t*>& getItemSet() { return itemset_list; }
	/**
	 * get item set.
	 * @param i Index number
	 * @return item set
	 */
	const std::vector<int>* getItemSet(int i) { return itemset_list[i]->item_list; }
	/**
	 * get transactions set
	 * @param i Index number
	 * @return transactions set
	 */
	const std::vector<int>* getTransactionSet(int i) { return itemset_list[i]->tran_list; }
	
private:
	double bound; /**< the value to use upper/lower bound */
	int total; /**< total number of item list. */
	std::vector<itemset_t*> itemset_list; /**< item set  */

};

#endif /* NODE_H */

