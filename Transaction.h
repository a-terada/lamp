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
 * File:   Transaction.h
 * Author: miura
 *
 * Created on 2016/01/05, 15:59
 */

#ifndef TRANSACTION_H
#define TRANSACTION_H

#include <iostream>
#include <limits>
#include <vector>
#include <string>

/** @file */
/** this class manage transaction data.
 */
class Transaction {
public:
	Transaction(const std::string& name);
	virtual ~Transaction();
	
	/**
	 * get Identifier
	 * @return 
	 */
	int getID() { return id; }
	/**
	 * set Identifier
	 * @param id
	 */
	void setID(int id) { this->id = id; }
	/**
	 * get Transaction name (Gene name)
	 * @return Transaction name (Gene name) 
	 */
	const std::string& getName() { return name; }
	void addItem(int item);
	/**
	 * get Items that is belonged to transaction (associated TFs set)
	 * @return Items
	 */
	const std::vector<int>& getItemset() { return itemset; }
	/**
	 * get value
	 * @return value
	 */
	double getValue() { return value; }
	void setValue(double value);
	
	static bool comparator(const Transaction* t1, const Transaction* t2);
	
	friend std::ostream& operator<<(std::ostream& os, const Transaction* t);

private:
	int id; /**< Identifier */
	std::string name; /**< Transaction name (Gene name) */
	std::vector<int> itemset; /**< Items that is belonged to transaction (associated TFs set) */
	double value; /**< Indicate this transaction related to feature.
				 * If fisher's exact test, the value is 1 or 0.
				 * If Mann-Whitney's u-test, the value takes any value.
				 */

};

#endif /* TRANSACTION_H */

