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
 * File:   Transaction.cpp
 * Author: miura
 * 
 * Define Class that indicate a transaction.
 * The transaction has name and item set
 * 
 * Created on 2016/01/05, 15:59
 */

#include "Transaction.h"

/**
 * Constructor
 * @param name Transaction name (Gene name)
 */
Transaction::Transaction(const std::string& name) {
	this->name = name;
	id = 0;
	value = std::numeric_limits<double>::quiet_NaN();
}

/**
 * Destructor
 */
Transaction::~Transaction() {
}

/**
 * Add item to this instance
 * @param item
 */
void Transaction::addItem(int item) {
	this->itemset.push_back(item);
}

/**
 * Set value
 * @param value
 */
void Transaction::setValue(double value) {
	this->value = value;
}

/**
 * This function is used for sort transaction list in u_test
 * @param t1
 * @param t2
 * @return 
 */
bool Transaction::comparator(const Transaction* t1, const Transaction* t2) {
	return t1->value < t2->value;
}

/**
 * Claas data output
 * @param os Output stream
 * @param t Output data
 * @return Output stream
 */
std::ostream& operator<<(std::ostream& os, const Transaction* t) {
	os << t->id << " " << t->name << " ";
	for (int i : t->itemset) {
		os << i << " ";
	}
	os << t->value << std::endl;
	
	return os;
}
