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
 * File:   FunctionsSuper.cpp
 * Author: miura
 * 
 * Created on 2016/01/06, 15:28
 */

#include "FunctionsSuper.h"

/**
 * Constructor
 * @param transaction_list list of transactions
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
FunctionsSuper::FunctionsSuper(const std::vector<Transaction*>& transaction_list, int alternative) :
		alternative(alternative), transaction_list(transaction_list)
{
}

/**
 * Destructor
 */
FunctionsSuper::~FunctionsSuper() {
}

/**
 * Total value of transaction
 * @return 
 */
double FunctionsSuper::sumValue() {
	double f_size = 0.0f;
	for (Transaction* t : transaction_list) {
		f_size += t->getValue();
	}
	return f_size;
}

/**
 * Calculate probability of standard normal distribution.
 * this function returns the probability of one-sided test.
 * @param x
 * @return 
 */
double FunctionsSuper::stdNorDistribution(double x) {
	double is_value = -1;
	double y = std::abs(x);
	double c = y * y;
	double p = 0.0f;
	double z = std::exp(-c * 0.5) * pi2;
	if (y < 2.5) {
		for (double i = 20.0; 0 < i; i -= 1.0f) {
			p = i * c / (i * 2.0 + 1.0 + is_value * p);
			is_value = -is_value;
		}
		p = 0.5 - z * y / (1.0 - p);
	} else {
		for (double i = 20.0; 0 < i; i -= 1.0f) {
			p = i / (y + p);
		}
		p = z / (y + p);
	}
	return p;
}

/**
 * Make the contingency table.
 * This function is used by fisher, chi-square test and exact logistic regression.
 * @param flag_transactions_id
 * @param total
 * @param total_col1
 * @param ovalues
 */
void FunctionsSuper::contingencyTable(const std::vector<int>& flag_transactions_id, int total, int total_col1, double (&ovalues)[2][2]) {
	double total_col2 = total - total_col1; // the number of all flag 0 transactio (n0)
	// count trahsaction which contains itemset and flag is 1. (This is indicate a of paper.)
	double total_row1 = flag_transactions_id.size(); // count all size that flag = 1 (x of paper)
	double total_row = 0.0f;
	for (int i : flag_transactions_id) {
		Transaction* t = transaction_list[i];
		// If t flag = 1, then sum_has_flag ++.
		total_row += t->getValue();
	}
	ovalues[0][0] = total_row;
	ovalues[0][1] = total_row1 - ovalues[0][0]; // the number of transaction which contains itemset and flag is 0 (This is indicate b of paper)
	ovalues[1][0] = total_col1 - ovalues[0][0];
	ovalues[1][1] = total_col2 - ovalues[0][1];
}
