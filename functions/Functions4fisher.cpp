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
 * File:   Functions4fisher.cpp
 * Author: miura
 * 
 * This class calculate function f that define in paper. (C(n1, x)/C(n0+n1, x))
 * 
 * Created on 2016/01/06, 16:23
 */

#include "Functions4fisher.h"

/**
 * Constructor
 * @param transaction_list list of transactions
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
Functions4fisher::Functions4fisher(const std::vector<Transaction*>& transaction_list, int alternative) : 
		FunctionsSuper(transaction_list, alternative)
{
	__t_size = transaction_list.size(); // all transaction size
	__f_size = sumValue(); // transaction size which have flag = 1 (n1)
	calTime = 0; // Total number of calculate P-value
	this->alternative = alternative; // alternative hypothesis. greater or less -> 1, two.sided -> 0.
	if (__f_size == 0) {
		throw std::string("Error: There is no up-regulate gene.");
	}
	// Check the transaction value.
	// If the value is not 1 or 0, raise error.
	// Because fisher's exact test does not handle numerical value.
	for (Transaction* t : transaction_list) {
		if (!(t->getValue() == 1.0 || t->getValue() == 0.0)) {
			throw std::string("Error: \"" + t->getName() + "\" value is " + std::to_string(t->getValue()) + ".\n" +
				"       But value is 1 or 0 if you test by fisher's exact test.");
		}
	}
}

/**
 * Destructor
 */
Functions4fisher::~Functions4fisher() {
}

/**
 * This function calculates the minimum p-value which support size is x.
 * @param x support size
 * @return 
 */
double Functions4fisher::funcF(int x) {
	int n1 = __f_size;
	int n1_n0 = __t_size; // the number of all genes
	// If x > n1, then print error and exit calculation.
	// (MASL may not follow the monotonic decrease.)
	int all_x = n1_n0 - n1;
	
	double ans = -1;
	// x < n1 <= n0
	if (x <= n1) {
		ans = __probability(x, x);
	}
	// n1 <= x <= n0
	else if (x <= all_x) {
		ans = __probability(x, n1);
	}
	// x > n0:
	else {
		throw std::string("Error: x > n1, n0. This code cannot consider this case.");
	}
	return ans;
}

/**
 * Calculate probability of occurrence probability about table.
 * @param x total of the first row (the number of targetting gene)
 * @param a the top-left of the table
 * @return 
 */
double Functions4fisher::__probability(int x, int a) {
	double p = __occrTable.getValue(x, a);
	if (p < 0) {
		int n = __t_size;
		int n1 = __f_size;
		int n0 = n - n1;
		int b = x - a;
		int i = 0;
		p = 1.0;
		while (i < a) {
			p = p * (n1 - i) / (a - i); // c(n1, a)
			p = p * (x  - i) / (n - i); // c(n1+n0, x)
			i = i + 1;
		}
		i = 0;
		while (i < b) {
			p = p * (n0 - i) / (b - i); // c(n0, b)
			double minus_denominator = a + i;
			p = p * (x - minus_denominator)/(n - minus_denominator); // c(n1+n0, x)
			i = i + 1;
		}
		__occrTable.putValue(x, a, p);
	}
	return p;
}

/**
 * Calculate p-value by using fisher's exact test.
 * @param flag_transactions_id Transactions which have items
 * @param score return value
 * @return 
 */
double Functions4fisher::calPValue(std::vector<int>& flag_transactions_id, double& score) {
	double ovalues[2][2] = {{0, 0},{0, 0}};
	contingencyTable( flag_transactions_id, __t_size, __f_size, ovalues);
	int total_col1 = __f_size;
	int total_row1 = ovalues[0][0] + ovalues[0][1];//sum( ovalues[0] );
	double p = __pvalTable.getValue( total_row1, ovalues[0][0] );
	if (p < 0) { // calculate P-value and save to the table
		double p0 = __probability(total_row1, ovalues[0][0]);
		p = p0;
		int pos_max = std::min( total_row1, total_col1 );
		// when the alternative hypothesis is "two.sided",
		// the lower case probability is cumulated. 
		if (alternative < 1) {
			int a = 0;
			// cumulate the lower case probability.
			while (a < ovalues[0][0]) {
				double pa = __probability( total_row1, a );
				if (1.E-16 < pa - p0) // pa > p0
						break;
				p = p + pa;
				a = a + 1;
			}
			// cumulate the upper case probability.
			a = pos_max;
			while ( ovalues[0][0] < a ) {
				double pa = __probability( total_row1, a );
				if (1.E-16 < pa - p0) // pa > p0
						break;
				p = p + pa;
				a = a - 1;
			}
		}
		// when the alternative hypothesis is "greater" or "less",
		// the higher/less case probability is cumulated.  
		else {
			int a = ovalues[0][0] + 1;
			while ( a <= pos_max ){
				double pa = __probability( total_row1, a );
				p = p + pa;
				a = a + 1;
			}
		}
		__pvalTable.putValue( total_row1, ovalues[0][0], p );
		calTime = calTime + 1;
	}
	score = ovalues[0][0];
	return p;
}
