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
 * File:   Functions4chi.cpp
 * Author: miura
 * 
 * This source includes calculate P-value and MASL of the chi-square test.
 * 
 * Created on 2016/01/06, 17:43
 */

#include <vector>

#include "Functions4chi.h"

/**
 * Constructor
 * @param transaction_list list of transactions
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
Functions4chi::Functions4chi(const std::vector<Transaction*>& transaction_list, int alternative) : 
		FunctionsSuper(transaction_list, alternative)
{
	__t_size = transaction_list.size(); // all transaction size
	__f_size = sumValue(); // transaction size which have flag = 1 (n1)
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
Functions4chi::~Functions4chi() {
}

/**
 * This function calculates the minimum p-value which support size is x.
 * @param x
 * @return 
 */
double Functions4chi::funcF(int x) {
	double p1 = 1.0f, p2 = 1.0f;
	double chi1 = 0.0f, chi2 = 0.0f;
	double total_row1 = __f_size;
	double total = __t_size;
	// when x < n_u
	if (x < total_row1) {
		double ovalues[2][2] = {{(double)x, 0}, {total_row1 - (double)x, total - total_row1}};
		chi1 = __probabilityTable( ovalues );
		p1 = __chi2pval( chi1 );
		ovalues[0][0] = 0;
		ovalues[0][1] = (double)x;
		ovalues[1][0] = total_row1;
		ovalues[1][1] = total - total_row1 - (double)x;
		chi2 = __probabilityTable( ovalues );
		p2 = __chi2pval( chi2 );
	}
	// when x >= n_u
	else {
		double ovalues[2][2] = {{total_row1, (double)x - total_row1}, {0, total - (double)x}};
		chi1 = __probabilityTable( ovalues );
		p1 = __chi2pval( chi1 );
		ovalues[0][0] = 0;
		ovalues[0][1] = (double)x;
		ovalues[1][0] = total_row1;
		ovalues[1][1] = total - total_row1 - (double)x;
		chi2 = __probabilityTable( ovalues );
		p2 = __chi2pval( chi2 );
	}
	if (alternative == 0) {
		p1 = std::min( p1 * 2.0, 1.0 );
		p2 = std::min( p2 * 2.0, 1.0 );
	}
	if (p1 < p2){
		return p1;
	} else {
		return p2;
	}
}

/**
 * Calculate probability of occurrence probability about table.
 * @param ovalues
 * @return 
 */
double Functions4chi::__probabilityTable(const double (&ovalues)[2][2]) {
	double means[2][2] = {{0, 0}, {0, 0}};
	__calMeans(ovalues, means); // calculate the exception value
	
	// Yate continuity correction
	double yate_corr = 0;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			if (means[i][j] < 5) {
				yate_corr = 0.5;
				break;
			}
		}
	}
	
	double chi = 0;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			double row = ovalues[i][j];
			double mean = means[i][j];
			double po = (std::abs(row - mean) - yate_corr);
			chi += (po * po) / mean;
		}
	}
	return chi;
}

/**
 * Calculate means
 * @param ovalues
 * @param means
 */
void Functions4chi::__calMeans(const double (&ovalues)[2][2], double (&means)[2][2]) {
	double total = __t_size;
	double total_col1 = __f_size; // the number of all flag 1 transaction (n1)
	double total_col2 = total - total_col1; // the number of all flag 0 transactio (n0)
	double total_row1 = ovalues[0][0] + ovalues[0][1];
	double total_row2 = ovalues[1][0] + ovalues[1][1];
	means[0][0] = (total_row1 * total_col1) / total;
	means[0][1] = (total_row1 * total_col2) / total;
	means[1][0] = (total_row2 * total_col1) / total;
	means[1][1] = (total_row2 * total_col2) / total;
}

/**
 * 
 * @param chi
 * @return 
 */
double Functions4chi::__chi2pval(double chi) {
	if (chi == 0.0) {
		return 1.0;
	}
	else { // dimension = 1
		return stdNorDistribution(std::pow(chi, 0.5));
	}
}

/**
 * Calculate p-value by using chi-square test.
 * @param flag_transactions_id Transactions which have items
 * @param score return value
 * @return 
 */
double Functions4chi::calPValue(std::vector<int>& flag_transactions_id, double& score) {
	double ovalues[2][2] = {{0, 0},{0, 0}};
	contingencyTable( flag_transactions_id, __t_size, __f_size, ovalues );
	double total_row1 = ovalues[0][0] + ovalues[0][1];//sum( ovalues[0] );
	double p = __pvalTable.getValue( total_row1, ovalues[0][0] );
	double chi = __chiTable.getValue( total_row1, ovalues[0][0] );
	if (p < 0) { // calculate P-value and save to the table
		chi = __probabilityTable( ovalues );
		p = __chi2pval( chi );
		if (0 < alternative) {
			if (ovalues[0][0] < (std::min((double)__f_size, total_row1) / 2.0))
				p = 1. - p;
		}
		// when the alternative hypothesis is "two.sided", 
		// the P-value is doubled. 
		else {
				p = std::min( p * 2., 1.0 );
		}
		__pvalTable.putValue( total_row1, ovalues[0][0], p );
		__chiTable.putValue( total_row1, ovalues[0][0], chi );
	}
	score = ovalues[0][0];
	return p;
}
