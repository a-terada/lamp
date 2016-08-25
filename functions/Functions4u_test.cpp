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
 * File:   Functions4u_test.cpp
 * Author: miura
 * 
 * This class calculate function f that means minimum p-value (MASL).
 * 
 * Created on 2016/01/06, 17:27
 */

#include <vector>

#include "Functions4u_test.h"

/**
 * Constructor
 * @param transaction_list list of transactions
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
Functions4u_test::Functions4u_test(const std::vector<Transaction*>& transaction_list, int alternative) :
		FunctionsSuper(transaction_list, alternative)
{
	__t_size = transaction_list.size(); // all transaction size
	this->alternative = alternative; // alternative hypothesis. greater -> 1, less -> -1, two.sided -> 0.
	calTime = 0; // Total number of calculate P-value
}

/**
 * Destructor
 */
Functions4u_test::~Functions4u_test() {
}

/**
 * This function calculates the minimum p-value which support size is x.
 * That is, calculates MASL.
 * The z-value that minimum p-value is mean/var
 * @param x
 * @return 
 */
double Functions4u_test::funcF(int x) {
	// calculate p-value if the group_x is consisted from max x transactions.
	
	// mean and variance which is used in U test.
	double size_x = x;
	if (__t_size < size_x)
		size_x = __t_size;		
	double size_y = __t_size - x;
	if (size_y < 0)
		size_y = 0;
	double mean_u = (size_x * size_y) / 2;
	double var_u = size_x * size_y * (size_x + size_y + 1) / 12;
	
	double min_z = mean_u / std::sqrt(var_u); // minimum z-value limited x.
	double p = stdNorDistribution(min_z); // p-value if transaction divided into max x and other.
	return p;
}

/**
 * calculate p-value by t test.
 * @param flag_transactions_id return value
 * @param score return value
 * @return 
 */
double Functions4u_test::calPValue(std::vector<int>& flag_transactions_id, double& score) {
	std::vector<Transaction*> in_t_list, out_t_list;
	__divideGroup(flag_transactions_id, in_t_list, out_t_list); // dvide transactions to itemset or not.
	double z_value;
	double p_value = __uTest(in_t_list, out_t_list, z_value);
	if (alternative == 0) {
		p_value = std::min( p_value * 2., 1.0 );
	}
	else {
		if (z_value < 0)
			p_value = 1. - p_value;
		if (alternative < 0)
			z_value = -z_value;
		calTime = calTime + 1;
	}
	score = z_value;
	return p_value;
}

/**
 * divide transaction_list to two groups.
 * One group is transactions which included itemset, the other is not.
 * @param frequent_itemset itemset of transaction
 * @param in_t_list return value
 * @param out_t_list return value
 */
void Functions4u_test::__divideGroup(std::vector<int>& frequent_itemset,
		std::vector<Transaction*>& in_t_list, std::vector<Transaction*>& out_t_list) {
	// If itemset of t contains test itemset, t puts in_t_list.
	// Else, t puts out_t_list
	for (int i = 0; i < (int)transaction_list.size(); i++) {
		Transaction* t = transaction_list[i];
		if (std::find( frequent_itemset.begin(), frequent_itemset.end() , i ) != frequent_itemset.end())
			in_t_list.push_back(t);
		else
			out_t_list.push_back(t);
	}
}

/**
 * Calculate p-value by using Mann-Whitney U test
 * @param tgroup_x test group 1. This is consisted of transactions.
 * @param tgroup_y test group 2. This is consisted of transactions.
 * @param z_value return value
 * @return 
 */
double Functions4u_test::__uTest(std::vector<Transaction*>& tgroup_x, std::vector<Transaction*>& tgroup_y, double& z_value) {
	double u_value = __uValue(tgroup_x, tgroup_y); // u-value of two groups.
	// z value of u-value
	//mean_u, var_u = __calStatValue(tgroup_x, tgroup_y);
	double size_x = tgroup_x.size();
	double size_y = tgroup_y.size();
	double mean_u = (size_x * size_y) / 2.0;
	double var_u = size_x * size_y * (size_x + size_y + 1.0) / 12.0;
	if (var_u == 0) {
		z_value = 0;
		return 1.0;
	}
	z_value = (u_value - mean_u) / std::sqrt(var_u);
	
	// calculate p-value from z_value
	// this value approximation of standard normal distribution
	return stdNorDistribution(z_value);
}

/**
 * Calculate u value which measurs difference rank sum of two groups.
 * @param tgroup_x test group 1. This is consisted of transactions.
 * @param tgroup_y test group 2. This is consisted of transactions.
 * tgroup_x and t_group_y already sorted by transaction value.
 * @return 
 */
double Functions4u_test::__uValue(std::vector<Transaction*>& tgroup_x, std::vector<Transaction*>& tgroup_y) {
	double u_value = 0.0;
	int previous_u_x_min = 0; // The rank of transaction in previous search
	int previous_u_x_max = 0; // The rank of transaction in previous search
	double previous_value = std::numeric_limits<double>::quiet_NaN(); // The previous expression value
	int left_index = 0; // The start point of searching value.
	int right_index = tgroup_y.size() - 1; // The end point of searching value.
	for (Transaction* t_x : tgroup_x) {
		// u_x_min: rank sum of transaction which the value < t_x in tgroup_y
		// u_x_max: rank sum of transaction which the value <= t_x in tgroup_y
		int u_x_min = std::numeric_limits<int>::quiet_NaN();
		int u_x_max = std::numeric_limits<int>::quiet_NaN();
		// If t_x.value is equal to previous one, u_x_min and u_x_max are also equals.
		if (t_x->getValue() == previous_value) {
			u_x_min = previous_u_x_min;
			u_x_max = previous_u_x_max;
		}
		// Caluclate u_value because tgroup_x value exists between tgroup_y range
		else {
			__binarySearch(t_x->getValue(), tgroup_y, left_index, right_index, u_x_min, u_x_max);
			left_index = u_x_max;
		}
		// Add rank of t_x to u_value
		u_value = u_value + ((double)u_x_min + (double)u_x_max) / 2.0;
		previous_u_x_min = u_x_min;
		previous_u_x_max = u_x_max;
		previous_value = t_x->getValue();
	}
	return u_value;
}

/**
 * Search the group which t.value less than threshold.
 * @param threshold search value
 * @param tgroup list to search threshold
 * @param left_index start index of tgroup to search
 * @param right_index end index of tgroup to search
 * @param u_x_min return value
 * @param u_x_max return value
 */
void Functions4u_test::__binarySearch(double threshold, std::vector<Transaction*>& tgroup,
		int left_index, int right_index, int& u_x_min, int& u_x_max) {
	if ((int)tgroup.size() <= left_index) {
		u_x_min = tgroup.size();
		u_x_max = tgroup.size();
		return;
	}
		
	// compare threshold to min and max value
	if (threshold < tgroup[left_index]->getValue()) {
		u_x_min = left_index;
		u_x_max = left_index;
		return;
	}
	if (tgroup[right_index]->getValue() < threshold) {
		u_x_min = right_index;
		u_x_max = right_index;
		return;
	}
		
	// serach the index which all larger indexes are more than threshold.
	int mid_index = -1;
	while (left_index <= right_index) {
		mid_index = (left_index + right_index) / 2;
		Transaction* mid_transaction = tgroup[mid_index];

		// When the mid.value = threshold, finish the search.
		if (mid_transaction->getValue() == threshold)
			break;
		// When the check value less than threshod, go to serach right.
		else if (mid_transaction->getValue() < threshold)
			left_index = mid_index + 1;
		// When the check value >= threshold, go to search left.
		else
			right_index = mid_index - 1;
	}
	
	// search the same range of the threshold
	Transaction* mid_transaction = tgroup[mid_index];
	if (mid_transaction->getValue() == threshold) {
		int min_index = mid_index;
		int max_index = mid_index;
		Transaction* min_transaction = tgroup[min_index];
		Transaction* max_transaction = tgroup[max_index];
		while (threshold <= min_transaction->getValue()) {
			min_index = min_index - 1;
			if (min_index < 0)
				break;
			min_transaction = tgroup[min_index];
		}
		while (max_transaction->getValue() <= threshold) {				
			max_index = max_index + 1;
			if ((int)tgroup.size() <= max_index)
				break;
			max_transaction = tgroup[max_index];
		}
		u_x_min = min_index + 1;
		u_x_max = max_index;
	}
	// not found the threshold in tgroup.
	// in this case, min_index > max_index
	else if (mid_transaction->getValue() < threshold) {
		while (mid_transaction->getValue() < threshold) {
			mid_index = mid_index + 1;
			mid_transaction = tgroup[mid_index];
		}
		u_x_min = mid_index;
		u_x_max = mid_index;
	}
	// not found the threshod in tgroup and the case of mid_transaction value > threshold
	// In this case, min_index > max_index
	else {
		while (threshold < mid_transaction->getValue()) {
			mid_index = mid_index - 1;
			mid_transaction = tgroup[mid_index];
		}
		u_x_min = mid_index + 1;
		u_x_max = mid_index + 1;
	}
	return;
}
