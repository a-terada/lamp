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
 * File:   Functions4u_test.h
 * Author: miura
 *
 * Created on 2016/01/06, 17:27
 */

#ifndef FUNCTIONS4U_TEST_H
#define FUNCTIONS4U_TEST_H

#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include "FunctionsSuper.h"
#include "../Transaction.h"

/** @file */
/** This class calculate function f that means minimum p-value (MASL).
 */
class Functions4u_test : public FunctionsSuper {
public:
	Functions4u_test(const std::vector<Transaction*>& transaction_list, int alternative);
	virtual ~Functions4u_test();
	
	double funcF(int x) override;
	double calPValue(std::vector<int>& flag_transactions_id, double& score) override;

private:
	void __divideGroup(std::vector<int>& frequent_itemset,
			std::vector<Transaction*>& in_t_list, std::vector<Transaction*>& out_t_list);
	double __uTest(std::vector<Transaction*>& tgroup_x, std::vector<Transaction*>& tgroup_y, double& z_value);
	double __uValue(std::vector<Transaction*>& tgroup_x, std::vector<Transaction*>& tgroup_y);
	void __binarySearch(double threshold, std::vector<Transaction*>& tgroup,
		int left_index, int right_index, int& u_x_min, int& u_x_max);
	
	int calTime; /**< Total number of calculate P-value */
};

#endif /* FUNCTIONS4U_TEST_H */

