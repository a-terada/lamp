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
 * File:   FunctionsSuper.h
 * Author: miura
 *
 * Created on 2016/01/06, 15:28
 */

#ifndef FUNCTIONSSUPER_H
#define FUNCTIONSSUPER_H

#include <cmath>
#include <vector>
#include "../Transaction.h"

/** @file */
/** this class is base to calculate P-value
 */
class FunctionsSuper {
public:
	FunctionsSuper(const std::vector<Transaction*>& transaction_list, int alternative);
	virtual ~FunctionsSuper();
	
	/**
	 * get N1
	 * @return N1
	 */
	int getN1() { return __f_size; }
	/**
	 * get size of all.
	 * @return size of all
	 */
	int getAllSize() { return __t_size; }
	
	double sumValue();
	double stdNorDistribution(double x);

	/**
	 * This function calculates the minimum p-value which support size is x.
	 * @param x support size
	 * @return  
	 */
	virtual double funcF(int x) = 0;
	/**
	 * calculate p-value
	 * @param flag_transactions_id
	 * @param score
	 * @return 
	 */
	virtual double calPValue(std::vector<int>& flag_transactions_id, double& score) = 0;
	
protected:
	void contingencyTable(const std::vector<int>& flag_transactions_id, int total, int total_col1, double (&ovalues)[2][2]);

	int __f_size; /**< transaction size which have flag = 1 (n1) */
	int __t_size; /**< all transaction size */
	int alternative; /**< alternative hypothesis. greater or less -> 1, two.sided -> 0. */
	const std::vector<Transaction*>& transaction_list; /**< transaction list */
	
private:
	const double pi2 = 0.398942280401432677940; /**< constant value */

};

#endif /* FUNCTIONSSUPER_H */

