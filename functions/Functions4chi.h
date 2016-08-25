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
 * File:   Functions4chi.h
 * Author: miura
 *
 * Created on 2016/01/06, 17:43
 */

#ifndef FUNCTIONS4CHI_H
#define FUNCTIONS4CHI_H

#include <cmath>
#include "FunctionsSuper.h"
#include "PvalTable.h"

/** @file */

/** This class calculate P-value and MASL of the chi-square test.
 */
class Functions4chi : public FunctionsSuper {
public:
	Functions4chi(const std::vector<Transaction*>& transaction_list, int alternative);
	virtual ~Functions4chi();

	double funcF(int x) override;
	double calPValue(std::vector<int>& flag_transactions_id, double& score) override;

private:
	double __probabilityTable(const double (&ovalues)[2][2]);
	void __calMeans(const double (&ovalues)[2][2], double (&means)[2][2]);
	double __chi2pval(double chi);

	PvalTable __pvalTable; /**< p-value by using chi-square test */
	PvalTable __chiTable; /**< chi-square-value by using chi-square test */
};

#endif /* FUNCTIONS4CHI_H */

