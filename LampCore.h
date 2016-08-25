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

#ifndef LAMPCORE_H
#define LAMPCORE_H

#include <cstdlib>
#include <ctime>
#include <climits>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include "ReadFile.h"
#include "Transaction.h"
#include "functions/FunctionsSuper.h"
#include "functions/Functions4fisher.h"
#include "functions/Functions4u_test.h"
#include "functions/Functions4chi.h"
#include "frepattern/LCM.h"

#define __LAMP_VER__ "3.0.0" /**< version of LAMP */

/** @file */
/** Execute LAMP function
 */
class LampCore {
public:
	/** This struct is storing the LAMP result
	 */
	typedef struct {
		const std::vector<int>* item_set; /**< item set */
		double p; /**< P-value */
		int len; /**< support value */
		double stat_score; /**< score(statistics) */
	} enrich_t;
	
public:
	LampCore();
	virtual ~LampCore();
	void clean();
	
	void run(std::string& transaction_file, std::string& flag_file, double threshold,
		std::string& set_method, int max_comb,
		std::string& log_file, int alternative);
	void print(std::string& transaction_file, std::string& flag_file, double threshold,
        std::string& set_method, int alternative);
	
	/**
	 * get result set
	 * @return result set
	 */
	const std::vector<enrich_t*>& getEnrich_lst() { return enrich_lst;}
	/**
	 * get value of K
	 * @return K
	 */
	int getK() { return k; }
	/**
	 * 
	 * @return 
	 */
	int getLam_star() { return lam_star; }
	/**
	 * get list about mapping column id to column name
	 * @return list about mapping column id to column name
	 */
	std::vector<std::string*>& getColumnid2name() { return readFile.getColumnid2name(); }
	/**
	 * get value of N1
	 * @return value of N1
	 */
	int getN1() { return func_f->getN1(); }
	/**
	 * get size of transaction set.
	 * @return size of transaction set.
	 */
	int getTransactionSize() { return readFile.getTransaction_list().size(); }
	
	/**
	 * Comparator of enrich_t
	 * @param x left object 
	 * @param y right object
	 * @return 
	 */
	static bool cmpEnrich(const enrich_t* x, const enrich_t* y) {
		return x->p < y->p;
	}
	
	/**
	 * get version number
	 * @return version string
	 */
	static std::string getVersion() { return __LAMP_VER__; }
	
	const std::vector<std::string> BINARY_METHODS{"fisher", "chi"}; /**< name list of Binary methods */
	void reverseValue(std::vector<Transaction*>& transaction_list, const std::string& set_method );
	
protected:
	int runMultTest(const std::vector<Transaction*>& transaction_list, std::string& trans4lcm,
		double threshold, std::string& set_method, int max_comb,
		std::ostream& outlog, int alternative, double max_lambda);
	double maxLambda(std::vector<Transaction*>& transaction_list);
	int depthFirst( const std::string& trans4lcm, int max_comb,
		int n1, double threshold, int p_mode );
	int breadthFirst( const std::string& trans4lcm, int max_comb,
		double threshold, int lam, std::ostream& outlog );
	double calBound(int min_sup );
	void fwerControll(double max_lambda, double threshold, std::ostream& outlog);
	void outputResult(std::string& transaction_file, std::string& flag_file,
		double threshold, std::string& set_method, std::vector<std::string*>& columnid2name,
		std::vector<Transaction*>& transaction_list,
		int alternative );
	
	FunctionsSuper *func_f = NULL; /**< Instance to perform the statistical test. */
	std::vector<enrich_t*> enrich_lst; /**< results */
	ReadFile readFile; /**< Query data */
	LCM *fre_pattern = NULL; /**< Instance to run LCM. */
	int k; /**< result */
	int lam_star; /**< result */
	std::chrono::high_resolution_clock::time_point starttime, /**< start time */
		correction_term_time, /**< termination time of correction */
		finish_test_time; /**< finish time */
};

#endif /* LAMPCORE_H */

