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
 * File:   flower_readfile.h
 * Author: ooi
 *
 * Created on 2016/01/12, 10:03
 */

#ifndef FLOWER_READFILE_H
#define FLOWER_READFILE_H

#include <map>
#include <string>
#include <vector>
#include "../functions/FunctionsSuper.h"
#include "../functions/Functions4fisher.h"
#include "../functions/Functions4chi.h"
#include "../functions/Functions4u_test.h"
#include "../ReadFile.h"
#include "../LampCore.h"
#include "../Transaction.h"

using namespace std;

/** @file */

/** Read result file of lamp
 */
class Flower_readfile {
public:
	vector <double> motifRpvalue;	/**< list of p-value of motif */
	vector <double> motifApvalue;	/**< list of adjusted p-value of motif */
	vector <int> motifNgenes;		/**< list of number of genes */
	vector <double> motifSscore;	/**< list of score of motif */
	vector <string> motifName;		/**< list of name of motif */
	vector <double> combiRank;		/**< list of rank of combination */
	vector <double> combiRpvalue;	/**< list of p-value of combination */
	vector <double> combiApvalue;	/**< list of adjusted p-value of combination */
	vector <int> combiNgenes;		/**< list of number of combination genes */
	vector <double> combiSscore;	/**< list of score of combination */
	vector<vector<string>> combiName;	/**< list of member of combination */
	double factor;			/**< factor */
	double threshFactor;	/**< threshhold of factor */
	double significance;	/**< significance */

public:
	Flower_readfile();
	virtual ~Flower_readfile();
	double readResult(string resfname, string csvfname, string tabfname);
	
private:
	double calPValue(const vector<Transaction*> &transaction_list,
		map<string, int>& colname2id, const vector<string>& itemset_str_lst,
		int &down_size);
	
	FunctionsSuper* pval_f;
};

#endif /* FLOWER_READFILE_H */

