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
 * File:   PvalTable.h
 * Author: miura
 *
 * Created on 2016/01/19, 18:18
 */

#ifndef PVALTABLE_H
#define PVALTABLE_H

#include <map>
#include <utility>

/** @file */
/** this class is cahe of P-value
 */
class PvalTable {
public:
	PvalTable();
//	virtual ~PvalTable();
	
	/**
	 * get stored value.
	 * @param row index of row.
	 * @param col index of column.
	 * @return stored value.
	 */
	double getValue(int row, int col ) {
		std::pair<int, int> k(row, col);
		if (t.find(k) == t.end())
			return -1;
		else
			return t[k];
	}

	/**
	 * store value
	 * @param row index of row.
	 * @param col index of column
	 * @param pval store value.
	 */
	void putValue(int row, int col, double pval ) {
		std::pair<int, int> k(row, col);
		t[k] = pval;
	}

private:
	std::map<std::pair<int, int>, double> t; /**< to store value */
};

#endif /* PVALTABLE_H */

