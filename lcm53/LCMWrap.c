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
 * LCM wrapeer
 * File:   LCMWrap.c
 * Author: miura
 * 
 * Created on 2016/01/06, 16:23
 */

#ifndef _lcmwrap_c_
#define _lcmwrap_c_

#define _NO_MAIN_

#include "lcm.c"
#include "LCMWrap.h"

/**
 * Run LCM-LAMP and return the optimal minimum support.
 * @param mode
 * @param n1
 * @param p_mode
 * @param input_file
 * @param sig_level
 * @param arity_limit
 * @return 
 */
int LCMWrap_LAMP (char mode, int n1, int p_mode, char* input_file, double sig_level, int arity_limit){
	PROBLEM PP;
	ITEMSET *II = &PP.II;
	TRSACT *TT = &PP.TT;
	
	PROBLEM_init (&PP);

	II->lamp_stat = 1;
	// mode = C
	// option : "C", "-LAMP", n1_str, "-LAMP_P", p_mode_str, input_file, sig_level
	if ( mode == 'C' ){
		PP.problem |= PROBLEM_CLOSED; PP.TT.flag2 |= TRSACT_INTSEC;
        if ( II->flag2 & ITEMSET_LAMP ){
			II->topk.base = n1;
        } else {
			II->flag2 |= ITEMSET_LAMP;
			II->flag |= ITEMSET_SC2;
			II->topk.end = n1;
        }
		// change the statistical test in LAMP mode
		II->lamp_stat = p_mode;
		PP.TT.fname = input_file;
		II->frq_lb = (WEIGHT)sig_level;
		if ( II->flag2 & ITEMSET_LAMP ){ 
			II->th = sig_level; 
			II->frq_lb = 1; 
			II->lamp_alpha = II->th;
			// check the ...
			printf ("II->lamp_stat: %d, II->lamp_alpha: %f\n", 	II->lamp_stat, II->lamp_alpha);
		}
	}
	// mode F
	// option : "F", "-LAMP", n1_str, "-LAMP_P", p_mode_str, "-u", arity_limit, input_file, sig_level
	else if ( mode == 'F' ){
		PP.problem |= PROBLEM_FREQSET; II->flag |= ITEMSET_ALL;
        if ( II->flag2 & ITEMSET_LAMP ){
			II->topk.base = n1;
        } else {
			II->flag2 |= ITEMSET_LAMP;
			II->flag |= ITEMSET_SC2;
			II->topk.end = n1;
        }
		II->lamp_stat = p_mode;
		II->ub = arity_limit;
		PP.TT.fname = input_file;
		II->frq_lb = (WEIGHT)sig_level;
		if ( II->flag2 & ITEMSET_LAMP ){ 
			II->th = sig_level; 
			II->frq_lb = 1; 
			II->lamp_alpha = II->th;
			// check the ...
			printf ("II->lamp_stat: %d, II->lamp_alpha: %f\n", 	II->lamp_stat, II->lamp_alpha);
		}
	}
	
	if ( ERROR_MES ) return (-1);
	TT->flag |= LOAD_PERM +LOAD_DECSORT +LOAD_RM_DUP;
	TT->flag2 |= TRSACT_FRQSORT +TRSACT_MAKE_NEW +TRSACT_DELIV_SC +TRSACT_ALLOC_OCC + ((II->flag & ITEMSET_TRSACT_ID)?0: (TRSACT_SHRINK+TRSACT_1ST_SHRINK));
	if ( II->flag&ITEMSET_RULE ) TT->w_lb = -WEIGHTHUGE; else TT->w_lb = II->frq_lb;
	PP.SG.flag =  LOAD_EDGE;
	PROBLEM_load (&PP);

	int frq = -1;
	if ( !ERROR_MES && TT->T.clms>0 ){
		LCM_init (&PP);
		if ( !ERROR_MES ) LCM (&PP, TT->T.clms, &PP.oo, TT->total_w_org, TT->total_pw_org);
		frq = II->topk_frq;
	}
	
	TT->sc = NULL;
	PROBLEM_end (&PP);
	return frq;
}

/**
 * Construct frequent patterns list.
 * @param mode
 * @param upper_sup
 * @param low_sup
 * @param arity_limit
 * @param input_file
 * @param out_file
 * @return 
 */
int LCMWrap_freq(char mode, double upper_sup, double low_sup, int arity_limit, char* input_file, char* out_file) {
	PROBLEM PP;
	ITEMSET *II = &PP.II;
	TRSACT *TT = &PP.TT;
	
	PROBLEM_init (&PP);

	II->lamp_stat = 1;
	II->flag |= ITEMSET_FREQ;
	II->flag |= ITEMSET_TRSACT_ID;
	
	// mode C
	// option : "CIf", "-U", upper_sup, input_file, low_sup, out_file 
	if ( mode == 'C' ){
		PP.problem |= PROBLEM_CLOSED; PP.TT.flag2 |= TRSACT_INTSEC;
		II->frq_ub = (WEIGHT)upper_sup;
	}
	// mode c
	// option : "CIf", input_file, low_sup, out_file 
	else if ( mode == 'c' ){
		PP.problem |= PROBLEM_CLOSED; PP.TT.flag2 |= TRSACT_INTSEC;
	}
	// mode F
	// option "FIf", "-U", upper_sup, "-u", arity_limit , input_file, low_sup, out_file
	else if ( mode == 'F' ){
		PP.problem |= PROBLEM_FREQSET; II->flag |= ITEMSET_ALL;
		II->frq_ub = (WEIGHT)upper_sup;
		II->ub = arity_limit;
	}
	// mode f
	// option : "FIf", "-u", arity_limit, input_file, low_sup, out_file
	else if ( mode == 'f' ) {
		PP.problem |= PROBLEM_FREQSET; II->flag |= ITEMSET_ALL;
		II->ub = arity_limit;
	}
	
	PP.TT.fname = input_file;
	II->frq_lb = (WEIGHT)low_sup;
	PP.output_fname = out_file;
	if ( II->flag2 & ITEMSET_LAMP ){
		II->th = low_sup;
		II->frq_lb = 1; 
		II->lamp_alpha = II->th;
		// check the ...
		printf ("II->lamp_stat: %d, II->lamp_alpha: %f\n", 	II->lamp_stat, II->lamp_alpha);
	}
	
	if ( ERROR_MES ) return (1);
	TT->flag |= LOAD_PERM +LOAD_DECSORT +LOAD_RM_DUP;
	TT->flag2 |= TRSACT_FRQSORT +TRSACT_MAKE_NEW +TRSACT_DELIV_SC +TRSACT_ALLOC_OCC + ((II->flag & ITEMSET_TRSACT_ID)?0: (TRSACT_SHRINK+TRSACT_1ST_SHRINK));
	if ( II->flag&ITEMSET_RULE ) TT->w_lb = -WEIGHTHUGE; else TT->w_lb = II->frq_lb;
	PP.SG.flag =  LOAD_EDGE;
	PROBLEM_load (&PP);

	if ( !ERROR_MES && TT->T.clms>0 ){
		LCM_init (&PP);
		if ( !ERROR_MES ) LCM (&PP, TT->T.clms, &PP.oo, TT->total_w_org, TT->total_pw_org);
		ITEMSET_last_output (II);
	}
	
	TT->sc = NULL;
	PROBLEM_end (&PP);
	return (ERROR_MES?1:0);
}

#endif // _lcmwrap_c_
