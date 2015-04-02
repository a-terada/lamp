/*  itemset search input/output common routines 
            Apr. 2, 2015  modified by LAMP development team */
/* This code is modified from LCM ver 5.3 
    which is downloaded from http://research.nii.ac.jp/~uno/codes.htm.
   For fast computation of LAMP, 
    ITEMSET_lamp and stdNorDistribution functions are added. 
*/

/* routines for itemset mining */

#ifndef _itemset_c_
#define _itemset_c_

#include"itemset.h"
#include"queue.c"
#include"aheap.c"
#include <math.h>

/* flush the write buffer, available for multi-core mode */
void ITEMSET_flush (ITEMSET *I, FILE2 *fp){
  if ( !(I->flag&ITEMSET_MULTI_OUTPUT) || (fp->buf-fp->buf_org) > FILE2_BUFSIZ/2 ){
    SPIN_LOCK(I->multi_core, I->lock_output);
    FILE2_flush (fp);
    SPIN_UNLOCK(I->multi_core, I->lock_output);
  }
}

/* Output information about ITEMSET structure. flag&1: print frequency constraint */
void ITEMSET_print (ITEMSET *I, int flag){
  if ( I->lb>0 || I->ub<INTHUGE ){
    if ( I->lb > 0 ) print_err ("%d <= ", I->lb);
    print_err ("itemsets ");
    if ( I->ub < INTHUGE ) print_err (" <= %d\n", I->ub);
    print_err ("\n");
  }
  if ( flag&1 ){
    if ( I->frq_lb > -WEIGHTHUGE ) print_err (WEIGHTF" <=", I->frq_lb);
    print_err (" frequency ");
    if ( I->frq_ub < WEIGHTHUGE ) print_err (" <="WEIGHTF, I->frq_ub);
    print_err ("\n");
  }
}

/* ITEMSET initialization */
void ITEMSET_init (ITEMSET *I){
  I->flag = I->flag2 = 0;
  I->progress = 0;
  I->iters = I->iters2 = I->iters3 = 0;
  I->solutions = I->solutions2 = I->max_solutions = I->outputs = I->outputs2 = 0;
  I->topk.end = 0;
  I->item_max = I->item_max_org = 0;
  I->ub = I->len_ub = I->gap_ub = INTHUGE;
  I->lb = I->len_lb = I->gap_lb = 0;
  I->frq = I->pfrq = I->total_weight = 0;
  I->ratio = I->prob = I->th = 0.0;
  I->posi_ub = I->nega_ub = I->frq_ub = WEIGHTHUGE;
  I->posi_lb = I->nega_lb = I->frq_lb = I->setrule_lb = -WEIGHTHUGE;
  I->dir = 0;
  I->target = INTHUGE;
  I->prob_ub = I->ratio_ub = I->rposi_ub = 1;
  I->prob_lb = I->ratio_lb = I->rposi_lb = 0;
  I->itemflag = NULL;
  I->perm = NULL;
  I->item_frq = NULL;
  I->sc = I->sc2 = NULL;
  I->X = NULL;
  I->fp = NULL;
  I->separator = ' ';
  I->topk = INIT_AHEAP;
  I->itemtopk = NULL;
  I->itemtopk_ary = NULL;
  I->itemtopk_item = I->itemtopk_item2 = I->itemtopk_end = 0;
  I->topk_k = I->topk_frq = 0;
  I->topk_sign = I->itemtopk_sign = 1;  // initialization ; max topk
  I->itemset = I->add = INIT_QUEUE;
  I->set_weight = NULL;
  I->set_occ = NULL;

  I->multi_iters = I->multi_iters2 = I->multi_iters3 = NULL;
  I->multi_outputs = I->multi_outputs2 = NULL;
  I->multi_solutions = I->multi_solutions2 = NULL;
  I->multi_fp = NULL;

  I->multi_core = 0;

  I->lamp_stat = 1; // statistical test, 1 -> Fisher, 2 -> Chi square
}


/* second initialization
   topk.end>0 => initialize heap for topk mining */
/* all pointers will be set to 0, but not for */
/* if topK mining, set topk.end to "K" */
void ITEMSET_alloc (ITEMSET *I, char *fname, PERM *perm, QUEUE_INT item_max, size_t item_max_org){
  LONG i, ii;
  size_t siz = (I->flag&ITEMSET_USE_ORG)?item_max_org+2: item_max+2;
  int j;

  I->prob = I->ratio = 1.0;
  I->frq = 0;
  I->perm = perm;
  I->topk.v = NULL;
  QUEUE_alloc (&I->itemset, (QUEUE_ID)siz); I->itemset.end = (QUEUE_ID)siz;
  if ( I->flag&ITEMSET_ADD ) QUEUE_alloc (&I->add, (QUEUE_ID)siz);
  calloc2 (I->sc, siz+2, goto ERR);  
  if ( I->flag&ITEMSET_SC2 ) calloc2 (I->sc2, I->frq_ub+2, goto ERR); // upper bound of frequency
  if ( I->flag2 & ITEMSET_LAMP ) I->topk_frq = I->frq_lb = 1;  // LAMP mode
  if ( I->topk_k > 0 ){  // allocate topk heap
    if (I->flag & ITEMSET_SC2){
      I->frq_lb = 1; I->topk_frq = 0;
      I->sc2[I->topk_frq] = I->topk_k;
    } else {
      AHEAP_alloc (&I->topk, I->topk_k);
      FLOOP (i, 0, I->topk_k) AHEAP_chg (&I->topk, (AHEAP_ID)i, -WEIGHTHUGE);
      I->frq_lb = -WEIGHTHUGE * I->topk_sign;
    }
  }
  if ( I->itemtopk_end > 0 ){ // allocate topkheap for each element
    calloc2 (I->itemtopk, I->itemtopk_end, goto ERR);
    if ( I->itemtopk_item2 > 0 )
        calloc2 (I->itemtopk_ary, I->itemtopk_end, goto ERR); // allocate itemary
    FLOOP (i, 0, I->itemtopk_end){
      AHEAP_alloc (&I->itemtopk[i], I->itemtopk_item);
      if ( I->itemtopk_item2 > 0 )
          calloc2 (I->itemtopk_ary[i], I->itemtopk_item, goto ERR); // allocate each itemary
      FLOOP (ii, 0, I->itemtopk_item)
          AHEAP_chg (&I->itemtopk[i], (AHEAP_ID)ii, -WEIGHTHUGE);
    }
  }
  
  if ( I->flag&ITEMSET_SET_RULE ){
    calloc2 (I->set_weight, siz, goto ERR);
    if ( I->flag&(ITEMSET_TRSACT_ID+ITEMSET_MULTI_OCC_PRINT) )
        calloc2 (I->set_occ, siz, goto ERR);
  }
  I->iters = I->iters2 = I->solutions = 0;
  I->item_max = item_max;
  I->item_max_org = (QUEUE_INT)item_max_org;
  if ( fname ){
#ifdef _FILE2_LOAD_FROM_MEMORY_
    I->fp = (FILE *)1;
#else 
    if ( strcmp (fname, "-") == 0 ) I->fp = stdout;
    else fopen2 (I->fp, fname, (I->flag&ITEMSET_APPEND)?"a":"w", goto ERR);
#endif
  } else I->fp = 0;
  if ( I->flag&ITEMSET_ITEMFRQ )
    malloc2 (I->item_frq, item_max+2, goto ERR);
  if ( I->flag&ITEMSET_RULE ){
    calloc2 (I->itemflag, item_max+2, goto ERR);
  }
  I->total_weight = 1;
  j = MAX(I->multi_core,1);
  calloc2 (I->multi_iters, j*7, goto ERR);
  I->multi_iters2 = I->multi_iters + j;
  I->multi_iters3 = I->multi_iters2 + j;
  I->multi_outputs = I->multi_iters3 + j;
  I->multi_outputs2 = I->multi_outputs + j;
  I->multi_solutions = I->multi_outputs2 + j;
  I->multi_solutions2 = I->multi_solutions + j;
  
  calloc2 (I->multi_fp, j, goto ERR);
  FLOOP (i, 0, j)
      FILE2_open_ (I->multi_fp[i], I->fp, goto ERR);
#ifdef MULTI_CORE
  if ( I->multi_core > 0 ){
    pthread_spin_init (&I->lock_counter, PTHREAD_PROCESS_PRIVATE);
    pthread_spin_init (&I->lock_sc, PTHREAD_PROCESS_PRIVATE);
    pthread_spin_init (&I->lock_output, PTHREAD_PROCESS_PRIVATE);
  }
#endif
  return;
  ERR:;
  ITEMSET_end (I);
  EXIT;
}

/* sum the counters computed by each thread */
void ITEMSET_merge_counters (ITEMSET *I){
  int i;
  FLOOP (i, 0, MAX(I->multi_core,1)){
    I->iters += I->multi_iters[i];
    I->iters2 += I->multi_iters2[i];
    I->iters3 += I->multi_iters3[i];
    I->outputs += I->multi_outputs[i];
    I->outputs2 += I->multi_outputs2[i];
    I->solutions += I->multi_solutions[i];
    I->solutions2 += I->multi_solutions2[i];
    if ( I->multi_fp[i].buf ) FILE2_flush_last (&I->multi_fp[i]);
  }
  ARY_FILL (I->multi_iters, 0, MAX(I->multi_core,1)*7, 0);
}

/*******************************************************************/
/* termination of ITEMSET */
/*******************************************************************/
void ITEMSET_end (ITEMSET *I){
  LONG i;
  QUEUE_end (&I->itemset);
  QUEUE_end (&I->add);
  AHEAP_end (&I->topk);
  FLOOP (i, 0, I->itemtopk_end){
    AHEAP_end (&I->itemtopk[i]);
    if ( I->itemtopk_ary ) free2 (I->itemtopk_ary[i]);
  }

#ifndef _FILE2_LOAD_FROM_MEMORY_
  fclose2 (I->fp);
#endif
  mfree (I->sc, I->sc2, I->item_frq, I->itemflag, I->perm, I->set_weight, I->set_occ, I->itemtopk_ary);

  if ( I->multi_fp )
      FLOOP (i, 0, MAX(I->multi_core,1)) free2 (I->multi_fp[i].buf_org);
  mfree (I->multi_iters, I->multi_fp);
#ifdef MULTI_CORE
  if ( I->multi_core>0 ){
    pthread_spin_destroy(&I->lock_counter);
    pthread_spin_destroy(&I->lock_sc);
    pthread_spin_destroy(&I->lock_output);
  }
#endif
  ITEMSET_init (I);
}

/*******************************************************************/
/* output at the termination of the algorithm */
/* print #of itemsets of size k, for each k */
/*******************************************************************/
void ITEMSET_last_output (ITEMSET *I){
  QUEUE_ID i;
  LONG n=0, nn=0;
  WEIGHT w;

  ITEMSET_merge_counters (I);
  if ( !(I->flag&SHOW_MESSAGE) ) return;  // "no message" is specified

  if ( I->flag2 & ITEMSET_LAMP ){
    printf ("frq= %lld ,#sol.= %lld\n", I->topk_frq, I->topk_k);
    print_err ("iters=" LONGF, I->iters);
    if ( I->flag&ITEMSET_ITERS2 ) print_err (", iters2=" LONGF, I->iters2);
    print_err ("\n");
    return;
  }

  if ( I->itemtopk_end > 0 ){  // output values of the kth-best solution for each item
    FLOOP (n, 0, I->itemtopk_end){
      FLOOP (nn, 0, I->itemtopk[n].end){
        i = AHEAP_findmin_head (&I->itemtopk[n]);
        w = AHEAP_H (I->itemtopk[n], i);
        if ( w == -WEIGHTHUGE*I->itemtopk_sign ) break;
        if ( I->itemtopk_ary ) printf ("%d ", I->itemtopk_ary[n][i]);
        fprint_WEIGHT (stdout, w);
        printf (" ");
        AHEAP_chg (&(I->itemtopk[n]), i, WEIGHTHUGE);
      }
      printf ("\n");
    }
    goto END;
  }

  if ( I->topk_k > 0 ){  // output value of the kth-best solution
    if ( I->topk.v ){
      i = AHEAP_findmin_head (&I->topk);
      fprint_WEIGHT (stdout, AHEAP_H (I->topk, i)*I->topk_sign);
    } else fprintf (stdout, LONGF, I->topk_frq);
    printf ("\n");

    goto END;
  }
  FLOOP (i, 0, I->itemset.end+1){
    n += I->sc[i];
    if ( I->sc[i] != 0 ) nn = i;
  }
  if ( n!=0 ){
    printf (LONGF "\n", n);
    FLOOP (i, 0, nn+1) printf (LONGF "\n", I->sc[i]);
  }
  
  END:;
  print_err ("iters=" LONGF, I->iters);
  if ( I->flag&ITEMSET_ITERS2 ) print_err (", iters2=" LONGF, I->iters2);
  print_err ("\n");
  
  if (I->flag & ITEMSET_SC2){ // count by frequency
    FLOOP (i, 0, I->frq_ub+1){
      if ( I->sc2[i] != 0 ) printf (QUEUE_INTF "," LONGF "\n", i, I->sc2[i]);
    }
  }
}

/* output frequency, coverage */
void ITEMSET_output_frequency (ITEMSET *I, int core_id){
  FILE2 *fp = &I->multi_fp[core_id];
  if ( I->flag&(ITEMSET_FREQ+ITEMSET_PRE_FREQ) ){
    if ( I->flag&ITEMSET_FREQ ) FILE2_putc (fp, ' ');
    FILE2_print_WEIGHT (fp, I->frq, 4, '(');
    FILE2_putc (fp, ')');
    if ( I->flag&ITEMSET_PRE_FREQ ) FILE2_putc (fp, ' ');
  }
  if ( I->flag&ITEMSET_OUTPUT_POSINEGA ){ // output positive sum, negative sum in the occurrence
    FILE2_putc (fp, ' ');
    FILE2_print_WEIGHT (fp, I->pfrq, 4, '(');
    FILE2_print_WEIGHT (fp, I->pfrq-I->frq, 4, ',');
    FILE2_print_WEIGHT (fp, I->pfrq/(2*I->pfrq-I->frq), 4, ',');
    FILE2_putc (fp, ')');
  }
}

/**
 * Calculate probability of standard normal distribution.
 * this function returns the probability of one-sided test.
 * x:  
 **/
double stdNorDistribution ( double x ){
  double pi2 = 0.398942280401432677940; 
  int is_value = -1; 
  double y = fabs(x);
  double c = y*y;
  double p = 0.0; 
  double z = exp(-c*0.5)*pi2; 
  if (y < 2.5){
	int i = -1;
	for (i = 20; i > 0; i--){
		p = i*c/(i*2+1+is_value*p);
		is_value = -is_value;
	}
	p = 0.5-z*y/(1.0-p);
  }
  else{
	int i = -1;
	for (i = 20; i > 0; i--){
	  p = i/(y+p); 
	}
	p = z/(y+p);
  }
  //  p = 2. * p; // double p-value because returens about two-sided test.
  return p;
}

/**
 * int p_mode: 1 -> Fisher's exact test, greater
 *           : 2 -> Chi-square test, greater
 */
// topk.end: #records, topk.base: #positive records, PP.th: \alpha, topk_k: #patterns found
void ITEMSET_lamp (ITEMSET *I, LONG s){
  //printf("I->lamp_stat: %d\n", I->lamp_stat);
  if ( I->frq >= I->topk_frq ){ // LAMP  histogram version
	int base0 = I->topk.base - I->topk.end;
	// topk_k: frequency, frq_lb: minimum support, th: alpha/f(lambd) (the upper bound for the frequency)
    I->topk_k += s;  // topk_k: frequency
    while ( I->topk_k >= I->th ){
      I->topk_k -= I->sc2[I->topk_frq]; I->sc2[I->topk_frq] = 0;
	  printf ("frq_lb: %d, topk_k: %lld,  th(%d): %f ->", (int)I->frq_lb, I->topk_k, (int)I->frq_lb - 1, I->th);
	  // update lower bound for the frequency.
	  // p_mode: 1 -> Fisher's exact test (greater)
	  if ( I->lamp_stat == 1 ){
		I->th = I->th * (I->topk.base - I->topk_frq + 1) / (I->topk.end - I->topk_frq + 1);
	  }
	  // p_mode: 2 -> Chi-square test (greater)
	  else{
		int nonfrq = I->topk.base - I->topk_frq;
		double means[2][2] = {{(double)(I->topk_frq) * I->topk.end/I->topk.base, 
							   (double)(I->topk_frq) * base0/I->topk.base}, 
							 {(double)(nonfrq) * I->topk.end/I->topk.base, 
							  (double)(nonfrq) * base0/I->topk.base}};
		//printf ("Mean: %f, %f, %f, %f\n", means[0][0], means[0][1], means[1][0], means[1][1] );
		double chi = pow( fabs((I->topk_frq) - means[0][0]) - YATE_CORR, 2.0)/means[0][0];
		chi += pow( fabs(0 - means[0][1]) - YATE_CORR, 2.0)/means[0][1];
		chi += pow( fabs(I->topk.end - (I->topk_frq) - means[1][0]) - YATE_CORR, 2.0)/means[1][0];
		chi += pow( fabs(base0 - means[1][1]) - YATE_CORR, 2.0)/means[1][1];
		double pval = 1.0;
		if (chi != 0.0){
		  pval = stdNorDistribution( sqrt( chi ) );
		}
		//printf ("chi^2: %f, p-value: %e, ", chi, pval );
		I->th = I->lamp_alpha / pval; 
	  }
	  printf ("th: %f\n", I->th);
	  I->topk_frq++; 
      I->frq_lb = I->topk_frq;
      if ( I->topk_frq == I->topk.end ) I->frq_lb = I->topk.base+1;
	  printf ("frq_lb: %f\n", I->frq_lb);
    }
  }
  return;
}


#ifdef _trsact_h_
void ITEMSET_output_occ (ITEMSET *I, QUEUE *occ, int core_id){
  QUEUE_ID i;
  QUEUE_INT *x;
  FILE2 *fp = &I->multi_fp[core_id];
  TRSACT *TT = (TRSACT *)(I->X);
  VEC_ID j, ee = TT->rows_org;
  int flag = I->flag&(ITEMSET_TRSACT_ID+ITEMSET_MULTI_OCC_PRINT);

  i=0; MQUE_FLOOP_ (*occ, x, TT->occ_unit){
    if ( (I->flag&ITEMSET_RM_DUP_TRSACT)==0 || *x != ee ){
      FILE2_print_int (fp, TT->trperm? TT->trperm[*x]: *x, I->separator);
      if (flag == ITEMSET_MULTI_OCC_PRINT ){
        FLOOP (j, 1, (VEC_ID)(TT->occ_unit/sizeof(QUEUE_INT)))
            FILE2_print_int (fp, *(x+j), I->separator);
      } else if ( flag == (ITEMSET_MULTI_OCC_PRINT+ITEMSET_TRSACT_ID) ){
         FILE2_print_int (fp, *(x+1), I->separator);
      }
    }
    ee = *x;
    if ( (++i)%256==0 ) ITEMSET_flush (I, fp);
  }
#ifdef _FILE2_LOAD_FROM_MEMORY_
  *((int *)__write_to_memory__) = INTHUGE;
  __write_to_memory__ = (char *)(((int *)__write_to_memory__) + 1);
#else
  FILE2_putc (fp, '\n');
#endif
}
#endif

/* output an itemset to the output file */
void ITEMSET_output_itemset (ITEMSET *I, QUEUE *occ, int core_id){
  QUEUE_ID i;
  QUEUE_INT e;
#ifdef _agraph_h_
  QUEUE_INT ee;
#endif

  FILE2 *fp = &I->multi_fp[core_id];
  
  I->multi_outputs[core_id]++;
  if ( (I->flag&SHOW_PROGRESS ) && (I->multi_outputs[core_id]%(ITEMSET_INTERVAL) == 0) )
      print_err ("---- " LONGF " solutions in " LONGF " candidates\n",
                  I->multi_solutions[core_id], I->multi_outputs[core_id]);
  if ( I->itemset.t < I->lb || I->itemset.t > I->ub ) return;
  if ( (I->flag&ITEMSET_IGNORE_BOUND)==0 && (I->frq < I->frq_lb || I->frq > I->frq_ub) ) return;
  if ( (I->flag&ITEMSET_IGNORE_BOUND)==0 && (I->pfrq < I->posi_lb || I->pfrq > I->posi_ub || (I->frq - I->pfrq) > I->nega_ub || (I->frq - I->pfrq) < I->nega_lb) ) return;

  I->multi_solutions[core_id]++;
  if ( I->max_solutions>0 && I->multi_solutions[core_id] > I->max_solutions ){
    ITEMSET_last_output (I);
    ERROR_MES = "reached to maximum number of solutions";
    EXIT;
  }

  I->sc[I->itemset.t]++;
  if (I->flag & ITEMSET_SC2) I->sc2[(QUEUE_INT)I->frq]++;  // LAMP mode

  if ( I->flag2 & ITEMSET_LAMP ){ ITEMSET_lamp (I, 1); return; }
  if ( I->itemtopk_end > 0 ){
    e = AHEAP_findmin_head (&(I->itemtopk[I->itemtopk_item]));
    if ( I->frq > AHEAP_H (I->itemtopk[I->itemtopk_item], e) ){
      AHEAP_chg (&(I->itemtopk[I->itemtopk_item]), e, I->frq * I->itemtopk_sign);
      if ( I->itemtopk_ary ) I->itemtopk_ary[I->itemtopk_item][e] = I->itemtopk_item2;
    }
    return;
  }

  if ( I->topk_k > 0 ){
    if ( I->topk.v ){
      e = AHEAP_findmin_head (&(I->topk));
      if ( I->frq * I->topk_sign > AHEAP_H (I->topk, e) ){
        AHEAP_chg (&(I->topk), e, I->frq * I->topk_sign);
        e = AHEAP_findmin_head (&(I->topk));
        I->frq_lb = AHEAP_H (I->topk, e) * I->topk_sign;
      }
    } else {  // histogram version
      if ( I->frq > I->topk_frq ){
        I->sc2[I->topk_frq]--;
        while (I->sc2[I->topk_frq]==0) I->topk_frq++;
        I->frq_lb = I->topk_frq+1;
      }
    }
    return;
  }
  
  if ( I->fp ){
    if ( I->flag&ITEMSET_PRE_FREQ ) ITEMSET_output_frequency (I, core_id);
    if ( (I->flag & ITEMSET_NOT_ITEMSET) == 0 ){
#ifdef _agraph_h_
      if ( I->flag&ITEMSET_OUTPUT_EDGE ){
        FLOOP (i, 0, I->itemset.t){
          e = I->itemset.v[i];
          ee = AGRAPH_INC_FROM(*((AGRAPH *)(I->X)), e, I->dir);
          FILE2_print_int (fp, I->perm? I->perm[ee]: ee, '(' );
          ee = AGRAPH_INC_TO(*((AGRAPH *)(I->X)), e, I->dir);
          FILE2_print_int (fp, I->perm? I->perm[ee]: ee, I->separator);
#ifdef _FILE2_LOAD_FROM_MEMORY_
          FILE2_putc (fp, ')');
#endif
          if ( i<I->itemset.t-1 ) FILE2_putc (fp, I->separator);
          if ( (i+1)%256==0 ) ITEMSET_flush (I, fp);
        }
        goto NEXT;
      }
#endif
      FLOOP (i, 0, I->itemset.t){
        e = I->itemset.v[i];
        FILE2_print_int (fp,  I->perm? I->perm[e]: e, i==0? 0: I->separator);
        if ( (i+1)%256==0 ) ITEMSET_flush (I, fp);
      }
#ifdef _agraph_h_
      NEXT:;
#endif
    }
    if ( !(I->flag&ITEMSET_PRE_FREQ) ) ITEMSET_output_frequency (I, core_id);
    if ( ((I->flag & ITEMSET_NOT_ITEMSET) == 0) || (I->flag&ITEMSET_FREQ) || (I->flag&ITEMSET_PRE_FREQ) ){
#ifdef _FILE2_LOAD_FROM_MEMORY_
  FILE2_WRITE_MEMORY (QUEUE_INT, FILE2_LOAD_FROM_MEMORY_END);
#else
      FILE2_putc (fp, '\n');
#endif
    }
#ifdef _trsact_h_
    if (I->flag&(ITEMSET_TRSACT_ID+ITEMSET_MULTI_OCC_PRINT)) ITEMSET_output_occ (I, occ, core_id);
#endif
    ITEMSET_flush (I, fp);
  }
}

/* output itemsets with adding all combination of "add"
   at the first call, i has to be "add->t" */
void ITEMSET_solution_iter (ITEMSET *I, QUEUE *occ, int core_id){
  QUEUE_ID t=I->add.t;
  if ( I->itemset.t > I->ub ) return;
  ITEMSET_output_itemset (I, occ, core_id);
if ( ERROR_MES ) return;
  BLOOP (I->add.t, I->add.t, 0){
    QUE_INS (I->itemset, I->add.v[I->add.t]);
    ITEMSET_solution_iter (I, occ, core_id);
if ( ERROR_MES ) return;
    I->itemset.t--;
  }
  I->add.t = t;
}

void ITEMSET_solution (ITEMSET *I, QUEUE *occ, int core_id){
  QUEUE_ID i;
  LONG s;
  if ( I->itemset.t > I->ub ) return;
  if ( I->flag & ITEMSET_ALL ){
    if ( I->fp || I->topk.v ) ITEMSET_solution_iter (I, occ, core_id);
    else {
      s=1; FLOOP (i, 0, I->add.t+1){
        I->sc[I->itemset.t+i] += s;
        s = s*(I->add.t-i)/(i+1);
      }
      if (I->flag & ITEMSET_SC2){
        s = 1<<I->add.t;
        I->sc2[(QUEUE_INT)I->frq] += s;
        if ( I->flag2 & ITEMSET_LAMP ) ITEMSET_lamp (I, s);  // LAMP mode
        else if ( I->topk_k > 0 && I->frq > I->topk_frq ){ // top-k histogram version
          while (1){
            if ( I->sc2[I->topk_frq] > s ){ I->sc2[I->topk_frq] -= s; break; }
            s -= I->sc2[I->topk_frq];
            I->sc2[I->topk_frq++] = 0; 
          }
          I->frq_lb = I->topk_frq+1;
        }
      }
    }
  } else {
    FLOOP (i, 0, I->add.t) QUE_INS (I->itemset, I->add.v[i]);
    ITEMSET_output_itemset (I, occ, core_id);
    I->itemset.t -= I->add.t;
  }
}

/*************************************************************************/
/* ourput a rule */
/*************************************************************************/
void ITEMSET_output_rule (ITEMSET *I, QUEUE *occ, double p1, double p2, size_t item, int core_id){
  FILE2 *fp = &I->multi_fp[core_id];
  if ( fp->fp && !(I->topk.v) ){
    FILE2_print_real (fp, p1, 4, '(');
    FILE2_print_real (fp, p2, 4, ',');
    FILE2_putc (fp, ')');
    FILE2_print_int (fp, I->perm[item], I->separator);
    FILE2_puts (fp, " <= ");
  }
  if ( I->flag & ITEMSET_RULE ){
    if ( I->flag & ITEMSET_RULE_ADD ) ITEMSET_solution (I, occ, core_id);
    else ITEMSET_output_itemset (I, occ, core_id);
  } else ITEMSET_solution (I, occ, core_id);
}
/*************************************************************************/
/* check all rules for a pair of itemset and item */
/*************************************************************************/
void ITEMSET_check_rule (ITEMSET *I, WEIGHT *w, QUEUE *occ, size_t item, int core_id){
  double p = w[item]/I->frq, pp, ff;
//  printf ("[ratio] %f, p=%f, (%f/ %f), %d(%d) <= ", I->ratio_lb, p, w[item], I->frq, I->perm[item], I->itemflag[item]);
  if ( I->itemflag[item]==1 ) return;
  if ( w[item] <= -WEIGHTHUGE ) p = 0;
  pp = p; ff = I->item_frq[item];
  if ( I->flag & ITEMSET_RULE_SUPP ){ pp = w[item]; ff *= I->total_weight; }

  if ( I->flag & (ITEMSET_RULE_FRQ+ITEMSET_RULE_INFRQ)){
    if ( (I->flag & ITEMSET_RULE_FRQ) && p < I->ratio_lb ) return;
    if ( (I->flag & ITEMSET_RULE_INFRQ) && p > I->ratio_ub ) return;
    ITEMSET_output_rule (I, occ, p, ff, item, core_id);
  } else if ( I->flag & (ITEMSET_RULE_RFRQ+ITEMSET_RULE_RINFRQ) ){
    if ( (I->flag & ITEMSET_RULE_RFRQ) && (1-p) > I->ratio_lb * (1-I->item_frq[item]) ) return;
    if ( (I->flag & ITEMSET_RULE_RINFRQ) && p > I->ratio_ub * I->item_frq[item] ) return;
    ITEMSET_output_rule (I, occ, pp, ff, item, core_id);
  }
}

/*************************************************************************/
/* check all rules for an itemset and all items */
/*************************************************************************/
void ITEMSET_check_all_rule (ITEMSET *I, WEIGHT *w, QUEUE *occ, QUEUE *jump, WEIGHT total, int core_id){
  QUEUE_ID i, t;
  QUEUE_INT e, f=0, *x;
  WEIGHT d = I->frq/total;

    // checking out of range for itemset size and (posi/nega) frequency
  if ( I->itemset.t+I->add.t < I->lb || I->itemset.t>I->ub || (!(I->flag&ITEMSET_ALL) && I->itemset.t+I->add.t>I->ub)) return;
  if ( !(I->flag&ITEMSET_IGNORE_BOUND) && (I->frq < I->frq_lb || I->frq > I->frq_ub) ) return;
  if ( !(I->flag&ITEMSET_IGNORE_BOUND) && (I->pfrq < I->posi_lb || I->pfrq > I->posi_ub || (I->frq - I->pfrq) > I->nega_ub || (I->frq - I->pfrq) < I->nega_lb) ) return;

  if ( I->flag&ITEMSET_SET_RULE ){  // itemset->itemset rule for sequence mining
    FLOOP (i, 0, I->itemset.t-1){
      if ( I->frq/I->set_weight[i] >= I->setrule_lb && I->fp ){
        I->sc[i]++;
        if (I->flag & ITEMSET_SC2) I->sc2[(QUEUE_INT)I->frq]++;
        if ( I->flag2 & ITEMSET_LAMP ) ITEMSET_lamp (I, 1);  // LAMP mode
        if ( I->flag&ITEMSET_PRE_FREQ ) ITEMSET_output_frequency (I, core_id);
        FLOOP (t, 0, I->itemset.t){
          FILE2_print_int (&I->multi_fp[core_id], I->itemset.v[t], t?I->separator:0);
          if ( t == i ){
            FILE2_putc (&I->multi_fp[core_id], ' ');
            FILE2_putc (&I->multi_fp[core_id], '=');
            FILE2_putc (&I->multi_fp[core_id], '>');
          }
        }
        if ( !(I->flag&ITEMSET_PRE_FREQ) ) ITEMSET_output_frequency ( I, core_id);
        FILE2_putc (&I->multi_fp[core_id], ' ');
        FILE2_print_real (&I->multi_fp[core_id], I->frq/I->set_weight[i], 4, '(');
        FILE2_putc (&I->multi_fp[core_id], ')');
#ifdef _FILE2_LOAD_FROM_MEMORY_
  FILE2_WRITE_MEMORY (QUEUE_INT, FILE2_LOAD_FROM_MEMORY_END);
#else
        FILE2_putc (&I->multi_fp[core_id], '\n');
#endif
#ifdef _trsact_h_
        if ( I->flag&(ITEMSET_TRSACT_ID+ITEMSET_MULTI_OCC_PRINT) ){
            ITEMSET_output_occ (I, I->set_occ[i], core_id);
        }
#endif
        ITEMSET_flush (I, &I->multi_fp[core_id]);
      }
    }
  }
    // constraint of relational frequency
  if ( ((I->flag&ITEMSET_RFRQ)==0 || d >= I->prob_lb * I->prob ) 
      && ((I->flag&ITEMSET_RINFRQ)==0 || d <= I->prob * I->prob_ub) ){
    if ( I->flag&ITEMSET_RULE ){  //  rule mining routines
      if ( I->itemset.t == 0 ) return;
      if ( I->target < I->item_max ){
        MQUE_FLOOP (*jump, x){
          if ( *x == I->target ){ 
              ITEMSET_check_rule (I, w, occ, *x, core_id);   if (ERROR_MES) return;
          }
        }
//        ITEMSET_check_rule (I, w, occ, I->target, core_id);    if (ERROR_MES) return;
      } else {
        if ( I->flag & (ITEMSET_RULE_FRQ + ITEMSET_RULE_RFRQ) ){
          if ( I->add.t>0 ){
//            if ( I->itemflag[I->add.v[0]] ) // for POSI_EQUISUPP (occ_w[e] may not be 100%, in the case)
            f = I->add.v[I->add.t-1]; t = I->add.t; I->add.t--;
            FLOOP (i, 0, t){
              e = I->add.v[i];
              I->add.v[i] = f;
              ITEMSET_check_rule (I, w, occ, e, core_id);    if (ERROR_MES) return;
              I->add.v[i] = e;
            }
            I->add.t++;
          }
          MQUE_FLOOP (*jump, x)
              ITEMSET_check_rule (I, w, occ, *x, core_id);   if (ERROR_MES) return;
        } else {
          if ( I->flag & (ITEMSET_RULE_INFRQ + ITEMSET_RULE_RINFRQ) ){
//          ARY_FLOOP ( *jump, i, e ) I->itemflag[e]--;
            FLOOP (i, 0, I->item_max){
              if ( I->itemflag[i] != 1 ){
                ITEMSET_check_rule (I, w, occ, i, core_id);     if (ERROR_MES) return;
              }
            }
//          ARY_FLOOP ( *jump, i, e ) I->itemflag[e]++;
//        } 
//        ARY_FLOOP ( *jump, i, e ) ITEMSET_check_rule (I, w, occ, e);
          }
        }
      }
    } else {  // usual mining (not rule mining)
      if ( I->fp && (I->flag&(ITEMSET_RFRQ+ITEMSET_RINFRQ))){
        FILE2_print_real (&I->multi_fp[core_id], d, 4, '[');
        FILE2_print_real (&I->multi_fp[core_id], I->prob, 4, ',');
        FILE2_putc (&I->multi_fp[core_id], ']');
      }
      ITEMSET_solution (I, occ, core_id);
    }
  }
}

#endif
