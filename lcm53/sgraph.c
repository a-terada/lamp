/*  graph library by array list
            12/Feb/2002    by Takeaki Uno
    homepage:   http://research.nii.ac.jp/~uno/index.html  */
/* This program is available for only academic use, basically.
   Anyone can modify this program, but he/she has to write down 
    the change of the modification on the top of the source code.
   Neither contact nor appointment to Takeaki Uno is needed.
   If one wants to re-distribute this code, please
    refer the newest code, and show the link to homepage of 
    Takeaki Uno, to notify the news about the codes for the users. */

#ifndef _sgraph_c_
#define _sgraph_c_

#include"sgraph.h"
#include"vec.c"

SGRAPH INIT_SGRAPH = {TYPE_SGRAPH,NULL,0,INIT_SETFAMILY_,INIT_SETFAMILY_,INIT_SETFAMILY_,0,NULL,NULL,NULL,NULL};

/*  initialization  */
void SGRAPH_alloc (SGRAPH *G, QUEUE_ID nodes, size_t edge_num, size_t arc_num){
  if ( edge_num > 0 ){
    SETFAMILY_alloc (&G->edge, nodes, NULL, nodes, edge_num);
    if ( G->flag&LOAD_EDGEW && (!ERROR_MES) ) SETFAMILY_alloc_weight (&G->edge, NULL);
  }
  if ( arc_num > 0 ){
    SETFAMILY_alloc (&G->in, nodes, NULL, nodes, arc_num);
    SETFAMILY_alloc (&G->out, nodes, NULL, nodes, arc_num);
    if ( G->flag&LOAD_EDGEW && (!ERROR_MES) ){
      SETFAMILY_alloc_weight (&G->in, NULL);
      SETFAMILY_alloc_weight (&G->out, NULL);
    }
  }
  if (G->flag&LOAD_NODEW) calloc2 (G->node_w, nodes, G->node_w=0);
  if ( ERROR_MES ){ SGRAPH_end (G); EXIT; }
}

/* copy graph G to graph G2. Underconstruction */
//void SGRAPH_cpy (SGRAPH *G2, SGRAPH *G){}

/* free graph object  */
void SGRAPH_end (SGRAPH *G){
  SETFAMILY_end (&G->edge);
  SETFAMILY_end (&G->in);
  SETFAMILY_end (&G->out);
  mfree (G->wbuf, G->perm);
  *G = INIT_SGRAPH;
}


/*  make an edge between u and v.
   If they are already connected, it will be a multiple edge */
void SGRAPH_edge_mk (SGRAPH *G, QUEUE_INT u, QUEUE_INT v, WEIGHT w){
  if ( G->edge.w ){
    G->edge.w[u][G->edge.v[u].t] = w;
    G->edge.w[v][G->edge.v[v].t] = w;
  }
  QUE_INS (G->edge.v[u], v);
  QUE_INS (G->edge.v[v], u);
  G->edge.eles += 2;
}

/*  make an arc between u and v.
   If they are already connected, it will be a multiple arc */
void SGRAPH_arc_mk (SGRAPH *G, QUEUE_INT u, QUEUE_INT v, WEIGHT w){
  if ( G->out.w ) G->out.w[u][G->out.v[u].t] = w;
  if ( G->in.w ) G->in.w[v][G->in.v[v].t] = w;
  QUE_INS (G->out.v[u], v);
  QUE_INS (G->in.v[v], u);
  G->in.eles++;
  G->out.eles++;
}

/* Delete the edge connecting u and v. If edge (u,v) does not exist, nothing will occur. */
void SGRAPH_edge_rm_iter (SETFAMILY *M, QUEUE_INT u, QUEUE_INT v){
  QUEUE_INT i;
  if ( (i = (QUEUE_INT)QUEUE_ele (&M->v[u], v)) >= 0 ){
    QUEUE_rm (&M->v[u], i);
    if ( M->w ) M->w[u][i] = M->w[u][M->v[u].t];
    M->eles--;
  }
}

/* Delete the edge connecting u and v. If edge (u,v) does not exist, nothing will occur. */
void SGRAPH_edge_rm (SGRAPH *G, QUEUE_INT u, QUEUE_INT v){
  SGRAPH_edge_rm_iter (&G->edge, u, v);
  SGRAPH_edge_rm_iter (&G->edge, v, u);
}

/* Delete the arc connecting u and v. If arc (u,v) does not exist, nothing will occur. */
void SGRAPH_arc_rm (SGRAPH *G, QUEUE_INT u, QUEUE_INT v){
  SGRAPH_edge_rm_iter (&G->out, u, v);
  SGRAPH_edge_rm_iter (&G->in, v, u);
}

/*  print graph by numbers  */
void SGRAPH_print (FILE *fp, SGRAPH *G){
  VEC_ID i, j;
  QUEUE_INT e;
  
  fprintf (fp, "#node "VEC_IDF" ,#edge %zd ,#arc %zd\n", SGRAPH_NODE_NUM(*G), G->edge.eles, G->in.eles);
  FLOOP (i, 0, SGRAPH_NODE_NUM(*G)){
    fprintf (fp, "NODE "VEC_IDF" ", i);
    if ( G->node_w ){ fputc ('(', fp); print_WEIGHT (G->node_w[i]); fputc (')', fp); }
    fprintf (fp, " >>\n");
    if ( G->edge.v && G->edge.v[i].t ){
      fprintf (fp, "    edge      : ");
      for (j=0; j<G->edge.v[i].t ; j++){
        e = G->edge.v[i].v[j];
        fprintf (fp, VEC_IDF, e);
        if ( G->edge.w ){ fputc ('(', fp); print_WEIGHT (G->edge.w[i][j]); fputc (')', fp); }
        fputc (',', fp);
      }
      fputc ('\n', fp);
    }
    if ( G->in.v ){
      if ( G->in.v[i].t ){
        fprintf (fp, "    in-arc      : ");
        for (j=0; j<G->in.v[i].t ; j++){
          e = G->in.v[i].v[j];
          fprintf (fp, VEC_IDF, e);
          if ( G->in.w ){ fputc ('(', fp); print_WEIGHT (G->in.w[i][j]); fputc (')', fp); }
          fputc (',', fp);
        }
        fputc ('\n', fp);
      }
    }
    if ( G->out.v ){
      if ( G->out.v[i].t ){
        fprintf (fp, "    out-arc      : ");
        for (j=0; j<G->out.v[i].t ; j++){
          e = G->out.v[i].v[j];
          fprintf (fp, VEC_IDF, e);
          if ( G->out.w ){ fputc ('(', fp); print_WEIGHT (G->out.w[i][j]); fputc (')', fp);}
          fputc (',', fp);
        }
        fputc ('\n', fp);
      }
    }
  }
}

/* Output a graph to file
  Vertices, edges, arcs less than node_num, edge_num, arc_num are written to the file. Input parameters are
  (graph) (file name) (flag)
  SGRAPH_READ_NODEW 512 // read node weight
  SGRAPH_READ_EDGEW 1024 // read edge weight
*/
/*
  format of file:(including notifications to make input file)
   
  the ith row corresponds to node i-1, and
    ID list of nodes adjacent to i, and having ID > i, for undirected graph
    ID list of nodes adjacent to i by out-going arc of i, for directed graph
   Separator is ",", but graph load routine accepts any letter for 
    separator but not a number.
   If the graph has both edges and arcs, write them in two lines separately,
    so a node then uses two lines, and #nodes = #lines/2.
  
    ==  Notifications to make input file ==
   Notice that if 0th line has node 2, and the 2nd line has 0, then there
    will be multiple edge (0,2) and (2,0).
   The read routine does not make error with multiple edges, it is allowed.

   The ID of nodes begin from 0. After reading graph, node_num is set to
    node_end.

   Input file example, without weights, E={(0,1),(0,2),(1,1),(1,3),(2,3)}
===========
   1,2
   1 3
   3
   
   [EOF]
=========
   Nodes are 0,1, and 2, both edges and arcs exist, with node/edge/arc weights)
   5000,1,30
   0,50,1,20,
   100,1,3
   2,20
   200
   
   [EOF]
=======
   where node weights are 5000, 100, 200, and edges and their weights are
    (0,1),30,   (1,1),3
    arcs and their weights are (0,0),50,   (0,1), 20,   (1,2), 20

    In the case of bipartite graph, write the adjacent-node lists only for 
     the node in node set one.
     
    
*/

/* graph load routine. Allocate memory as much as the size of input file.
   parameters are, 
   (graph) (file name) 
 LOAD_EDGE // read undirected edge from file
 LOAD_ARC // read directed arc from file
 LOAD_BIPARTITE // load bipartite graph
 LOAD_NODEW // read node weight
 LOAD_EDGEW // read edge weight
*/
/* In the bipartite case, even if the IDs of node set 2 begin from 0, i.e.,
   overlaps with node 1, the routine automatically correct them. */
/* Directed bipartite graph, all arcs are considered to be from node set 1
 to node set 2. If both directions exist, read as a general graph, and set
  node1_num later in some way. */
/* The routine compares the maximum node index and #lines, and set #node
  to the larger one. However, if node weight exists, weights will be included 
  in the candidates of maximum index, thus in this case we fix #node := #lines.
  In the case of bipartite graph, the routine compares, but the weights of 
   non-existing lines will be -1. */


/* load edges/arcs (determined by G->flag) from file */
void SGRAPH_load (SGRAPH *G){
  VEC_ID i;
  QUEUE_ID *c, j;
  QUEUE_INT e;
  SETFAMILY *F1, *F2;
  WEIGHT *ww;
  PERM *p;
  QUEUE Q;
  
//  if ( G->flag&LOAD_EDGE ){ F1 = F2 = &G->edge; F1->flag |= LOAD_DBLBUF; }
  if ( G->flag&LOAD_EDGE ) F1 = F2 = &G->edge;
  else {
    F1 = &G->in; F2 = &G->out;
    F1->flag |= LOAD_ARC;
    if ( G->flag & LOAD_TPOSE ) F1->flag |= LOAD_TPOSE;
  }
  F1->flag |= G->flag & (LOAD_ELE + LOAD_EDGEW + LOAD_EDGE + LOAD_RC_SAME + LOAD_ID1 + LOAD_NUM + LOAD_GRAPHNUM);
  if ( !(G->flag&LOAD_BIPARTITE) ) F1->flag |= LOAD_RC_SAME;
  F1->fname = G->fname; F1->wfname = G->wfname;
  SETFAMILY_load (F1);

  if ( G->nwfname ){
#ifdef WEIGHT_DOUBLE
    ARY_LOAD (G->node_w, double, i, G->nwfname, 1, EXIT);
#else
    ARY_LOAD (G->node_w, int, i, G->nwfname, 1, EXIT);
#endif
    reallocx_ (G->node_w, i, SGRAPH_NODE_NUM(*G)+1, 0, EXIT);
  }

  FLOOP (i, 0, F1->t) F1->v[i].v[F1->v[i].t] = F1->t; // set endmark

    // adjast so that #rows and #colums are the same

  if ( !(G->flag&LOAD_EDGE) ){  // make opposite-direction arc
    calloc2 (c, F1->t, EXIT);
    QUEUE_delivery (NULL, c, NULL, F1->v, NULL, F1->t, F1->t);  // comp. size of each adjacent list

    F2->t = F1->clms; F2->clms = F1->t;
    SETFAMILY_alloc (F2, F1->t, c, F1->t, 0);
    if ( F1->w ) SETFAMILY_alloc_weight (F2, c);
    FLOOP (i, 0, F1->t) c[i] = F1->v[i].t;
    FLOOP (i, 0, F1->t){
      if ( F2->rw ) F2->rw[i] = F1->rw[i];
      FLOOP (j, 0, c[i]){
        e = F1->v[i].v[j];
        if ( F2->w ) F2->w[e][F2->v[e].t] = F1->w[i][j];
        QUE_INS (F2->v[e], i);
      }
    }
    free (c);
    F2->clms = F2->t; FLOOP (i, 0, F2->t) F2->v[i].v[F2->v[i].t] = F2->t; // set endmark
  }
  
    // sort the nodes
  F1->flag |= G->flag; F1->rw = G->node_w;
  SETFAMILY_sort (F1);
  F1->rw = NULL; G->perm = F1->rperm; F1->rperm = NULL;
  if ( F1 != F2 ){
    F2->flag |= G->flag; BITRM (F2->flag, LOAD_SIZSORT+LOAD_WSORT);
    SETFAMILY_sort (F2);
    if ( G->flag & (LOAD_SIZSORT+LOAD_WSORT) ){
      ARY_INV_PERM (p, G->perm, F2->t, EXIT);
      if ( F2->w ) ARY_INVPERMUTE (F2->w, p, ww, F2->t, EXIT);
      ARY_INVPERMUTE_ (F2->v, p, Q, F2->t);
      free2 (p);
    }
  }
  print_mes (G, "sgraph: %s ,#nodes %d ,#edges %zd ,#arcs %zd", G->fname, SGRAPH_NODE_NUM(*G), G->edge.eles/2,  G->in.eles);
  if ( G->wfname ) print_mes (G, " ,weight file: %s", G->wfname);
  if ( G->nwfname ) print_mes (G, " ,node weight file: %s", G->nwfname);
  print_mes (G, "\n");
}

/* replace node i by perm[i] */
void SGRAPH_replace_index (SGRAPH *G, PERM *perm, PERM *invperm){
  QUEUE_INT *x;
  VEC_ID i;
  QUEUE Q;
  WEIGHT *w, ww;

  if ( G->edge.v ){
    FLOOP (i, 0, SGRAPH_NODE_NUM(*G))
        MQUE_FLOOP (G->edge.v[i], x) *x = perm[*x];
    ARY_INVPERMUTE (G->edge.v, invperm, Q, SGRAPH_NODE_NUM(*G), EXIT);
  }
  if ( G->in.v ){
    FLOOP (i, 0, SGRAPH_NODE_NUM(*G))
        MQUE_FLOOP (G->in.v[i], x) *x = perm[*x];
    ARY_INVPERMUTE (G->in.v, invperm, Q, SGRAPH_NODE_NUM(*G), EXIT);
  }

  if ( G->out.v ){
    FLOOP (i, 0, SGRAPH_NODE_NUM(*G))
        MQUE_FLOOP (G->out.v[i], x) *x = perm[*x];
    ARY_INVPERMUTE (G->out.v, invperm, Q, SGRAPH_NODE_NUM(*G), EXIT);
  }

  if ( G->edge.w ) ARY_INVPERMUTE (G->edge.w, invperm, w, SGRAPH_NODE_NUM(*G), EXIT);
  if ( G->in.w ) ARY_INVPERMUTE (G->in.w, invperm, w, SGRAPH_NODE_NUM(*G), EXIT);
  if ( G->out.w ) ARY_INVPERMUTE (G->out.w, invperm, w, SGRAPH_NODE_NUM(*G), EXIT);
  if ( G->node_w ) ARY_INVPERMUTE (G->node_w, invperm, ww, SGRAPH_NODE_NUM(*G), EXIT);
  G->perm = perm;

}

/* sort the nodes according to the permutation *tmp */
void SGRAPH_perm_node (SGRAPH *G, PERM *tmp){
  VEC_ID c1=0, c2=G->node1_num, i;
  PERM *perm;
  malloc2 (perm, SGRAPH_NODE_NUM(*G), {free(tmp);EXIT;});
  FLOOP (i, 0, SGRAPH_NODE_NUM(*G))
      if ( tmp[i]<G->node1_num ) perm[tmp[i]] = c1++; else perm[tmp[i]] = c2++;
  ARY_INV_PERM_ (tmp, perm, SGRAPH_NODE_NUM(*G));
  SGRAPH_replace_index (G, perm, tmp);
//  free2 (tmp);
}

/* sort the nodes by Q->t, increasing if flag=1, decreasing if flag=-1 */
void SGRAPH_sort_node (SGRAPH *G, int flag){
  PERM *tmp;
  tmp = qsort_perm_VECt ((VEC *)(G->edge.v), G->edge.t, flag==1?sizeof(QUEUE):-sizeof(QUEUE));
  SGRAPH_perm_node (G, tmp);
  free2 (tmp);
}

/* remove all selfloops */
void SGRAPH_rm_selfloop (SGRAPH *G){
  QUEUE_ID i, j, jj;
  QUEUE_INT x;
  FLOOP (i, 0, SGRAPH_NODE_NUM(*G)){
    if ( G->edge.v ){
      jj = 0;
      FLOOP (j, 0, G->edge.v[i].t){
        x = G->edge.v[i].v[j];
        if ( x != i ){
          if ( j != jj ){
            G->edge.v[i].v[jj] = G->edge.v[i].v[j];
            if ( G->edge.w ) G->edge.w[i][jj] = G->edge.w[i][j];
          }
          jj++;
        }
      }
      G->edge.v[i].t = jj;
    }
    if ( G->in.v ){
      jj = 0;
      FLOOP (j, 0, G->in.v[i].t){
        x = G->in.v[i].v[j];
        if ( x != i ){
          if ( j != jj ){
            G->in.v[i].v[jj] = G->in.v[i].v[j];
            if ( G->in.w ) G->in.w[i][jj] = G->in.w[i][j];
          }
          jj++;
        }
      }
      G->in.v[i].t = jj;
    }
    if ( G->out.v ){
      jj = 0;
      FLOOP (j, 0, G->out.v[i].t){
        x = G->out.v[i].v[j];
        if ( x != i ){
          if ( j != jj ){
            G->out.v[i].v[jj] = G->out.v[i].v[j];
            if ( G->out.w ) G->out.w[i][jj] = G->out.w[i][j];
          }
          jj++;
        }
      }
      G->out.v[i].t = jj;
    }
  }
}
#endif
