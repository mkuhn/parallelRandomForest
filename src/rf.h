/*******************************************************************
   Copyright (C) 2001-7 Leo Breiman, Adele Cutler and Merck & Co., Inc.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*******************************************************************/
#ifndef RF_H
#define RF_H

/* #define RF_DEBUG */

/* test if the bit at position pos is turned on */
#define isBitOn(x,pos) (((x) & (1 << (pos))) > 0)
/* swap two integers */
#define swapInt(a, b) ((a ^= b), (b ^= a), (a ^= b))
/*
void classRF(double *x, int *dimx, int *cl, int *ncl, int *cat, int *maxcat,
	int *sampsize, int *Options, int *ntree, int *nvar,
	int *ipi, double *pi, double *cut, int *nodesize,
        int *outcl, int *counttr, double *prox,
	double *imprt, double *, double *impmat, int *nrnodes, int *ndbigtree,
	int *nodestatus, int *bestvar, int *treemap, int *nodeclass,
	double *xbestsplit, double *pid, double *errtr,
	int *testdat, double *xts, int *clts, int *nts, double *countts,
	int *outclts, int *labelts, double *proxts, double *errts);
*/

void normClassWt(int *cl, const int nsample, const int nclass,
                 const int useWt, double *classwt, int *classFreq);

extern "C" void classForest(int *mdim, int *ntest, int *nclass, int *maxcat,
                 int *nrnodes, int *jbt, double *xts, double *xbestsplit,
                 double *pid, double *cutoff, double *countts, int *treemap,
                 int *nodestatus, int *cat, int *nodeclass, int *jts,
                 int *jet, int *bestvar, int *nodexts, int *ndbigtree,
                 int *keepPred, int *prox, double *proxmatrix, int *nodes);

template <typename T> void regTree(T *x, double *y, int *sampling, int mdim, size_t full_nsample, int nsample,
             int *lDaughter, int *rDaughter, T *upper, double *avnode,
             int *nodestatus, int nrnodes, int *treeSize, int nthsize,
             int mtry, int *mbest, int *cat, double *tgini, int *varUsed);

template <typename T> void detectSplit(T *x, T *xt, double *yl, int ndstart, int ndend, int nodecnt, int nsample,
   T max_x, double sumnode, double critParent,
   double *critvar, double *ubestt);

template <typename T> void findBestSplit(T *x, int *sampling, int *jdex, double *y, int mdim, size_t full_nsample, int nsample,
                   int ndstart, int ndend, int *msplit, double *decsplit,
                   double *ubest, int *ndendl, int *jstat, int mtry,
                   double sumnode, int nodecnt, int *cat);

template <typename T> void predictRegTree(T *x, int nsample, int mdim,
                    int *lDaughter, int *rDaughter, int *nodestatus,
                    double *ypred, T *split, double *nodepred,
                    int *splitVar, int treeSize, int *cat, int maxcat,
                    int *nodex);

void predictClassTree(double *x, int n, int mdim, int *treemap,
                      int *nodestatus, double *xbestsplit,
                      int *bestvar, int *nodeclass,
                      int ndbigtree, int *cat, int nclass,
                      int *jts, int *nodex, int maxcat);

double pack(const int l, const int *icat);
void unpack(const double pack, const int nBits, int *icat);


void zeroInt(int *x, int length);
void zeroDouble(double *x, int length);
void createClass(double *x, int realN, int totalN, int mdim);
void prepare(int *cl, const int nsample, const int nclass, const int ipi,
             double *pi, double *pid, int *nc, double *wtt);
void makeA(double *x, const int mdim, const int nsample, int *cat, int *a,
           int *b);
void modA(int *a, int *nuse, const int nsample, const int mdim, int *cat,
          const int maxcat, int *ncase, int *jin);
void Xtranslate(double *x, int mdim, int nrnodes, int nsample,
                int *bestvar, int *bestsplit, int *bestsplitnext,
                double *xbestsplit, int *nodestatus, int *cat, int treeSize);
void computeProximity(double *prox, int oobprox, int *node, int *inbag,
                      int *oobpair, int n);


// Version for rotated matrix, still used by classification RF
template <typename T> void permuteOOB_legacy(int m, T *x, int *in, int nsample, int mdim) {
    /* Permute the OOB part of a variable in x.
     * Have to keep template function in header
     * Argument:
     *   m: the variable to be permuted
     *   x: the data matrix (variables in rows)
     *   in: vector indicating which case is OOB
     *   nsample: number of cases in the data
     *   mdim: number of variables in the data
     */
    T *tp, tmp;
    int i, last, k, nOOB = 0;

    tp = (T *) Calloc(nsample, T);

    for (i = 0; i < nsample; ++i) {
        /* make a copy of the OOB part of the data into tp (for permuting) */
        if (in[i] == 0) {
            tp[nOOB] = x[m + i*mdim];
            nOOB++;
        }
    }
    /* Permute tp */
    last = nOOB;
    for (i = 0; i < nOOB; ++i) {
        k = (int) last * unif_rand();
        tmp = tp[last - 1];
        tp[last - 1] = tp[k];
        tp[k] = tmp;
        last--;
    }

    /* Copy the permuted OOB data back into x. */
    nOOB = 0;
    for (i = 0; i < nsample; ++i) {
        if (in[i] == 0) {
            x[m + i*mdim] = tp[nOOB];
            nOOB++;
        }
    }
    Free(tp);
}

template <typename T> void permuteOOB(int m, T *x, int *in, int nsample, int mdim) {
    /* Permute the OOB part of a variable in x.
     * Have to keep template function in header
     * Argument:
     *   m: the variable to be permuted
     *   x: the data matrix (variables in rows)
     *   in: vector indicating which case is OOB
     *   nsample: number of cases in the data
     *   mdim: number of variables in the data
     */
    T *tp, tmp;
    int i, last, k, nOOB = 0;

    tp = (T *) Calloc(nsample, T);

    for (i = 0; i < nsample; ++i) {
        /* make a copy of the OOB part of the data into tp (for permuting) */
        if (in[i] == 0) {
            tp[nOOB] = x[m*nsample + i];
            nOOB++;
        }
    }
    /* Permute tp */
    last = nOOB;
    for (i = 0; i < nOOB; ++i) {
        k = (int) last * unif_rand();
        tmp = tp[last - 1];
        tp[last - 1] = tp[k];
        tp[k] = tmp;
        last--;
    }

    /* Copy the permuted OOB data back into x. */
    nOOB = 0;
    for (i = 0; i < nsample; ++i) {
        if (in[i] == 0) {
            x[m*nsample + i] = tp[nOOB];
            nOOB++;
        }
    }
    Free(tp);
}

/* Template of Fortran subroutines to be called from the C wrapper */
extern "C" void F77_NAME(buildtree)(int *a, int *b, int *cl, int *cat,
                                int *maxcat, int *mdim, int *nsample,
                                int *nclass, int *treemap, int *bestvar,
                                int *bestsplit, int *bestsplitnext,
                                double *tgini, int *nodestatus, int *nodepop,
                                int *nodestart, double *classpop,
                                double *tclasspop, double *tclasscat,
                                int *ta, int *nrnodes, int *,
                                int *, int *, int *, int *, int *, int *,
                                double *, double *, double *,
                                int *, int *, int *);


/* maximum number of categories allowed in categorical predictors */
#define MAX_CAT 53

/* Node status */
#define NODE_TERMINAL -1
#define NODE_TOSPLIT  -2
#define NODE_INTERIOR -3

#endif /* RF_H */
