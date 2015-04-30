/*******************************************************************
   Copyright (C) 2001-2012 Leo Breiman, Adele Cutler and Merck & Co., Inc.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*******************************************************************/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "rf.h"


void simpleLinReg(int nsample, double *x, double *y, double *coef, double *mse, int *hasPred);





template <typename T> void regRF(T *x, double *y, int *xdim, int *sampsize,
           int *nthsize, int *nrnodes, int *nTree, int *mtry, int *imp,
           int *cat, int *maxcat, int *jprint, int *doProx, int *oobprox,
           int *biasCorr, double *yptr, double *errimp, double *impmat,
           double *impSD, double *prox, int *treeSize, int *nodestatus,
           int *lDaughter, int *rDaughter, double *avnode, int *mbest,
           T *upper, double *mse, int *keepf, int *replace,
           int *testdat, T *xts, int *nts, double *yts, int *labelts,
           double *yTestPred, double *proxts, double *msets, double *coef,
           int *nout, int *inbag) {
    /*************************************************************************
    Input:
    mdim=number of variables in data set
    nsample=number of cases

    nthsize=number of cases in a node below which the tree will not split,
    setting nthsize=5 generally gives good results.

    nTree=number of trees in run.  200-500 gives pretty good results

    mtry=number of variables to pick to split on at each node.  mdim/3
    seems to give genrally good performance, but it can be
    altered up or down

    imp=1 turns on variable importance.  This is computed for the
    mth variable as the percent rise in the test set mean sum-of-
    squared errors when the mth variable is randomly permuted.

    *************************************************************************/

    double errts = 0.0, averrb, meanY, meanYts, varY, varYts, r, xrand,
           errb = 0.0, resid=0.0, ooberr, ooberrperm, delta, *resOOB;

    double *yb, *ytr, *ytree, *tgini;
    T *xtmp, *xb;

    int k, m, mr, n, nOOB, j, jout, idx, ntest, last, ktmp, nPerm,
        nsample, mdim, keepF, keepInbag;
    int *oobpair, varImp, localImp, *varUsed;

    int *in, *nind, *nodex, *nodexts;

    int *sampling;

    nsample = xdim[0];
    mdim = xdim[1];
    ntest = *nts;
    varImp = imp[0];
    localImp = imp[1];
    nPerm = imp[2];
    keepF = keepf[0];
    keepInbag = keepf[1];

    if (*jprint == 0) *jprint = *nTree + 1;

    yb         = (double *) S_alloc(*sampsize, sizeof(double));
    ytr        = (double *) S_alloc(nsample, sizeof(double));
    xtmp       = (T *) S_alloc(nsample, sizeof(T));
    resOOB     = (double *) S_alloc(nsample, sizeof(double));

    sampling = (int *) S_alloc(*sampsize, sizeof(int));

    in        = (int *) S_alloc(nsample, sizeof(int));
    nodex      = (int *) S_alloc(nsample, sizeof(int));
    varUsed    = (int *) S_alloc(mdim, sizeof(int));
    nind = *replace ? NULL : (int *) S_alloc(nsample, sizeof(int));

    if (*testdat) {
        ytree      = (double *) S_alloc(ntest, sizeof(double));
        nodexts    = (int *) S_alloc(ntest, sizeof(int));
    }
    oobpair = (*doProx && *oobprox) ?
              (int *) S_alloc(nsample * nsample, sizeof(int)) : NULL;

    /* If variable importance is requested, tgini points to the second
       "column" of errimp, otherwise it's just the same as errimp. */
    tgini = varImp ? errimp + mdim : errimp;

    averrb = 0.0;
    meanY = 0.0;
    varY = 0.0;

    zeroDouble(yptr, nsample);
    zeroInt(nout, nsample);
    for (n = 0; n < nsample; ++n) {
        varY += n * (y[n] - meanY)*(y[n] - meanY) / (n + 1);
        meanY = (n * meanY + y[n]) / (n + 1);
    }
    varY /= nsample;

    varYts = 0.0;
    meanYts = 0.0;
    if (*testdat) {
        for (n = 0; n < ntest; ++n) {
            varYts += n * (yts[n] - meanYts)*(yts[n] - meanYts) / (n + 1);
            meanYts = (n * meanYts + yts[n]) / (n + 1);
        }
        varYts /= ntest;
    }

    if (*doProx) {
        zeroDouble(prox, nsample * nsample);
        if (*testdat) zeroDouble(proxts, ntest * (nsample + ntest));
    }

    if (varImp) {
        zeroDouble(errimp, mdim * 2);
        if (localImp) zeroDouble(impmat, nsample * mdim);
    } else {
        zeroDouble(errimp, mdim);
    }
    if (*labelts) zeroDouble(yTestPred, ntest);

    /* print header for running output */
    if (*jprint <= *nTree) {
        Rprintf("     |      Out-of-bag   ");
        if (*testdat) Rprintf("|       Test set    ");
        Rprintf("|\n");
        Rprintf("Tree |      MSE  %%Var(y) ");
        if (*testdat) Rprintf("|      MSE  %%Var(y) ");
        Rprintf("|\n");
    }
    GetRNGstate();
    /*************************************
     * Start the loop over trees.
     *************************************/
    for (j = 0; j < *nTree; ++j) {
        idx = keepF ? j * *nrnodes : 0;
        zeroInt(in, nsample);
        zeroInt(varUsed, mdim);
        /* Draw a random sample for growing a tree. */
        if (*replace) { /* sampling with replacement */
            for (n = 0; n < *sampsize; ++n) {
                xrand = unif_rand();
                k = xrand * nsample;
                in[k] = 1;
                yb[n] = y[k];
                sampling[n] = k;
            }
        } else { /* sampling w/o replacement */
            for (n = 0; n < nsample; ++n) nind[n] = n;
            last = nsample - 1;
            for (n = 0; n < *sampsize; ++n) {
                ktmp = (int) (unif_rand() * (last+1));
                k = nind[ktmp];
                swapInt(nind[ktmp], nind[last]);
                last--;
                in[k] = 1;
                yb[n] = y[k];
                sampling[n] = k;
            }
        }
        if (keepInbag) {
            for (n = 0; n < nsample; ++n) inbag[n + j * nsample] = in[n];
        }
        /* grow the regression tree */
        regTree(x, yb, sampling, mdim, nsample, *sampsize, lDaughter + idx, rDaughter + idx,
                upper + idx, avnode + idx, nodestatus + idx, *nrnodes,
                treeSize + j, *nthsize, *mtry, mbest + idx, cat, tgini,
                varUsed);
        /* predict the OOB data with the current tree */
        /* ytr is the prediction on OOB data by the current tree */
        predictRegTree(x, nsample, mdim, lDaughter + idx,
                       rDaughter + idx, nodestatus + idx, ytr, upper + idx,
                       avnode + idx, mbest + idx, treeSize[j], cat, *maxcat,
                       nodex);
        /* yptr is the aggregated prediction by all trees grown so far */
        errb = 0.0;
        ooberr = 0.0;
        jout = 0; /* jout is the number of cases that has been OOB so far */
        nOOB = 0; /* nOOB is the number of OOB samples for this tree */
        for (n = 0; n < nsample; ++n) {
            if (in[n] == 0) {
                nout[n]++;
                nOOB++;
                yptr[n] = ((nout[n]-1) * yptr[n] + ytr[n]) / nout[n];
                resOOB[n] = ytr[n] - y[n];
                ooberr += resOOB[n] * resOOB[n];
            }
            if (nout[n]) {
                jout++;
                errb += (y[n] - yptr[n]) * (y[n] - yptr[n]);
            }
        }
        errb /= jout;
        /* Do simple linear regression of y on yhat for bias correction. */
        if (*biasCorr) simpleLinReg(nsample, yptr, y, coef, &errb, nout);

        /* predict testset data with the current tree */
        if (*testdat) {
            predictRegTree(xts, ntest, mdim, lDaughter + idx,
                           rDaughter + idx, nodestatus + idx, ytree,
                           upper + idx, avnode + idx,
                           mbest + idx, treeSize[j], cat, *maxcat, nodexts);
            /* ytree is the prediction for test data by the current tree */
            /* yTestPred is the average prediction by all trees grown so far */
            errts = 0.0;
            for (n = 0; n < ntest; ++n) {
                yTestPred[n] = (j * yTestPred[n] + ytree[n]) / (j + 1);
            }
            /* compute testset MSE */
            if (*labelts) {
                for (n = 0; n < ntest; ++n) {
                    resid = *biasCorr ?
                            yts[n] - (coef[0] + coef[1]*yTestPred[n]) :
                            yts[n] - yTestPred[n];
                    errts += resid * resid;
                }
                errts /= ntest;
            }
        }
        /* Print running output. */
        if ((j + 1) % *jprint == 0) {
            Rprintf("%4d |", j + 1);
            Rprintf(" %8.4g %8.2f ", errb, 100 * errb / varY);
            if(*labelts == 1) Rprintf("| %8.4g %8.2f ",
                                          errts, 100.0 * errts / varYts);
            Rprintf("|\n");
        }
        mse[j] = errb;
        if (*labelts) msets[j] = errts;

        /*  DO PROXIMITIES */
        if (*doProx) {
            computeProximity(prox, *oobprox, nodex, in, oobpair, nsample);
            /* proximity for test data */
            if (*testdat) {
                /* In the next call, in and oobpair are not used. */
                computeProximity(proxts, 0, nodexts, in, oobpair, ntest);
                for (n = 0; n < ntest; ++n) {
                    for (k = 0; k < nsample; ++k) {
                        if (nodexts[n] == nodex[k]) {
                            proxts[n + ntest * (k+ntest)] += 1.0;
                        }
                    }
                }
            }
        }

        /* Variable importance */
        if (varImp) {
            for (mr = 0; mr < mdim; ++mr) {
             if (varUsed[mr]) { /* Go ahead if the variable is used */
                 /* make a copy of the m-th variable into xtmp */
                 for (n = 0; n < nsample; ++n)
                     xtmp[n] = x[nsample * mr + n];
                 ooberrperm = 0.0;
                 for (k = 0; k < nPerm; ++k) {
                     permuteOOB(mr, x, in, nsample, mdim);
                     predictRegTree(x, nsample, mdim, lDaughter + idx,
                                    rDaughter + idx, nodestatus + idx, ytr,
                                    upper + idx, avnode + idx, mbest + idx,
                                    treeSize[j], cat, *maxcat, nodex);
                     for (n = 0; n < nsample; ++n) {
                         if (in[n] == 0) {
                             r = ytr[n] - y[n];
                             ooberrperm += r * r;
                             if (localImp) {
                                 impmat[mr + n * mdim] +=
                                     (r*r - resOOB[n]*resOOB[n]) / nPerm;
                             }
                         }
                     }
                 }
                 delta = (ooberrperm / nPerm - ooberr) / nOOB;
                 errimp[mr] += delta;
                 impSD[mr] += delta * delta;
                 /* copy original data back */
                 for (n = 0; n < nsample; ++n)
                     x[nsample * mr + n] = xtmp[n];
             }
         }
        }
    }
    PutRNGstate();
    /* end of tree iterations=======================================*/

    if (*biasCorr) {  /* bias correction for predicted values */
        for (n = 0; n < nsample; ++n) {
            if (nout[n]) yptr[n] = coef[0] + coef[1] * yptr[n];
        }
        if (*testdat) {
            for (n = 0; n < ntest; ++n) {
                yTestPred[n] = coef[0] + coef[1] * yTestPred[n];
            }
        }
    }

    if (*doProx) {
        for (n = 0; n < nsample; ++n) {
            for (k = n + 1; k < nsample; ++k) {
                prox[nsample*k + n] /= *oobprox ?
                    (oobpair[nsample*k + n] > 0 ? oobpair[nsample*k + n] : 1) :
                    *nTree;
                prox[nsample * n + k] = prox[nsample * k + n];
            }
            prox[nsample * n + n] = 1.0;
        }
        if (*testdat) {
            for (n = 0; n < ntest; ++n)
                for (k = 0; k < ntest + nsample; ++k)
                    proxts[ntest*k + n] /= *nTree;
        }
    }

    if (varImp) {
        for (m = 0; m < mdim; ++m) {
            errimp[m] = errimp[m] / *nTree;
            impSD[m] = sqrt( ((impSD[m] / *nTree) -
                              (errimp[m] * errimp[m])) / *nTree );
            if (localImp) {
                for (n = 0; n < nsample; ++n) {
                    impmat[m + n * mdim] /= nout[n];
                }
            }
        }
    }
    for (m = 0; m < mdim; ++m) tgini[m] /= *nTree;
}

extern "C" SEXP callRegRFRaw(SEXP x, SEXP y, SEXP xdim, SEXP sampsize,
               SEXP nthsize, SEXP nrnodes, SEXP nTree, SEXP mtry, SEXP imp,
               SEXP cat, SEXP maxcat, SEXP jprint, SEXP doProx, SEXP oobprox,
               SEXP biasCorr, SEXP yptr, SEXP errimp, SEXP impmat,
               SEXP impSD, SEXP prox, SEXP treeSize, SEXP nodestatus,
               SEXP lDaughter, SEXP rDaughter, SEXP avnode, SEXP mbest,
               SEXP upper, SEXP mse, SEXP keepf, SEXP replace,
               SEXP testdat, SEXP xts, SEXP nts, SEXP yts, SEXP labelts,
               SEXP yTestPred, SEXP proxts, SEXP msets, SEXP coef,
               SEXP nout, SEXP inbag) {


    regRF(RAW(x), REAL(y), INTEGER(xdim), INTEGER(sampsize),
          INTEGER(nthsize), INTEGER(nrnodes), INTEGER(nTree), INTEGER(mtry), INTEGER(imp),
          INTEGER(cat), INTEGER(maxcat), INTEGER(jprint), INTEGER(doProx), INTEGER(oobprox),
          INTEGER(biasCorr), REAL(yptr), REAL(errimp), REAL(impmat),
          REAL(impSD), REAL(prox), INTEGER(treeSize), INTEGER(nodestatus),
          INTEGER(lDaughter), INTEGER(rDaughter), REAL(avnode), INTEGER(mbest),
          RAW(upper), REAL(mse), INTEGER(keepf), INTEGER(replace),
          INTEGER(testdat), RAW(xts), INTEGER(nts), REAL(yts), INTEGER(labelts),
          REAL(yTestPred), REAL(proxts), REAL(msets), REAL(coef),
          INTEGER(nout), INTEGER(inbag)
         );

    return R_NilValue;

}

extern "C" SEXP callRegRFDouble(SEXP x, SEXP y, SEXP xdim, SEXP sampsize,
               SEXP nthsize, SEXP nrnodes, SEXP nTree, SEXP mtry, SEXP imp,
               SEXP cat, SEXP maxcat, SEXP jprint, SEXP doProx, SEXP oobprox,
               SEXP biasCorr, SEXP yptr, SEXP errimp, SEXP impmat,
               SEXP impSD, SEXP prox, SEXP treeSize, SEXP nodestatus,
               SEXP lDaughter, SEXP rDaughter, SEXP avnode, SEXP mbest,
               SEXP upper, SEXP mse, SEXP keepf, SEXP replace,
               SEXP testdat, SEXP xts, SEXP nts, SEXP yts, SEXP labelts,
               SEXP yTestPred, SEXP proxts, SEXP msets, SEXP coef,
               SEXP nout, SEXP inbag) {


    regRF(REAL(x), REAL(y), INTEGER(xdim), INTEGER(sampsize),
          INTEGER(nthsize), INTEGER(nrnodes), INTEGER(nTree), INTEGER(mtry), INTEGER(imp),
          INTEGER(cat), INTEGER(maxcat), INTEGER(jprint), INTEGER(doProx), INTEGER(oobprox),
          INTEGER(biasCorr), REAL(yptr), REAL(errimp), REAL(impmat),
          REAL(impSD), REAL(prox), INTEGER(treeSize), INTEGER(nodestatus),
          INTEGER(lDaughter), INTEGER(rDaughter), REAL(avnode), INTEGER(mbest),
          REAL(upper), REAL(mse), INTEGER(keepf), INTEGER(replace),
          INTEGER(testdat), REAL(xts), INTEGER(nts), REAL(yts), INTEGER(labelts),
          REAL(yTestPred), REAL(proxts), REAL(msets), REAL(coef),
          INTEGER(nout), INTEGER(inbag)
         );

    return R_NilValue;

}


/*----------------------------------------------------------------------*/
template <typename T> void regForest(T *x, double *ypred, int *mdim, int *n,
               int *ntree, int *lDaughter, int *rDaughter,
               int *nodestatus, int *nrnodes, T *xsplit,
               double *avnodes, int *mbest, int *treeSize, int *cat,
               int *maxcat, int *keepPred, double *allpred, int *doProx,
               double *proxMat, int *nodes, int *nodex) {
    int i, j, *junk;
    size_t idx1, idx2;
    double *ytree;

    junk = NULL;
    ytree = (double *) S_alloc(*n, sizeof(double));
    if (*nodes) {
        zeroInt(nodex, *n * *ntree);
    } else {
        zeroInt(nodex, *n);
    }
    if (*doProx) zeroDouble(proxMat, *n * *n);
    if (*keepPred) zeroDouble(allpred, *n * *ntree);
    idx1 = 0;
    idx2 = 0;
    for (i = 0; i < *ntree; ++i) {

        zeroDouble(ytree, *n);
        predictRegTree(x, *n, *mdim, lDaughter + idx1, rDaughter + idx1,
                       nodestatus + idx1, ytree, xsplit + idx1,
                       avnodes + idx1, mbest + idx1, treeSize[i], cat, *maxcat,
                       nodex + idx2);

        for (j = 0; j < *n; ++j) ypred[j] += ytree[j];
        if (*keepPred) {
            for (j = 0; j < *n; ++j) allpred[j + i * *n] = ytree[j];
        }

        /* if desired, do proximities for this round */
        if (*doProx) computeProximity(proxMat, 0, nodex + idx2, junk,
                                          junk, *n);
        idx1 += *nrnodes; /* increment the offset */
        if (*nodes) idx2 += *n;
    }
    for (i = 0; i < *n; ++i) ypred[i] /= *ntree;
    if (*doProx) {
        for (i = 0; i < *n; ++i) {
            for (j = i + 1; j < *n; ++j) {
                proxMat[i + j * *n] /= *ntree;
                proxMat[j + i * *n] = proxMat[i + j * *n];
            }
            proxMat[i + i * *n] = 1.0;
        }
    }
}

extern "C" SEXP callRegForestRaw(SEXP x, SEXP ypred, SEXP mdim, SEXP n,
               SEXP ntree, SEXP lDaughter, SEXP rDaughter,
               SEXP nodestatus, SEXP nrnodes, SEXP xsplit,
               SEXP avnodes, SEXP mbest, SEXP treeSize, SEXP cat,
               SEXP maxcat, SEXP keepPred, SEXP allpred, SEXP doProx,
               SEXP proxMat, SEXP nodes, SEXP nodex)
{

    regForest(RAW(x), REAL(ypred), INTEGER(mdim), INTEGER(n),
               INTEGER(ntree), INTEGER(lDaughter), INTEGER(rDaughter),
               INTEGER(nodestatus), INTEGER(nrnodes), RAW(xsplit),
               REAL(avnodes), INTEGER(mbest), INTEGER(treeSize), INTEGER(cat),
               INTEGER(maxcat), INTEGER(keepPred), REAL(allpred), INTEGER(doProx),
               REAL(proxMat), INTEGER(nodes), INTEGER(nodex));

    return R_NilValue;
}

extern "C" SEXP callRegForestDouble(SEXP x, SEXP ypred, SEXP mdim, SEXP n,
               SEXP ntree, SEXP lDaughter, SEXP rDaughter,
               SEXP nodestatus, SEXP nrnodes, SEXP xsplit,
               SEXP avnodes, SEXP mbest, SEXP treeSize, SEXP cat,
               SEXP maxcat, SEXP keepPred, SEXP allpred, SEXP doProx,
               SEXP proxMat, SEXP nodes, SEXP nodex)
{

    regForest(REAL(x), REAL(ypred), INTEGER(mdim), INTEGER(n),
               INTEGER(ntree), INTEGER(lDaughter), INTEGER(rDaughter),
               INTEGER(nodestatus), INTEGER(nrnodes), REAL(xsplit),
               REAL(avnodes), INTEGER(mbest), INTEGER(treeSize), INTEGER(cat),
               INTEGER(maxcat), INTEGER(keepPred), REAL(allpred), INTEGER(doProx),
               REAL(proxMat), INTEGER(nodes), INTEGER(nodex));

    return R_NilValue;
}


void simpleLinReg(int nsample, double *x, double *y, double *coef,
                  double *mse, int *hasPred) {
    /* Compute simple linear regression of y on x, returning the coefficients,
       the average squared residual, and the predicted values (overwriting y). */
    int i, nout = 0;
    double sxx=0.0, sxy=0.0, xbar=0.0, ybar=0.0;
    double dx = 0.0, dy = 0.0, py=0.0;

    for (i = 0; i < nsample; ++i) {
        if (hasPred[i]) {
            nout++;
            xbar += x[i];
            ybar += y[i];
        }
    }
    xbar /= nout;
    ybar /= nout;

    for (i = 0; i < nsample; ++i) {
        if (hasPred[i]) {
            dx = x[i] - xbar;
            dy = y[i] - ybar;
            sxx += dx * dx;
            sxy += dx * dy;
        }
    }
    coef[1] = sxy / sxx;
    coef[0] = ybar - coef[1] * xbar;

    *mse = 0.0;
    for (i = 0; i < nsample; ++i) {
        if (hasPred[i]) {
            py = coef[0] + coef[1] * x[i];
            dy = y[i] - py;
            *mse += dy * dy;
            /* y[i] = py; */
        }
    }
    *mse /= nout;
    return;
}



template <typename T> void regTree(T *x, double *y, int *sampling, int mdim, size_t full_nsample, int nsample, int *lDaughter,
             int *rDaughter,
             T *upper, double *avnode, int *nodestatus, int nrnodes,
             int *treeSize, int nthsize, int mtry, int *mbest, int *cat,
             double *tgini, int *varUsed) {
    int i, j, k, m, ncur, *jdex, *nodestart, *nodepop;
    int ndstart, ndend, ndendl, nodecnt, jstat, msplit;
    double d, av, decsplit, ubest, sumnode;

    nodestart = (int *) Calloc(nrnodes, int);
    nodepop   = (int *) Calloc(nrnodes, int);

    /* initialize some arrays for the tree */
    zeroInt(nodestatus, nrnodes);
    zeroInt(nodestart, nrnodes);
    zeroInt(nodepop, nrnodes);
    zeroDouble(avnode, nrnodes);

    jdex = (int *) Calloc(nsample, int);
    for (i = 1; i <= nsample; ++i) jdex[i-1] = i;

    ncur = 0;
    nodestart[0] = 0;
    nodepop[0] = nsample;
    nodestatus[0] = NODE_TOSPLIT;

    /* compute mean and sum of squares for Y */
    av = 0.0;
    for (i = 0; i < nsample; ++i) {
        d = y[jdex[i] - 1];
        av = (i * av + d) / (i + 1);
    }
    avnode[0] = av;

    /* start main loop */
    for (k = 0; k < nrnodes - 2; ++k) {
        if (k > ncur || ncur >= nrnodes - 2) break;
        /* skip if the node is not to be split */
        if (nodestatus[k] != NODE_TOSPLIT) continue;

        /* initialize for next call to findbestsplit */
        ndstart = nodestart[k];
        ndend = ndstart + nodepop[k] - 1;
        nodecnt = nodepop[k];
        sumnode = nodecnt * avnode[k];
        jstat = 0;
        decsplit = 0.0;

        findBestSplit(x, sampling, jdex, y, mdim, full_nsample, nsample, ndstart, ndend, &msplit,
                      &decsplit, &ubest, &ndendl, &jstat, mtry, sumnode,
                      nodecnt, cat);
        if (jstat == 1) {
            /* Node is terminal: Mark it as such and move on to the next. */
            nodestatus[k] = NODE_TERMINAL;
            continue;
        }
        /* Found the best split. */
        mbest[k] = msplit;
        varUsed[msplit - 1] = 1;
        upper[k] = ubest;
        tgini[msplit - 1] += decsplit;
        nodestatus[k] = NODE_INTERIOR;

        /* leftnode no.= ncur+1, rightnode no. = ncur+2. */
        nodepop[ncur + 1] = ndendl - ndstart + 1;
        nodepop[ncur + 2] = ndend - ndendl;
        nodestart[ncur + 1] = ndstart;
        nodestart[ncur + 2] = ndendl + 1;

        /* compute mean and sum of squares for the left daughter node */
        av = 0.0;
        for (j = ndstart; j <= ndendl; ++j) {
            d = y[jdex[j]-1];
            m = j - ndstart;
            av = (m * av + d) / (m+1);
        }
        avnode[ncur+1] = av;
        nodestatus[ncur+1] = NODE_TOSPLIT;
        if (nodepop[ncur + 1] <= nthsize) {
            nodestatus[ncur + 1] = NODE_TERMINAL;
        }

        /* compute mean and sum of squares for the right daughter node */
        av = 0.0;
        for (j = ndendl + 1; j <= ndend; ++j) {
            d = y[jdex[j]-1];
            m = j - (ndendl + 1);
            av = (m * av + d) / (m + 1);
        }
        avnode[ncur + 2] = av;
        nodestatus[ncur + 2] = NODE_TOSPLIT;
        if (nodepop[ncur + 2] <= nthsize) {
            nodestatus[ncur + 2] = NODE_TERMINAL;
        }

        /* map the daughter nodes */
        lDaughter[k] = ncur + 1 + 1;
        rDaughter[k] = ncur + 2 + 1;
        /* Augment the tree by two nodes. */
        ncur += 2;
    }
    *treeSize = nrnodes;
    for (k = nrnodes - 1; k >= 0; --k) {
        if (nodestatus[k] == 0) (*treeSize)--;
        if (nodestatus[k] == NODE_TOSPLIT) {
            nodestatus[k] = NODE_TERMINAL;
        }
    }
    Free(nodestart);
    Free(jdex);
    Free(nodepop);
}

/*--------------------------------------------------------------
Rely on the fact that we have only relatively few values (e.g. 0/1/2),
so don't sort the values but rather scan through each split
--------------------------------------------------------------*/
template <> void detectSplit(unsigned char *x, unsigned char *xt, double *yl, int ndstart, int ndend, int nodecnt, int nsample,
   unsigned char max_x, double sumnode, double critParent,
   double *critvar, double *ubestt)
{
    if (max_x == 0) return;

    double suml = 0.0;
    int npopl = 0;

    /* gather all cases where x is 0 */
    for (size_t j = ndstart; j <= ndend; ++j) {
        unsigned char flag = xt[j] == 0;
        suml += flag * yl[j];
        npopl += flag;
    }

    unsigned char last_split = 0;

    /* try all possible splits */
    for (unsigned char split = 1; split <= max_x; split++) {

        /* compute effect of split here, but don't split now as we don't know the next potential value */
        int npopr = nodecnt - npopl;
        double sumr = sumnode - suml;
        double crit = (suml * suml / npopl) + (sumr * sumr / npopr) - critParent;

        /* we don't need to compute anything for the last split, as this would only set "sumr" to 0 */
        if (split < max_x) {
            int last_npopl = npopl;

            for (size_t j = ndstart; j <= ndend; ++j) {
                unsigned char flag = xt[j] == split;
                suml += flag * yl[j];
                npopl += flag;
            }

            if (last_npopl == npopl) continue;
        }

        if (crit > *critvar) {
            /* be careful not to cause an overflow as split and last_split are chars */
            *ubestt = split / 2.0 + last_split / 2.0;
            *critvar = crit;
        }

        last_split = split;
    }
}

/*--------------------------------------------------------------
Use QuickSort to find the best split for continuous data
--------------------------------------------------------------*/
template <> void detectSplit(double *x, double *xt, double *yl, int ndstart, int ndend, int nodecnt, int nsample,
   double max_x, double sumnode, double critParent,
   double *critvar, double *ubestt)
{
    double *v = (double *) Calloc(nsample, double);
    int *ncase = (int *) Calloc(nsample, int);

    /* copy the x data in this node. */

    for (size_t j = ndstart; j <= ndend; ++j) v[j] = xt[j];
    for (size_t j = 1; j <= nsample; ++j) ncase[j - 1] = j;
    R_qsort_I(v, ncase, ndstart + 1, ndend + 1);
    if (v[ndstart] >= v[ndend]) return;

    double suml = 0;
    int npopl = 0;

    /* Search through the "gaps" in the x-variable. */
    for (size_t j = ndstart; j <= ndend - 1; ++j) {
        suml += yl[ncase[j] - 1];
        npopl++;
        if (v[j] < v[j+1]) {
            int npopr = nodecnt - npopl;
            double sumr = sumnode - suml;
            double crit = (suml * suml / npopl) + (sumr * sumr / npopr) - critParent;
            if (crit > *critvar) {
                *ubestt = (v[j] + v[j+1]) / 2.0;
                *critvar = crit;
            }
        }
    }

    Free(v);
    Free(ncase);
}

/*--------------------------------------------------------------
findBestSplit for discrete and continuous variables
--------------------------------------------------------------*/
template <typename T> void findBestSplit(T *x, int *sampling, int *jdex, double *y, int mdim, size_t full_nsample, int nsample,
                   int ndstart, int ndend, int *msplit, double *decsplit,
                   double *ubest, int *ndendl, int *jstat, int mtry,
                   double sumnode, int nodecnt, int *cat) {
    int last, ncat[MAX_CAT], icat[MAX_CAT], lc, nl, nr;
    size_t i, j, kv, l;
    int *mind, *ncase;
    T *xt;
    T max_x;
    double *ut, *yl, sumcat[MAX_CAT], avcat[MAX_CAT], tavcat[MAX_CAT], ubestt;
    double critmax, critvar, critParent;

    ut = (double *) Calloc(nsample, double);
    xt = (T *) Calloc(nsample, T);
    yl = (double *) Calloc(nsample, double);
    mind  = (int *) Calloc(mdim, int);
    zeroDouble(avcat, MAX_CAT);
    zeroDouble(tavcat, MAX_CAT);

    /* START BIG LOOP */
    *msplit = -1;
    *decsplit = 0.0;
    critmax = 0.0;
    ubestt = 0.0;
    for (i=0; i < mdim; ++i) mind[i] = i;

    last = mdim - 1;
    for (i = 0; i < mtry; ++i) {
        critvar = 0.0;
        j = (int) (unif_rand() * (last+1));
        kv = mind[j];
        swapInt(mind[j], mind[last]);
        last--;

        lc = cat[kv];

        max_x = 0;

        if (lc == 1) {
            /* numeric variable */
            for (j = ndstart; j <= ndend; ++j) {
                T current_x = x[full_nsample * kv + sampling[jdex[j] - 1]];
                xt[j] = current_x;
                if (current_x > max_x) max_x = current_x;
                yl[j] = y[jdex[j] - 1];
            }
        } else {
            /* categorical variable */
            zeroInt(ncat, MAX_CAT);
            zeroDouble(sumcat, MAX_CAT);
            for (j = ndstart; j <= ndend; ++j) {
                l = (int) x[full_nsample * kv + sampling[jdex[j] - 1]];
                sumcat[l - 1] += y[jdex[j] - 1];
                ncat[l - 1] ++;
            }
            /* Compute means of Y by category. */
            for (j = 0; j < lc; ++j) {
                avcat[j] = ncat[j] ? sumcat[j] / ncat[j] : 0.0;
            }
            /* Make the category mean the `pseudo' X data. */
            for (j = 0; j < nsample; ++j) {
                T current_x = avcat[(int) x[full_nsample * kv + sampling[jdex[j] - 1]] - 1];
                xt[j] = current_x;
                if (current_x > max_x) max_x = current_x;
                yl[j] = y[jdex[j] - 1];
            }
        }

        critParent = sumnode * sumnode / nodecnt;

        detectSplit(x, xt, yl, ndstart, ndend, nodecnt, nsample, max_x, sumnode, critParent, &critvar, &ubestt);

        if (critvar > critmax) {
            *ubest = ubestt;
            *msplit = kv + 1;
            critmax = critvar;
            for (j = ndstart; j <= ndend; ++j) {
                ut[j] = xt[j];
            }
            if (cat[kv] > 1) {
                for (j = 0; j < cat[kv]; ++j) tavcat[j] = avcat[j];
            }
        }
    }
    *decsplit = critmax;

    /* If best split can not be found, set to terminal node and return. */
    if (*msplit != -1) {
        ncase = (int *) Calloc(nsample, int);
        nl = ndstart;
        for (j = ndstart; j <= ndend; ++j) {
            if (ut[j] <= *ubest) {
                nl++;
                ncase[nl-1] = jdex[j];
            }
        }
        *ndendl = imax2(nl - 1, ndstart);
        nr = *ndendl + 1;
        for (j = ndstart; j <= ndend; ++j) {
            if (ut[j] > *ubest) {
                if (nr >= nsample) break;
                nr++;
                ncase[nr - 1] = jdex[j];
            }
        }
        if (*ndendl >= ndend) *ndendl = ndend - 1;
        for (j = ndstart; j <= ndend; ++j) jdex[j] = ncase[j];

        lc = cat[*msplit - 1];
        if (lc > 1) {
            for (j = 0; j < lc; ++j) {
                icat[j] = (tavcat[j] < *ubest) ? 1 : 0;
            }
            *ubest = pack(lc, icat);
        }
        Free(ncase);
    } else *jstat = 1;

    Free(mind);
    Free(yl);
    Free(xt);
    Free(ut);
}


/*====================================================================*/
template <typename T> void predictRegTree(T *x, int nsample, int mdim,
                    int *lDaughter, int *rDaughter, int *nodestatus,
                    double *ypred, T *split, double *nodepred,
                    int *splitVar, int treeSize, int *cat, int maxcat,
                    int *nodex) {
    size_t i, j, k, m;
    int *cbestsplit;
    unsigned int npack;

    /* decode the categorical splits */
    if (maxcat > 1) {
        cbestsplit = (int *) Calloc(maxcat * treeSize, int);
        zeroInt(cbestsplit, maxcat * treeSize);
        for (i = 0; i < treeSize; ++i) {
            if (nodestatus[i] != NODE_TERMINAL && cat[splitVar[i] - 1] > 1) {
                npack = (unsigned int) split[i];
                /* unpack `npack' into bits */
                for (j = 0; npack; npack >>= 1, ++j) {
                    cbestsplit[j + i*maxcat] = npack & 1;
                }
            }
        }
    }

    for (i = 0; i < nsample; ++i) {
        k = 0;
        while (nodestatus[k] != NODE_TERMINAL) { /* go down the tree */
            m = splitVar[k] - 1;
            if (cat[m] == 1) {
                k = (x[nsample*m + i] <= split[k]) ?
                    lDaughter[k] - 1 : rDaughter[k] - 1;
            } else {
                /* Split by a categorical predictor */
                k = cbestsplit[(int) x[nsample*m + i] - 1 + k * maxcat] ?
                    lDaughter[k] - 1 : rDaughter[k] - 1;
            }
        }
        /* terminal node: assign prediction and move on to next */
        ypred[i] = nodepred[k];
        nodex[i] = k + 1;
    }
    if (maxcat > 1) Free(cbestsplit);
}



