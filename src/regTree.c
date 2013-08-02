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

/******************************************************************
 * buildtree and findbestsplit routines translated from Leo's
 * original Fortran code.
 *
 *      copyright 1999 by leo Breiman
 *      this is free software and can be used for any purpose.
 *      It comes with no guarantee.
 *
 ******************************************************************/
#include <Rmath.h>
#include <R.h>
#include "rf.h"

void regTree(unsigned char *x, double *y, int *sampling, int mdim, size_t full_nsample, int nsample, int *lDaughter,
             int *rDaughter,
             unsigned char *upper, double *avnode, int *nodestatus, int nrnodes,
             int *treeSize, int nthsize, int mtry, int *mbest, int *cat,
             double *tgini, int *varUsed) {
    int i, j, k, m, ncur, *jdex, *nodestart, *nodepop;
    int ndstart, ndend, ndendl, nodecnt, jstat, msplit;
    double d, ss, av, decsplit, ubest, sumnode;

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
    ss = 0.0;
    for (i = 0; i < nsample; ++i) {
        d = y[jdex[i] - 1];
        ss += i * (av - d) * (av - d) / (i + 1);
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
        ss = 0.0;
        for (j = ndstart; j <= ndendl; ++j) {
            d = y[jdex[j]-1];
            m = j - ndstart;
            ss += m * (av - d) * (av - d) / (m + 1);
            av = (m * av + d) / (m+1);
        }
        avnode[ncur+1] = av;
        nodestatus[ncur+1] = NODE_TOSPLIT;
        if (nodepop[ncur + 1] <= nthsize) {
            nodestatus[ncur + 1] = NODE_TERMINAL;
        }

        /* compute mean and sum of squares for the right daughter node */
        av = 0.0;
        ss = 0.0;
        for (j = ndendl + 1; j <= ndend; ++j) {
            d = y[jdex[j]-1];
            m = j - (ndendl + 1);
            ss += m * (av - d) * (av - d) / (m + 1);
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

/*--------------------------------------------------------------*/
void findBestSplit(unsigned char *x, int *sampling, int *jdex, double *y, int mdim, size_t full_nsample, int nsample,
                   int ndstart, int ndend, int *msplit, double *decsplit,
                   double *ubest, int *ndendl, int *jstat, int mtry,
                   double sumnode, int nodecnt, int *cat) {
    int last, ncat[32], icat[32], lc, nl, nr, npopl, npopr, last_npopl;
    size_t i, j, kv, l;
    int *split_nodes;
    double *split_sum;
    int *mind, *ncase;
    unsigned char max_x, split, current_x, last_split;
    double *ut, sumcat[32], avcat[32], tavcat[32], ubestt;
    double crit, critmax, critvar, suml, sumr, d, critParent, current_y;

    ut = (double *) Calloc(nsample, double);
    mind  = (int *) Calloc(mdim, int);

    split_sum = (double *) Calloc(256, double);
    split_nodes = (int *) Calloc(256, int);

    zeroDouble(avcat, 32);
    zeroDouble(tavcat, 32);

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
        max_x = 0;

        zeroDouble(split_sum, 256);
        zeroInt(split_nodes, 256);

        lc = cat[kv];
        if (lc == 1) {
            /* numeric variable */
            for (j = ndstart; j <= ndend; ++j) {
                current_x = x[full_nsample * kv + sampling[jdex[j] - 1]];
                split_sum[current_x] += y[jdex[j] - 1];
                ++split_nodes[current_x];
                if (current_x > max_x) max_x = current_x;
            }
        } else {
            /* categorical variable */
            zeroInt(ncat, 32);
            zeroDouble(sumcat, 32);
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
                current_x = avcat[(int) x[full_nsample * kv + sampling[jdex[j] - 1]] - 1];
                split_sum[current_x] += y[jdex[j] - 1];
                ++split_nodes[current_x];
                if (current_x > max_x) max_x = current_x;
            }
        }

        if (max_x == 0) continue;

        /* changed implementation: rely on the fact that we have only relatively few values (e.g. 0/1/2),
           so don't sort the values but rather scan through each split
         */
        critParent = sumnode * sumnode / nodecnt;
        suml = split_sum[0];
        npopl = split_nodes[0];
        crit = 0.0;

        last_split = 0;

        /* try all possible splits */
        for (split = 1; split <= max_x; split++) {

            /* compute effect of split here, but don't split now as we don't know the next potential value */
            sumr = sumnode - suml;
            npopr = nodecnt - npopl;
            crit = (suml * suml / npopl) + (sumr * sumr / npopr) - critParent;

            /* we don't need to compute anything for the last split, as this would only set "sumr" to 0 */
            if (split < max_x) {
                last_npopl = npopl;

                npopl += split_nodes[split];
                suml += split_sum[split];

                if (last_npopl == npopl) continue;
            }

            if (crit > critvar) {
                /* be careful not to cause an overflow as split and last_split are chars */
                ubestt = split / 2.0 + last_split / 2.0;
                critvar = crit;
            }

            last_split = split;
        }

        if (critvar > critmax) {
            *ubest = ubestt;
            *msplit = kv + 1;
            critmax = critvar;
            if (lc == 1) {
                for (j = ndstart; j <= ndend; ++j) {
                    ut[j] = x[full_nsample * kv + sampling[jdex[j] - 1]];
                }
            } else {
                for (j = ndstart; j <= ndend; ++j) {
                    ut[j] = avcat[(int) x[full_nsample * kv + sampling[jdex[j] - 1]] - 1];
                }
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

    Free(split_nodes);
    Free(split_sum);
    Free(mind);
    Free(ut);
}

/*====================================================================*/
void predictRegTree(unsigned char *x, int nsample, int mdim,
                    int *lDaughter, int *rDaughter, int *nodestatus,
                    double *ypred, unsigned char *split, double *nodepred,
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
