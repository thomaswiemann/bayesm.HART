/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#ifndef GUARD_info_h
#define GUARD_info_h
#include "common.h"
//data
class dinfo {
public:
   dinfo() {p=0;n=0;x=0;y=0;}
   size_t p;  //number of vars
   size_t n;  //number of observations
   double *x; // jth var of ith obs is *(x + p*i+j)
   double *y; // ith y is *(y+i) or y[i]
};
//prior and mcmc
//
// nmin / ess_min control the minimum-leaf-size guard used by birth/death MH
// (see heterbd.cpp).  Defaults preserve original behavior:
//
//   nmin    = 5   - hard lower bound on raw observation counts per leaf,
//                   matching the historical hard-coded threshold in heterbd.
//   ess_min = 0.0 - no effective-sample-size guard; turn on (e.g. ess_min=5)
//                   for varying-coefficient (phi-tree) heterbart updates where
//                   raw counts no longer measure information about the leaf
//                   parameter.  See discussions/2026-05-12-phibart-vs-heterbart.md.
class pinfo
{
public:
   pinfo(): pbd(1.0),pb(.5),alpha(.95),mybeta(2.0),tau(1.0),
            nmin(5),ess_min(0.0) {}
//mcmc info
   double pbd; //prob of birth/death
   double pb;  //prob of birth
//prior info
   double alpha;
   double mybeta;
   double tau;
//birth/death leaf-size guards
   size_t nmin;     // raw min observations per leaf
   double ess_min;  // min effective sample size = sum_{i in leaf} w_i
   void pr() {
      cout << "pbd,pb: " << pbd << ", " << pb << std::endl;
      cout << "alpha,beta,tau: " << alpha <<
             ", " << mybeta << ", " << tau << std::endl;
      cout << "nmin,ess_min: " << nmin << ", " << ess_min << std::endl;
   }
};

#endif
