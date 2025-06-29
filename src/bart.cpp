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

#include "bart.h"
#include "logging_bart.h"

//--------------------------------------------------
//constructor
bart::bart():m(200),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false) {}
bart::bart(size_t im):m(im),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false) {}
bart::bart(const bart& ib):m(ib.m),t(m),pi(ib.pi),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false)
{
   this->t = ib.t;
}
bart::~bart()
{
   if(allfit) delete[] allfit;
   if(r) delete[] r;
   if(ftemp) delete[] ftemp;
}

//--------------------------------------------------
//operators
bart& bart::operator=(const bart& rhs)
{
   if(&rhs != this) {

      this->t = rhs.t;
      this->m = t.size();

      this->pi = rhs.pi;

      p=0;n=0;x=0;y=0;
      xi.clear();

      if(allfit) {delete[] allfit; allfit=0;}
      if(r) {delete[] r; r=0;}
      if(ftemp) {delete[] ftemp; ftemp=0;}

   }
   return *this;
}
//--------------------------------------------------
//get,set
void bart::setm(size_t m)
{
   t.resize(m);
   this->m = t.size();

   if(allfit && (xi.size()==p)) predict(p,n,x,allfit);
}

//--------------------------------------------------
void bart::setxinfo(xinfo& _xi)
{
   size_t p=_xi.size();
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
     size_t nc=_xi[i].size();
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = _xi[i][j];
   }
}
//--------------------------------------------------
void bart::setdata(size_t p, size_t n, double *x, double *y, size_t numcut)
{
  int* nc = new int[p];
  for(size_t i=0; i<p; ++i) nc[i]=numcut;
  this->setdata(p, n, x, y, nc);
  delete [] nc;
}

void bart::setdata(size_t p, size_t n, double *x, double *y, int *nc)
{
   this->p=p; this->n=n; this->x=x; this->y=y;
   if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);

   if(allfit) delete[] allfit;
   allfit = new double[n];
   predict(p,n,x,allfit);

   if(r) delete[] r;
   r = new double[n];

   if(ftemp) delete[] ftemp;
   ftemp = new double[n];

   di.n=n; di.p=p; di.x = &x[0]; di.y=r;
   for(size_t j=0;j<p;j++){
     nv.push_back(0);
     pv.push_back(1/(double)p);
   }
}
//--------------------------------------------------
void bart::predict(size_t p, size_t n, double *x, double *fp)
//uses: m,t,xi
{
    double* fptemp = new double[n];
    //printf("predict - m: %zu\n", m);  // Print number of trees

    for (size_t j = 0; j < n; j++) fp[j] = 0.0;
    for (size_t j = 0; j < m; j++) {
        fit(t[j], xi, p, n, x, fptemp);
        // For first tree, print range and sample values
        if (j == 0) {
            double min_val = fptemp[0], max_val = fptemp[0];
            for (size_t k = 1; k < n; k++) {
                if (fptemp[k] < min_val) min_val = fptemp[k];
                if (fptemp[k] > max_val) max_val = fptemp[k];
            }
            //printf("predict - First tree: min=%f, max=%f, first 3 values: %f %f %f\n",
                //min_val, max_val, fptemp[0], fptemp[1], fptemp[2]);
        }
        for (size_t k = 0; k < n; k++) fp[k] += fptemp[k];
    }

    // Also print range for final predictions
    double min_fp = fp[0], max_fp = fp[0];
    for (size_t k = 1; k < n; k++) {
        if (fp[k] < min_fp) min_fp = fp[k];
        if (fp[k] > max_fp) max_fp = fp[k];
    }
    //printf("predict - Final: min=%f, max=%f, first 3 values: %f %f %f\n",
        //min_fp, max_fp, fp[0], fp[1], fp[2]);

    delete[] fptemp;
}
//--------------------------------------------------


void bart::draw(double sigma, rn& gen) {
	std::vector<double> tree_preds;

	for (size_t j = 0; j < m; j++) {
		//if (debug && debug_obs >= 0) {
		//	tree::npv nodes;
		//	t[j].getnodes(nodes);
		//	bart_logging::log_tree_structure(*this, j, "Tree structure before update", nodes);
		//}

		fit(t[j], xi, p, n, x, ftemp);

		//if (debug && debug_obs >= 0) {
		//	double* xx = x + debug_obs * p;
		//	tree::tree_p node = t[j].bn(xx, xi);
		//	bart_logging::log_draw_state(*this, j, "Before removing", debug_obs, allfit[debug_obs]);
		//	bart_logging::log_draw_state(*this, j, "Current ftemp", debug_obs, ftemp[debug_obs]);
		//	bart_logging::log_node_details(*this, j, node->nid(), node->gettheta());
		//	tree_preds.push_back(ftemp[debug_obs]);
		//}

		for (size_t k = 0; k < n; k++) {
			allfit[k] = allfit[k] - ftemp[k];
			r[k] = y[k] - allfit[k];
		}

		//if (debug && debug_obs >= 0) {
		//	bart_logging::log_draw_state(*this, j, "After removing tree", debug_obs, allfit[debug_obs]);
		//}

		bd(t[j], xi, di, pi, sigma, nv, pv, aug, gen);
		drmu(t[j], xi, di, pi, sigma, gen);

		//if (debug && debug_obs >= 0) {
		//	tree::npv nodes;
		//	t[j].getnodes(nodes);
		//	bart_logging::log_tree_structure(*this, j, "Tree structure after bd/drmu", nodes);
		//}

		// Before fit
		//if (debug && debug_obs >= 0) {
		//	double* xx = x + debug_obs * p;
		//	tree::tree_p node = t[j].bn(xx, xi);
		//	bart_logging::log_tree_traversal(*this, j, debug_obs, xx, xi, node);
		//}

		fit(t[j], xi, p, n, x, ftemp);

		// After fit  
		//if (debug && debug_obs >= 0) {
		//	double* xx = x + debug_obs * p;
		//	tree::tree_p node = t[j].bn(xx, xi);
		//	bart_logging::log_tree_traversal(*this, j, debug_obs, xx, xi, node);

		//	if (fabs(node->gettheta() - ftemp[debug_obs]) > 1e-10) {
		//		bart_logging::log_prediction_mismatch(*this, j, debug_obs,
		//			node->gettheta(), ftemp[debug_obs], xi);
		//	}
		//}

		//if (debug && debug_obs >= 0) {
		//	double* xx = x + debug_obs * p;
		//	tree::tree_p node = t[j].bn(xx, xi);
		//	bart_logging::log_draw_state(*this, j, "Before adding updated tree", debug_obs, ftemp[debug_obs]);
		//	bart_logging::log_tree_prediction(*this, j, debug_obs, node->gettheta(), node->nid(), ftemp[debug_obs]);
		//	bart_logging::log_node_details(*this, j, node->nid(), node->gettheta());
		//}

		for (size_t k = 0; k < n; k++) allfit[k] += ftemp[k];

		//if (debug && debug_obs >= 0) {
		//	bart_logging::log_draw_state(*this, j, "After adding updated tree", debug_obs, allfit[debug_obs]);
		//	tree::npv nodes;
		//	t[j].getnodes(nodes);
		//	bart_logging::log_tree_structure(*this, j, "Final tree structure", nodes);
		//}
	}

	//if (debug && debug_obs >= 0) {
	//	bart_logging::log_prediction_details(*this, debug_obs, tree_preds, allfit[debug_obs]);
	//}

	if (dartOn) {
		draw_s(nv, lpv, theta, gen);
		draw_theta0(const_theta, theta, lpv, a, b, rho, gen);
		for (size_t j = 0; j < p; j++) pv[j] = ::exp(lpv[j]);
	}
}

//void bart::draw(double sigma, rn& gen) {
//	std::vector<double> tree_preds;  // Store predictions for logging
//
//	for (size_t j = 0; j < m; j++) {
//		// Remove current tree's contribution
//		fit(t[j], xi, p, n, x, ftemp);
//
//		if (debug && debug_obs >= 0) {
//			bart_logging::log_draw_state(*this, j, "Before removing", debug_obs, allfit[debug_obs]);
//			bart_logging::log_draw_state(*this, j, "Current ftemp", debug_obs, ftemp[debug_obs]);
//			tree_preds.push_back(ftemp[debug_obs]);  // Store prediction for logging
//		}
//
//		for (size_t k = 0; k < n; k++) {
//			allfit[k] = allfit[k] - ftemp[k];
//			r[k] = y[k] - allfit[k];
//		}
//
//		if (debug && debug_obs >= 0) {
//			bart_logging::log_draw_state(*this, j, "After removing tree", debug_obs, allfit[debug_obs]);
//		}
//
//		// Do BART updates
//		bd(t[j], xi, di, pi, sigma, nv, pv, aug, gen);
//		drmu(t[j], xi, di, pi, sigma, gen);
//
//		// Add back updated tree's contribution
//		fit(t[j], xi, p, n, x, ftemp);
//
//		if (debug && debug_obs >= 0) {
//			bart_logging::log_draw_state(*this, j, "Before adding updated tree", debug_obs, ftemp[debug_obs]);
//
//			// Get detailed node information
//			double* xx = x + debug_obs * p;
//			tree::tree_p bottom_node = t[j].bn(xx, xi);
//			bart_logging::log_tree_prediction(*this, j, debug_obs,
//				bottom_node->gettheta(),
//				bottom_node->nid(),
//				ftemp[debug_obs]);
//		}
//
//		for (size_t k = 0; k < n; k++) allfit[k] += ftemp[k];
//
//		if (debug && debug_obs >= 0) {
//			bart_logging::log_draw_state(*this, j, "After adding updated tree", debug_obs, allfit[debug_obs]);
//		}
//	}
//
//	// Log summary of all tree predictions at end of iteration
//	if (debug && debug_obs >= 0) {
//		bart_logging::log_prediction_details(*this, debug_obs, tree_preds, allfit[debug_obs]);
//	}
//
//	if (dartOn) {
//		draw_s(nv, lpv, theta, gen);
//		draw_theta0(const_theta, theta, lpv, a, b, rho, gen);
//		for (size_t j = 0; j < p; j++) pv[j] = ::exp(lpv[j]);
//	}
//}

//--------------------------------------------------
void bart::updateY(double* newy) {

	y = newy;    // Update main y pointer

    if (r) delete[] r; // update residuals
    r = new double[n];
    di.y = r;

}
//--------------------------------------------------
//public functions
void bart::pr() //print to screen
{
   cout << "*****bart object:\n";
   cout << "m: " << m << std::endl;
   cout << "t[0]:\n " << t[0] << std::endl;
   cout << "t[m-1]:\n " << t[m-1] << std::endl;
   cout << "prior and mcmc info:\n";
   pi.pr();
   if(dart){
     cout << "*****dart prior (On):\n";
     cout << "a: " << a << std::endl;
     cout << "b: " << b << std::endl;
     cout << "rho: " << rho << std::endl;
     cout << "augmentation: " << aug << std::endl;
   }
   else cout << "*****dart prior (Off):\n";
   if(p) cout << "data set: n,p: " << n << ", " << p << std::endl;
   else cout << "data not set\n";
}

//void print_pred(size_t i) const {
//    if (i >= n) {
//        printf("Error: Index %zu out of bounds (n=%zu)\n", i, n);
//        return;
//    }
//    printf("Current prediction for obs %zu: %f\n", i, allfit[i]);
//}