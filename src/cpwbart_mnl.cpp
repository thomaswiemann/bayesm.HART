#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "tree.h"
#include "treefuns.h"
#include "bayesm.HART.h"

typedef std::vector<tree> vtree;

//[[Rcpp::export]]
Rcpp::NumericMatrix cpwbart_mnl(
	std::string const& trees,
	Rcpp::List const& cutpoints,
	arma::mat const& X
) {
	using namespace arma;
	using namespace Rcpp;

	// Process cutpoints into xinfo
	size_t p = cutpoints.size();
	xinfo xi;
	xi.resize(p);
	for (size_t i = 0; i < p; i++) {
		Rcpp::NumericVector cp = cutpoints[i];
		xi[i].resize(cp.size());
		std::copy(cp.begin(), cp.end(), xi[i].begin());
	}

	// Process tree strings
	std::stringstream treess(trees);

	// Read header
	size_t nd, m, p_check;
	treess >> nd >> m >> p_check;

	// Verify dimensions match
	if (p_check != p) {
		stop("Dimension mismatch between trees and cutpoints");
	}

	// Process new x matrix
	size_t np = X.n_rows;

	// Create a non-const copy of X for prediction
	double* px = new double[X.n_elem];
	std::memcpy(px, X.memptr(), X.n_elem * sizeof(double));

	if (X.n_cols != p) {
		delete[] px;
		stop("Dimension mismatch between new data and cutpoints");
	}

	// Storage for trees
	std::vector<vtree> tmat(nd);
	for (size_t i = 0; i < nd; i++) {
		tmat[i].resize(m);
	}

	// Read in all trees
	for (size_t i = 0; i < nd; i++) {
		for (size_t j = 0; j < m; j++) {
			treess >> tmat[i][j];
		}
	}

	// Storage for predictions and temporary values
	// mat yhat = zeros<mat>(nd, np);
	Rcpp::NumericMatrix yhat(nd, np);
	double* fptemp = new double[np];

	// Make predictions for each draw
	for (size_t i = 0; i < nd; i++) {
		for (size_t j = 0; j < m; j++) {
			fit(tmat[i][j], xi, p, np, px, fptemp);
			for (size_t k = 0; k < np; k++) {
				yhat(i, k) += fptemp[k];
			}
		}
		// Divide by number of trees to get average prediction for this draw
		for (size_t k = 0; k < np; k++) {
			yhat(i, k) /= m;
		}
	}

	// Clean up
	delete[] fptemp;
	delete[] px;

	return yhat;
		//Rcpp::List::create(
	//	Rcpp::Named("yhat.test") = yhat
	//);
}