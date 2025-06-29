#include "bayesm.HART.h"

//Used in rhierLinearModel, rhierLinearMixture and rhierMnlRWMixture------------------------------------------------------
mat drawDelta(mat const& x, mat const& y, vec const& z, List const& comps, vec const& deltabar, mat const& Ad) {

	// Wayne Taylor 10/01/2014

	// delta = vec(D)
	//  given z and comps (z[i] gives component indicator for the ith observation, 
	//   comps is a list of mu and rooti)
	// y is n x p
	// x is n x k
	// y = xD' + U , rows of U are indep with covs Sigma_i given by z and comps

	int p = y.n_cols;
	int k = x.n_cols;
	int ncomp = comps.length();
	mat xtx = zeros<mat>(k * p, k * p);
	mat xty = zeros<mat>(p, k); //this is the unvecced version, reshaped after the sum

	//Create the index vectors, the colAll vectors are equal to span::all but with uvecs (as required by .submat)
	uvec colAlly(p), colAllx(k);
	for (int i = 0; i < p; i++) colAlly(i) = i;
	for (int i = 0; i < k; i++) colAllx(i) = i;

	//Loop through the components
	for (int compi = 0; compi < ncomp; compi++) {

		//Create an index vector ind, to be used like y[ind,]
		uvec ind = find(z == (compi + 1));

		//If there are observations in this component
		if (ind.size() > 0) {
			mat yi = y.submat(ind, colAlly);
			mat xi = x.submat(ind, colAllx);

			List compsi = comps[compi];
			rowvec mui = as<rowvec>(compsi[0]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
			mat rootii = trimatu(as<mat>(compsi[1])); //trimatu interprets the matrix as upper triangular
			yi.each_row() -= mui; //subtracts mui from each row of yi
			mat sigi = rootii * trans(rootii);
			xtx = xtx + kron(trans(xi) * xi, sigi);
			xty = xty + (sigi * (trans(yi) * xi));
		}
	}
	xty.reshape(xty.n_rows * xty.n_cols, 1);

	//vec(t(D)) ~ N(V^{-1}(xty + Ad*deltabar),V^{-1}) where V = (xtx+Ad)
	// compute the inverse of xtx+Ad
	mat ucholinv = solve(trimatu(chol(xtx + Ad)), eye(k * p, k * p)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
	mat Vinv = ucholinv * trans(ucholinv);

	return(Vinv * (xty + Ad * deltabar) + trans(chol(Vinv)) * as<vec>(rnorm(deltabar.size())));
}