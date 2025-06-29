#include "bayesm.HART.h"

std::vector<bart> initializeBART(
	double* pZt,  // Pointer to transposed Z matrix
	size_t nlgt,  // Number of observations
	size_t nz,    // Number of Z variables
	std::vector<double*> const& pstd_oldbetas_cols,  // Pointers to standardized beta columns
	List const& bart_params,
	arn gen)
{
	// Extract BART parameters
	size_t num_trees = bart_params["num_trees"];
	double power = bart_params["power"];
	double base = bart_params["base"];
	double tau = bart_params["tau"];
	size_t numcut = bart_params["numcut"];

	// Extract Dirichlet prior parameters
	bool sparse = bart_params["sparse"];
	double theta = bart_params["theta"];
	double omega = bart_params["omega"];
	double a = bart_params["a"];
	double b = bart_params["b"];
	double rho = bart_params["rho"];
	bool aug = bart_params["aug"];

	// Initialize models vector (p is length of pstd_oldbetas_cols)
	size_t p = pstd_oldbetas_cols.size();
	std::vector<bart> models(p);

	// Initialize each model
	for (int i = 0; i < p; i++) {
		models[i] = bart(num_trees);
		models[i].setprior(base, power, tau);
		models[i].setdata(nz, nlgt, pZt, pstd_oldbetas_cols[i], numcut);

		// Use dart
		models[i].setdart(a, b, rho, aug, sparse, theta, omega);
	}

	return models;
}

void update_stdoldbetas(mat const& oldbetas, std::vector<double*>& pstd_oldbetas_cols,
	vec const& ind, List const& comps) {

	// Get dimensions
	int n = oldbetas.n_rows;
	int p = oldbetas.n_cols;
	int ncomp = comps.length();

	// Create column indices
	uvec colAll(p);
	for (int i = 0; i < p; i++) colAll(i) = i;

	// Standardize coefficients by mixture component
	for (int compi = 0; compi < ncomp; compi++) {
		uvec indx = find(ind == (compi + 1));

		if (indx.size() > 0) {
			// Get component parameters
			List compsi = comps[compi];
			rowvec mui = as<rowvec>(compsi[0]);
			mat rootii = trimatu(as<mat>(compsi[1]));

			// Get submatrix for this component
			mat betai = oldbetas.submat(indx, colAll);

			// Center and standardize
			mat betai_centered = betai.each_row() - mui;
			mat betai_std = betai_centered * rootii;
			//mat betai_std = (rootii * betai_centered.t()).t();


			// Update values through pointers
			for (int j = 0; j < p; j++) {
				for (size_t k = 0; k < indx.size(); k++) {
					pstd_oldbetas_cols[j][indx(k)] = betai_std(k, j);
				}
			}
		}
	}
}