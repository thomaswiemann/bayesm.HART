#include "bayesm.HART.h"

std::vector<bart> initializeBART(
	double* pZt,  // Z^T, column-major
	size_t nlgt,
	size_t nz,
	std::vector<double*> const& pstd_oldbetas_cols,  // one response buffer per coordinate
	List const& bart_params,
	arn gen)
{
	size_t num_trees = bart_params["num_trees"];
	double power = bart_params["power"];
	double base = bart_params["base"];
	double tau = bart_params["tau"];
	size_t numcut = bart_params["numcut"];

	bool sparse = bart_params["sparse"];
	double theta = bart_params["theta"];
	double omega = bart_params["omega"];
	double a = bart_params["a"];
	double b = bart_params["b"];
	double rho = bart_params["rho"];
	bool aug = bart_params["aug"];

	size_t nmin = as<size_t>(bart_params["nmin"]);

	size_t p = pstd_oldbetas_cols.size();
	std::vector<bart> models(p);

	for (size_t i = 0; i < p; i++) {
		models[i] = bart(num_trees);
		// `tau` is the ensemble scale: leaf coefficients are N(0, tau^2 / H).
		models[i].setprior(base, power, tau / std::sqrt((double)num_trees));
		models[i].setdata(nz, nlgt, pZt, pstd_oldbetas_cols[i], numcut);

		models[i].setdart(a, b, rho, aug, sparse, theta, omega);
		models[i].setleafguard(nmin);
	}

	return models;
}
