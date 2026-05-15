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
	for (size_t i = 0; i < p; i++) {
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

//--------------------------------------------------
// Heteroscedastic standardization: tilde_theta_i = D(Z_i)^{-1/2} L(Z_i) (theta_i - mu).
// Single function handles both the diagonal case (Phase 1; pass empty
// phi_models) and the full Cholesky case (Phase 2; pass populated phi_models)
// by routing every per-unit operation through cov::mahalanobis.
void update_stdoldbetas_het(
	mat const& oldbetas,
	std::vector<double*>& pstd_oldbetas_cols,
	vec const& mu,
	std::vector<varbart> const& d_models,
	std::vector<std::vector<heterbart>> const& phi_models)
{
	const size_t nlgt = oldbetas.n_rows;
	const size_t p    = oldbetas.n_cols;
	if (mu.n_elem != p)            stop("update_stdoldbetas_het: mu has wrong length");
	if (d_models.size() != p)      stop("update_stdoldbetas_het: d_models has wrong length");
	if (pstd_oldbetas_cols.size() != p) stop("update_stdoldbetas_het: pstd_oldbetas_cols has wrong length");

	for (size_t i = 0; i < nlgt; i++) {
		cov::cov_eval ev{d_models, phi_models, i};
		vec ttheta = cov::mahalanobis(oldbetas.row(i).t(), mu, ev);
		for (size_t j = 0; j < p; j++) {
			pstd_oldbetas_cols[j][i] = ttheta[j];
		}
	}
}

//--------------------------------------------------
// Construct D variance-tree ensembles (one per dimension).
//
// var_params carries:
//   num_trees   = m'  (trees per dimension)
//   power       = mybeta in pinfo
//   base        = alpha  in pinfo
//   nu          = baseline degrees of freedom (NOT per-tree)
//   lambda      = baseline scale^2            (NOT per-tree)
//   numcut      = (optional, default 100) cutpoints per Z covariate
//   sparse,a,b,rho,aug,theta,omega = optional DART hyperparameters
//
// Pratola per-tree calibration is applied internally:
//   nu'      = 2 / (1 - (1 - 2/nu)^{1/m'})
//   lambda'  = lambda^{1/m'}
//
// peta_sq_cols[j] points to a length-nlgt buffer holding the current
// squared innovations eta_i^{(j)2} for dimension j.  The buffer is owned
// by the calling MCMC loop, which mutates it in place each iteration.
std::vector<varbart> initializeVarBART(
	double* pZt, size_t nlgt, size_t nz,
	std::vector<double*> const& peta_sq_cols,
	List const& var_params, arn /*gen*/)
{
	size_t num_trees = as<size_t>(var_params["num_trees"]);
	double power     = as<double>(var_params["power"]);
	double base      = as<double>(var_params["base"]);
	double nu        = as<double>(var_params["nu"]);
	double lambda    = as<double>(var_params["lambda"]);
	size_t numcut    = var_params.containsElementNamed("numcut")
	                   ? as<size_t>(var_params["numcut"]) : 100;

	// Optional DART parameters (mirroring initializeBART defaults via R wrapper).
	bool   sparse = var_params.containsElementNamed("sparse") ? as<bool>(var_params["sparse"]) : false;
	double dart_a = var_params.containsElementNamed("a")      ? as<double>(var_params["a"])    : 0.5;
	double dart_b = var_params.containsElementNamed("b")      ? as<double>(var_params["b"])    : 1.0;
	double dart_r = var_params.containsElementNamed("rho")    ? as<double>(var_params["rho"])  : (double)nz;
	bool   aug    = var_params.containsElementNamed("aug")    ? as<bool>(var_params["aug"])    : false;
	double dart_t = var_params.containsElementNamed("theta")  ? as<double>(var_params["theta"]): 0.0;
	double dart_o = var_params.containsElementNamed("omega")  ? as<double>(var_params["omega"]): 1.0;

	// Pratola per-tree calibration (see discussion §1.6).
	if (nu <= 2.0)        stop("initializeVarBART: nu must be > 2");
	if (lambda <= 0.0)    stop("initializeVarBART: lambda must be > 0");
	if (num_trees == 0)   stop("initializeVarBART: num_trees must be >= 1");

	double inv_m   = 1.0 / static_cast<double>(num_trees);
	double nu_p    = 2.0 / (1.0 - std::pow(1.0 - 2.0 / nu, inv_m));
	double lambda_p = std::pow(lambda, inv_m);

	const size_t D = peta_sq_cols.size();
	std::vector<varbart> models(D);
	for (size_t j = 0; j < D; j++) {
		models[j] = varbart(num_trees);
		models[j].setprior(base, power, /*tau unused for varbart*/ 1.0);
		models[j].setdata(nz, nlgt, pZt, peta_sq_cols[j], numcut);
		models[j].setvarprior(nu_p, lambda_p);
		models[j].setdart(dart_a, dart_b, dart_r, aug, sparse, dart_t, dart_o);
	}
	return models;
}

//--------------------------------------------------
// Construct D(D-1)/2 phi-tree heterbart ensembles for the modified-Cholesky
// off-diagonal entries phi_{jk}(.), j > k.  Returned object is jagged:
// phi_models[j].size() == j (so phi_models[0] is empty, phi_models[1] has
// length 1, etc.).
//
// Per discussions/2026-05-12-phibart-vs-heterbart.md, phi_jk(z) is fit via a
// varying-coefficient -> heteroscedastic mean-regression transform that lets
// us reuse heterbart unchanged for the leaf draw and birth/death MH.  The
// only adjustment is the leaf-size guard: raw counts are not informative
// about phi when c_i^{(k)} is small, so we add an effective-sample-size
// floor (default ess_min = 5; configurable via phi_params$ess_min).
std::vector<std::vector<heterbart>> initializePhiBART(
	double* pZt, size_t nlgt, size_t nz,
	std::vector<std::vector<double*>> const& py_cols,
	List const& phi_params, arn /*gen*/)
{
	size_t num_trees = as<size_t>(phi_params["num_trees"]);
	double power     = as<double>(phi_params["power"]);
	double base      = as<double>(phi_params["base"]);
	double tau       = as<double>(phi_params["tau"]);
	size_t numcut    = phi_params.containsElementNamed("numcut")
	                   ? as<size_t>(phi_params["numcut"]) : 100;

	bool   sparse = phi_params.containsElementNamed("sparse") ? as<bool>(phi_params["sparse"]) : false;
	double dart_a = phi_params.containsElementNamed("a")      ? as<double>(phi_params["a"])    : 0.5;
	double dart_b = phi_params.containsElementNamed("b")      ? as<double>(phi_params["b"])    : 1.0;
	double dart_r = phi_params.containsElementNamed("rho")    ? as<double>(phi_params["rho"])  : (double)nz;
	bool   aug    = phi_params.containsElementNamed("aug")    ? as<bool>(phi_params["aug"])    : false;
	double dart_t = phi_params.containsElementNamed("theta")  ? as<double>(phi_params["theta"]): 0.0;
	double dart_o = phi_params.containsElementNamed("omega")  ? as<double>(phi_params["omega"]): 1.0;

	// Leaf-size guards: nmin defaults to 2 (looser than mean-tree default of
	// 5 to allow finer phi splits), ess_min defaults to 5.0.
	size_t nmin    = phi_params.containsElementNamed("nmin")    ? as<size_t>(phi_params["nmin"])    : 2;
	double ess_min = phi_params.containsElementNamed("ess_min") ? as<double>(phi_params["ess_min"]) : 5.0;

	const size_t D = py_cols.size();
	std::vector<std::vector<heterbart>> models(D);
	for (size_t j = 1; j < D; j++) {
		if (py_cols[j].size() != j)
			stop("initializePhiBART: py_cols[j] must have length j");
		models[j] = std::vector<heterbart>(j);
		for (size_t k = 0; k < j; k++) {
			models[j][k] = heterbart(num_trees);
			models[j][k].setprior(base, power, tau);
			models[j][k].setdata(nz, nlgt, pZt, py_cols[j][k], numcut);
			models[j][k].setdart(dart_a, dart_b, dart_r, aug, sparse, dart_t, dart_o);
			models[j][k].setleafguard(nmin, ess_min);
		}
	}
	return models;
}

//--------------------------------------------------
// Conjugate normal draw of the single posterior mean mu under heteroscedastic
// Sigma(Z_i).  Routes every Sigma(Z_i)^{-1} access through cov::rootpi so the
// diagonal (phi_models empty) and full-Cholesky cases share one code path.
//
// Math:
//   V_post^{-1} = Amu + sum_i Sigma(Z_i)^{-1}
//   m_post     = V_post (Amu mubar0 + sum_i Sigma(Z_i)^{-1} theta_i)
//   mu | .     ~ N(m_post, V_post).
//
// Implementation accumulates Sigma(Z_i)^{-1} via Sigma(Z_i)^{-1} = R_i R_i^T
// where R_i = cov::rootpi(.) (upper-triangular).
vec drawMuHeterCov(
	mat const& oldbetas,
	vec const& mubar0,
	mat const& Amu,
	std::vector<varbart> const& d_models,
	std::vector<std::vector<heterbart>> const& phi_models,
	arn& gen)
{
	const size_t nlgt = oldbetas.n_rows;
	const size_t p    = oldbetas.n_cols;
	if (mubar0.n_elem != p)                stop("drawMuHeterCov: mubar0 has wrong length");
	if (d_models.size() != p)              stop("drawMuHeterCov: d_models has wrong length");

	// Accept either a full p x p prior precision Amu, or a 1 x 1 scalar
	// precision (the default convention used by rhierMnlRwMixture's wrapper);
	// in the scalar case we expand to scalar * I_p.
	mat Amu_full;
	if (Amu.n_rows == p && Amu.n_cols == p) {
		Amu_full = Amu;
	} else if (Amu.n_rows == 1 && Amu.n_cols == 1) {
		Amu_full = Amu(0, 0) * eye<mat>(p, p);
	} else {
		stop("drawMuHeterCov: Amu must be p x p or 1 x 1 (scalar precision)");
	}

	mat V_inv = Amu_full;
	vec m_num = Amu_full * mubar0;

	for (size_t i = 0; i < nlgt; i++) {
		cov::cov_eval ev{d_models, phi_models, i};
		mat R = cov::rootpi(ev);                  // upper-tri, R R^T = Sigma(Z_i)^{-1}
		mat Sinv_i = R * R.t();                   // p x p, symmetric PD
		V_inv += Sinv_i;
		m_num += Sinv_i * oldbetas.row(i).t();
	}

	// V_inv = R^T R  with R upper-triangular  =>  V_post = R^{-1} R^{-T}.
	mat R = chol(V_inv);                          // upper-triangular
	// m_post = V_post m_num = R^{-1} R^{-T} m_num   (two triangular solves).
	vec y      = solve(trimatl(R.t()), m_num);    // R^T y = m_num
	vec m_post = solve(trimatu(R),     y);        // R m_post = y

	// Sample mu ~ N(m_post, V_post).  Use R^{-1} as a square root of V_post:
	//   (R^{-1})(R^{-1})^T = R^{-1} R^{-T} = V_post.  z ~ N(0, I)  =>
	//   mu = m_post + R^{-1} z.
	vec z(p);
	for (size_t j = 0; j < p; j++) z[j] = gen.normal();
	vec mu = m_post + solve(trimatu(R), z);
	return mu;
}