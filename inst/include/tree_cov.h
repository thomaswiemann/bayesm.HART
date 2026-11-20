#ifndef TREE_COV_H
#define TREE_COV_H

#include <sstream>
#include "direct_cov.h"
#include "cov_helpers.h"
#include "varbart.h"
#include "vcbart.h"

// Covariate-dependent Sigma(Z) through the modified Cholesky
// Sigma_i^{-1} = T_i^T V_i^{-1} T_i (appendix_sampler, app:hs_model):
//
//   v_j(Z)    = product of H_v trees,    s ~ IG(a_v, b_vj),
//   phi_jk(Z) = sum of H_phi trees,      leaves N(0, tau_phi^2 / H_phi).
//
// draw() runs the blocks in sweep order: residuals r = T e, the variance
// ensembles, the off-diagonal ensembles, then the per-unit precision.
class TreeCov : public DirectCov {
public:
    TreeCov(Rcpp::List const& var_params, Rcpp::List const& phi_params);

    void init(int n, int D, int nz, double* pZt, int nkeep) override;
    void start(arma::mat const& E) override;
    void start_dart(int rep) override;
    void draw(arma::mat const& E, rn& gen) override;

    std::vector<arma::mat> const& precisions() const override { return Q_; }
    arma::mat const& root(int i) const override { return rootpi_[i]; }

    void store(int mkeep, arma::vec const& mu, bool store_trees,
               UniqueZCache& zc) override;
    void pack(Rcpp::List& out, bool store_trees) override;
    bool covariate_sigma() const override { return true; }

protected:
    Rcpp::List var_params_, phi_params_;
    int        n_, D_, nz_;
    bool       full_cov_;                    // initialized before Hphi_
    size_t     Hv_, Hphi_;
    bool       sparse_v_, sparse_phi_;
    int        burn_v_, burn_phi_;

    std::vector<varbart>             v_;     // v_j; y pointers into r_sq_
    std::vector<std::vector<vcbart>> phi_;   // phi_[j][k], k < j; y into Y_phi_
    arma::mat  e_, r_sq_, vmat_;             // n x D
    arma::mat  Y_phi_;                       // n x D(D-1)/2
    std::vector<arma::mat> Q_, rootpi_;      // n, each D x D

    arma::cube var_varcount_, var_varprob_;  // nz x D x nkeep
    std::vector<std::stringstream> var_treess_;
    std::vector<std::vector<std::stringstream>> phi_treess_;

    void refresh_residuals(arma::mat const& E);
    void draw_variance_trees(rn& gen);
    void draw_phi_trees(rn& gen);
    void refresh_precision();

    cov::cov_eval eval(int i) const {
        return cov::cov_eval{v_, phi_, (size_t)i};
    }
    static size_t pair_index(size_t j, size_t k) { return j * (j - 1) / 2 + k; }
};

#endif
