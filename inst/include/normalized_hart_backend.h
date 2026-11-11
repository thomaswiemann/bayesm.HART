#ifndef NORMALIZED_HART_BACKEND_H
#define NORMALIZED_HART_BACKEND_H

#include <sstream>
#include "hart_backend.h"
#include "unique_z_cache.h"

// Normalized HART (model A):  theta_i = mu + L^{-1} { f*(Z_i) + eps_i }.
//
//   mu ~ N(mubar, Amu^{-1}),   Sigma = (L^T L)^{-1} ~ IW(nu, Psi),
//   f*_d(Z) = sum of H scalar trees per coordinate on the whitened scale,
//   leaf sd tau / sqrt(H), fixed.
//
// One iteration (appendix_models blocks 2-4; the driver runs block 1):
//   draw_mu()          mu | L, f*, theta                  (sigma_gibbs.h)
//   draw_cholesky()    L  | mu, f*, theta                 (sigma_gibbs.h)
//   draw_trees()       theta~ = L (theta - mu); backfit each f*_d
//   refresh_kernels()  delta_Z = L^{-1} f*, prior means, rooti = L^T
//
// Output keeps bayesm's one-component `nmix` schema (probdraw, compdraw with
// mu and rooti), so prediction code is shared with the bayesm paths.
class NormalizedHartBackend : public HartBackend {
public:
    NormalizedHartBackend(arma::mat const& theta, arma::mat const& Z,
                          arma::vec const& mubar, arma::mat const& Amu,
                          double nu, arma::mat const& Psi,
                          int R, int keep, Rcpp::List const& bart_params);

    void draw_iteration(int rep) override;
    arma::mat const& precision_root(int) const override { return rooti_; }
    void store(int mkeep) override;
    BackendPackPayload pack() override;
    BackendCapabilities capabilities() const override {
        return {true, false, true};
    }

protected:
    // Priors
    arma::vec  mubar_;
    arma::mat  Amu_, Psi_;
    double     nu_;
    Rcpp::List bart_params_;
    bool       sparse_;
    int        burn_;

    // State
    arma::vec  mu_;
    arma::mat  L_, rooti_;              // L lower, L^T L = Sigma^{-1}; rooti = L^T
    std::vector<bart> f_;               // D ensembles; y pointers into theta_std_
    arma::mat  theta_std_;              // n x D whitened responses (stable)
    std::vector<double*> theta_std_cols_;
    arma::mat  fstar_, delta_Z_;        // n x D

    // Draws
    arma::mat  probdraw_;
    Rcpp::List compdraw_;
    std::vector<std::stringstream> treess_;
    arma::cube varcount_, varprob_;
    UniqueZCache zcache_;

    void draw_mu();
    void draw_cholesky();
    void draw_trees(int rep);
    void refresh_kernels();
};

#endif
