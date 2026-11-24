#ifndef DIRECT_HART_BACKEND_H
#define DIRECT_HART_BACKEND_H

#include <memory>
#include <sstream>
#include "hart_backend.h"
#include "unique_z_cache.h"
#include "vecbart.h"
#include "direct_cov.h"

// Direct HART (model B):  theta_i = mu + f(Z_i) + eps_i, eps_i ~ N(0, Sigma_i)
// (appendix_sampler, app:hs_model).  The mean model lives here:
//
//   mu ~ N(mubar, Amu^{-1}),
//   f(Z) = sum of H_f vector-leaf trees,  m_hg ~ N_D(0, diag(tau^2) / H_f),
//   tau_d^2 ~ IG(a_tau, b_tau_d),         b_tau_d = (a_tau - 1) tau0_d^2.
//
// Sigma_i belongs to the DirectCov strategy: TreeCov for a covariate-dependent
// Sigma(Z_i), ConstantCov for one Sigma shared by all units.
//
// One iteration (app:hs_sweep blocks 2-4; the driver runs block 1):
//   draw_mu()            mu | theta, f, Q                (eq:hs_mu)
//   cov_->draw(E)        Sigma | E,  E = theta - mu - f  (block 3)
//   draw_mean_trees()    f | theta - mu, Q, tau^2        (eq:hs_mean_leaf)
//   draw_mean_scales()   tau_d^2 | f                     (eq:hs_scale_update)
//   refresh_kernels()    prior means mu + f_i
class DirectHartBackend : public HartBackend {
public:
    DirectHartBackend(arma::mat const& theta, arma::mat const& Z,
                      arma::vec const& mubar, arma::mat const& Amu,
                      int R, int keep, Rcpp::List const& bart_params,
                      std::unique_ptr<DirectCov> cov);

    void draw_iteration(int rep) override;
    arma::mat const& precision_root(int i) const override {
        return cov_->root(i);
    }
    void store(int mkeep) override;
    BackendPackPayload pack() override;
    Rcpp::List checkpoint() override;
    void restore(Rcpp::List const& state) override;
    BackendCapabilities capabilities() const override {
        bool cz = cov_->covariate_sigma();
        return {!cz, cz, true};
    }

protected:
    // Priors
    arma::vec  mubar_;
    arma::mat  Amu_;
    double     a_tau_;
    arma::vec  b_tau_;
    size_t     Hf_;
    bool       sparse_f_;
    int        burn_f_;

    // State
    std::unique_ptr<DirectCov> cov_;
    arma::vec  mu_, tau2_;
    vecbart    f_;
    arma::mat  E_;                           // n x D residuals theta - mu - f

    // Draws
    arma::mat  mu_draw_, scale_draw_;        // nkeep x D
    arma::mat  varcount_, varprob_;          // nz x nkeep
    std::stringstream treess_;
    UniqueZCache zcache_;

    // Mean blocks, each draw switchable (the test probe pins blocks).
    void sweep(bool mu, bool cov, bool f, bool tau);
    void draw_mu();
    void refresh_residuals();
    void draw_mean_trees();
    void draw_mean_scales();
    void refresh_kernels();
};

#endif
