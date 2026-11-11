#ifndef DIRECT_HART_BACKEND_H
#define DIRECT_HART_BACKEND_H

#include <sstream>
#include "hart_backend.h"
#include "unique_z_cache.h"
#include "vecbart.h"
#include "cov_helpers.h"

// Direct HART (model B):  theta_i = mu + f(Z_i) + S_i eps_i,
// Sigma_i = T_i^{-1} V_i T_i^{-T} (appendix_sampler, app:hs_model).
//
//   mu ~ N(mubar, Amu^{-1}),
//   f(Z) = sum of H_f vector-leaf trees,  m_hg ~ N_D(0, diag(tau^2) / H_f),
//   tau_d^2 ~ IG(a_tau, b_tau_d),         b_tau_d = (a_tau - 1) tau0_d^2,
//   v_j(Z) = product of H_v trees,        s ~ IG(a_v, b_vj),
//   phi_jk(Z) = sum of H_phi trees,       leaves N(0, tau_phi^2 / H_phi).
//
// One iteration (app:hs_sweep blocks 2-4; the driver runs block 1):
//   draw_mu()               mu | theta, f, Q               (eq:hs_mu)
//   refresh_residuals()     e = theta - mu - f,  r = T e    (eq:hs_residuals)
//   draw_variance_trees()   v_j | r_j                       (eq:hs_variance_leaf)
//   draw_phi_trees()        phi_jk | e, v_j                 (eq:hs_phi_leaf)
//   refresh_precision()     Q_i = T_i^T V_i^{-1} T_i, rootpi_i = T_i^T V_i^{-1/2}
//   draw_mean_trees()       f | theta - mu, Q, tau^2        (eq:hs_mean_leaf)
//   draw_mean_scales()      tau_d^2 | f                     (eq:hs_scale_update)
//   refresh_kernels()       prior means mu + f_i
class DirectHartBackend : public HartBackend {
public:
    DirectHartBackend(arma::mat const& theta, arma::mat const& Z,
                      arma::vec const& mubar, arma::mat const& Amu,
                      int R, int keep, Rcpp::List const& bart_params,
                      Rcpp::List const& var_params,
                      Rcpp::List const& phi_params);

    void draw_iteration(int rep) override;
    arma::mat const& precision_root(int i) const override { return rootpi_[i]; }
    void store(int mkeep) override;
    BackendPackPayload pack() override;
    BackendCapabilities capabilities() const override {
        return {false, true, true};
    }

protected:
    // Priors
    arma::vec  mubar_;
    arma::mat  Amu_;
    double     a_tau_;
    arma::vec  b_tau_;
    size_t     Hf_, Hv_, Hphi_;
    bool       full_cov_;
    bool       sparse_f_, sparse_v_, sparse_phi_;
    int        burn_f_, burn_v_, burn_phi_;

    // State
    arma::vec  mu_, tau2_;
    vecbart    f_;
    std::vector<varbart>             v_;     // v_j; y pointers into r_sq_
    std::vector<std::vector<vcbart>> phi_;   // phi_[j][k], k < j; y into Y_phi_
    arma::mat  e_, r_sq_, vmat_;             // n x D
    arma::mat  Y_phi_;                       // n x D(D-1)/2
    std::vector<arma::mat> Q_, rootpi_;      // n, each D x D

    // Draws
    arma::mat  mu_draw_, scale_draw_;        // nkeep x D
    arma::mat  varcount_, varprob_;          // nz x nkeep
    arma::cube var_varcount_, var_varprob_;  // nz x D x nkeep
    std::stringstream treess_;
    std::vector<std::stringstream> var_treess_;
    std::vector<std::vector<std::stringstream>> phi_treess_;
    UniqueZCache zcache_;

    // Blocks 2-4 with each draw switchable (the test probe pins blocks).
    void sweep(bool mu, bool v, bool phi, bool f, bool tau);
    void draw_mu();
    void refresh_residuals();
    void draw_variance_trees();
    void draw_phi_trees();
    void refresh_precision();
    void draw_mean_trees();
    void draw_mean_scales();
    void refresh_kernels();

    cov::cov_eval eval(int i) const {
        return cov::cov_eval{v_, phi_, (size_t)i};
    }
    static size_t pair_index(size_t j, size_t k) { return j * (j - 1) / 2 + k; }
    template <class B>
    static void setdart_from(B& model, Rcpp::List const& p) {
        model.setdart(Rcpp::as<double>(p["a"]), Rcpp::as<double>(p["b"]),
                      Rcpp::as<double>(p["rho"]), Rcpp::as<bool>(p["aug"]),
                      Rcpp::as<bool>(p["sparse"]), Rcpp::as<double>(p["theta"]),
                      Rcpp::as<double>(p["omega"]));
    }
};

#endif
