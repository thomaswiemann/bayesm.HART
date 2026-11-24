#ifndef CONSTANT_COV_H
#define CONSTANT_COV_H

#include "direct_cov.h"

// One Sigma ~ IW(nu, Psi) shared by all units, drawn from its conjugate
// conditional Sigma | E ~ IW(nu + n, Psi + E^T E).  Output keeps bayesm's
// one-component `nmix` schema so the prediction path is shared.
class ConstantCov : public DirectCov {
public:
    ConstantCov(double nu, arma::mat const& Psi);

    void init(int n, int D, int nz, double* pZt, int nkeep) override;
    void start(arma::mat const& E) override;
    void start_dart(int) override {}
    void draw(arma::mat const& E, rn& gen) override;

    std::vector<arma::mat> const& precisions() const override { return Q_; }
    arma::mat const& root(int) const override { return rooti_; }

    void store(int mkeep, arma::vec const& mu, bool store_trees,
               UniqueZCache& zc) override;
    void pack(Rcpp::List& out, bool store_trees) override;
    Rcpp::List checkpoint() override;
    void restore(Rcpp::List const& state) override;
    bool covariate_sigma() const override { return false; }

protected:
    double     nu_;
    arma::mat  Psi_;
    int        n_, D_;
    arma::mat  Sigma_, rooti_;           // rooti rooti^T = Sigma^{-1}, upper
    std::vector<arma::mat> Q_;           // n copies of Sigma^{-1}
    arma::mat  probdraw_;
    Rcpp::List compdraw_;

    void set_precision();
    void refresh_root();
};

#endif
