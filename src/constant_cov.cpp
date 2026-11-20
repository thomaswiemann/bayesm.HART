#include "constant_cov.h"
#include "bayesm.HART.h"

using namespace arma;
using namespace Rcpp;

ConstantCov::ConstantCov(double nu, mat const& Psi)
    : nu_(nu), Psi_(Psi), n_(0), D_(0)
{}

void ConstantCov::init(int n, int D, int /*nz*/, double* /*pZt*/, int nkeep) {
    n_ = n; D_ = D;
    if ((int)Psi_.n_rows != D_ || (int)Psi_.n_cols != D_)
        stop("ConstantCov: V must be nvar x nvar.");
    if (nu_ <= (double)(D_ - 1))
        stop("ConstantCov: nu must exceed nvar - 1.");
    Q_.assign(n_, mat(D_, D_, fill::zeros));
    probdraw_.ones(nkeep, 1);
    compdraw_ = List(nkeep);
}

void ConstantCov::start(mat const& /*E*/) {
    Sigma_ = Psi_ / nu_;                                    // prior mode scale
    refresh_root();
    set_precision();
}

// Sigma | theta, mu, f ~ IW(nu + n, Psi + E^T E),  E = theta - mu - f.
void ConstantCov::draw(mat const& E, rn& /*gen*/) {
    mat S = Psi_ + E.t() * E;
    List w = rwishart(nu_ + (double)n_, inv_sympd(S));
    Sigma_ = as<mat>(w["IW"]);
    refresh_root();
    set_precision();
}

void ConstantCov::refresh_root() {
    rooti_ = solve(trimatu(chol(Sigma_)), eye<mat>(D_, D_));
}

void ConstantCov::set_precision() {
    mat Q = rooti_ * rooti_.t();
    for (int i = 0; i < n_; ++i) Q_[i] = Q;
}

void ConstantCov::store(int mkeep, vec const& mu, bool /*store_trees*/,
                        UniqueZCache& /*zc*/) {
    compdraw_[mkeep - 1] = List::create(List::create(
        Named("mu")    = NumericVector(mu.begin(), mu.end()),
        Named("rooti") = rooti_));
}

void ConstantCov::pack(List& out, bool /*store_trees*/) {
    out["nmix"] = List::create(Named("probdraw") = probdraw_,
                               Named("zdraw")    = R_NilValue,
                               Named("compdraw") = compdraw_);
}

List ConstantCov::checkpoint() {
    return List::create(Named("Sigma") = Sigma_);
}

void ConstantCov::restore(List const& state) {
    Sigma_ = as<mat>(state["Sigma"]);
    if ((int)Sigma_.n_rows != D_ || (int)Sigma_.n_cols != D_)
        stop("ConstantCov: state dimension mismatch.");
    refresh_root();
    set_precision();
}
