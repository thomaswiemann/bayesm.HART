#include "mnl_theta_updater.h"

using namespace arma;
using namespace Rcpp;

MnlThetaUpdaterAdapter::MnlThetaUpdaterAdapter(
    mat& oldbetas,
    std::vector<moments> const& lgtdata,
    vec const& SignRes,
    double s, int R)
    : oldbetas_(oldbetas), lgtdata_(lgtdata), SignRes_(SignRes),
      nacceptbeta_(0), s_(s),
      nlgt_(oldbetas.n_rows), nvar_(oldbetas.n_cols), R_(R)
{
    // Initialize log-likelihoods (deterministic, no RNG impact)
    oldll_.zeros(nlgt_);
    hess_plus_prec_work_.set_size(nvar_, nvar_);
    incroot_work_.set_size(nvar_, nvar_);
    eye_work_.eye(nvar_, nvar_);
    for (int i = 0; i < nlgt_; i++) {
        oldll_[i] = llmnl_con(vectorise(oldbetas_(i, span::all)),
                               lgtdata_[i].y, lgtdata_[i].X, SignRes_);
    }
}

void MnlThetaUpdaterAdapter::update_unit(
    int i, mat const& rootpi_i, vec const& prior_mean_i)
{
    hess_plus_prec_work_ = lgtdata_[i].hess + rootpi_i * trans(rootpi_i);
    incroot_work_ = solve(trimatu(chol(hess_plus_prec_work_)), eye_work_);
    hess_plus_prec_work_ = incroot_work_ * trans(incroot_work_);
    incroot_work_ = chol(hess_plus_prec_work_);

    mnlMetropOnceOut out = mnlMetropOnce_con(
        lgtdata_[i].y, lgtdata_[i].X,
        vectorise(oldbetas_(i, span::all)),
        oldll_[i], s_, incroot_work_, prior_mean_i, rootpi_i, SignRes_);

    if (out.stay == 0) nacceptbeta_++;
    oldbetas_(i, span::all) = trans(out.betadraw);
    oldll_[i] = out.oldll;
}

void MnlThetaUpdaterAdapter::store(int /* mkeep */) {}

UpdaterPackPayload MnlThetaUpdaterAdapter::pack() {
    return {List::create()};
}

AcceptanceMetrics MnlThetaUpdaterAdapter::acceptance_metrics() const {
    AcceptanceMetrics am;
    am.acceptrbeta = nacceptbeta_ / (R_ * nlgt_ * 1.0) * 100;
    return am;
}

double MnlThetaUpdaterAdapter::current_loglike() const {
    return sum(oldll_);
}
