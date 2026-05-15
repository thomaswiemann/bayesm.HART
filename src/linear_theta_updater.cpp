#include "linear_theta_updater.h"

using namespace arma;
using namespace Rcpp;

LinearThetaUpdaterAdapter::LinearThetaUpdaterAdapter(
    mat& oldbetas,
    std::vector<moments> const& regdata,
    vec tau, double nu_e, vec const& ssq,
    int R, int keep)
    : oldbetas_(oldbetas), regdata_(regdata),
      tau_(tau), nu_e_(nu_e), ssq_(ssq),
      nreg_(oldbetas.n_rows), nvar_(oldbetas.n_cols), R_(R), keep_(keep)
{
    taudraw_.zeros(R / keep, nreg_);
}

void LinearThetaUpdaterAdapter::update_unit(
    int i, mat const& rootpi_i, vec const& prior_mean_i)
{
    mat Abeta    = rootpi_i * trans(rootpi_i);
    vec Abetabar = Abeta * prior_mean_i;

    unireg out = runiregG(regdata_[i].y, regdata_[i].X,
                           regdata_[i].XpX, regdata_[i].Xpy,
                           tau_[i], Abeta, Abetabar, nu_e_, ssq_[i]);
    oldbetas_(i, span::all) = trans(out.beta);
    tau_[i] = out.sigmasq;
}

void LinearThetaUpdaterAdapter::store(int mkeep) {
    taudraw_(mkeep - 1, span::all) = trans(tau_);
}

UpdaterPackPayload LinearThetaUpdaterAdapter::pack() {
    return {List::create(Named("taudraw") = taudraw_)};
}

AcceptanceMetrics LinearThetaUpdaterAdapter::acceptance_metrics() const {
    return {};
}

double LinearThetaUpdaterAdapter::current_loglike() const {
    double ll = 0.0;
    for (int r = 0; r < nreg_; r++) {
        vec resid = regdata_[r].y - regdata_[r].X * trans(oldbetas_.row(r));
        int ni = resid.n_elem;
        ll += -0.5 * ni * log(2.0 * M_PI * tau_[r]) - 0.5 * dot(resid, resid) / tau_[r];
    }
    return ll;
}
