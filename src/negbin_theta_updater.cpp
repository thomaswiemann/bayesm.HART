#include "negbin_theta_updater.h"

using namespace arma;
using namespace Rcpp;

NegbinThetaUpdaterAdapter::NegbinThetaUpdaterAdapter(
    mat& oldbetas,
    std::vector<moments> const& regdata,
    double sbeta, double alpha, double alphacroot,
    double a_prior, double b_prior, bool fixalpha,
    int R, int keep)
    : oldbetas_(oldbetas), regdata_(regdata),
      nacceptbeta_(0), nacceptalpha_(0),
      sbeta_(sbeta), alpha_(alpha), alphacroot_(alphacroot),
      a_prior_(a_prior), b_prior_(b_prior), fixalpha_(fixalpha),
      nreg_(oldbetas.n_rows), nvar_(oldbetas.n_cols), R_(R), keep_(keep)
{
    oldlpostbeta_.zeros(nreg_);
    clpostbeta_.zeros(nreg_);
    alphadraw_.zeros(R / keep);
    A_work_.set_size(nvar_, nvar_);
    U_work_.set_size(nvar_, nvar_);
    eye_work_.eye(nvar_, nvar_);
    oldbeta_work_.set_size(nvar_);
    betac_work_.set_size(nvar_);
}

void NegbinThetaUpdaterAdapter::update_unit(
    int i, mat const& rootpi_i, vec const& prior_mean_i)
{
    A_work_ = regdata_[i].hess + rootpi_i * trans(rootpi_i);
    U_work_ = sbeta_ * solve(A_work_, eye_work_);
    U_work_ = trans(chol(U_work_));
    oldbeta_work_ = vectorise(oldbetas_(i, span::all));
    betac_work_ = oldbeta_work_ + U_work_ * as<vec>(rnorm(nvar_));

    oldlpostbeta_[i] = lpostbeta(alpha_, oldbeta_work_,
                                  regdata_[i].X, regdata_[i].y, prior_mean_i, rootpi_i);
    clpostbeta_[i]   = lpostbeta(alpha_, betac_work_,
                                  regdata_[i].X, regdata_[i].y, prior_mean_i, rootpi_i);
    double ldiff = clpostbeta_[i] - oldlpostbeta_[i];
    double acc = exp(ldiff);
    if (acc > 1) acc = 1;
    double unif = (acc < 1) ? runif(1)[0] : 0;
    if (unif <= acc) {
        oldbetas_(i, span::all) = trans(betac_work_);
        nacceptbeta_++;
    }
}

void NegbinThetaUpdaterAdapter::draw_alpha() {
    if (fixalpha_) return;
    double logalphac = log(alpha_) + alphacroot_ * rnorm(1)[0];
    double oldlpostalpha = llnegbinpooled(regdata_, oldbetas_, alpha_) +
                           a_prior_ * log(alpha_) - b_prior_ * alpha_;
    double clpostalpha   = llnegbinpooled(regdata_, oldbetas_, exp(logalphac)) +
                           a_prior_ * logalphac - b_prior_ * exp(logalphac);
    double ldiff = clpostalpha - oldlpostalpha;
    double acc = exp(ldiff);
    if (acc > 1) acc = 1;
    double unif = (acc < 1) ? runif(1)[0] : 0;
    if (unif <= acc) {
        alpha_ = exp(logalphac);
        nacceptalpha_++;
    }
}

void NegbinThetaUpdaterAdapter::store(int mkeep) {
    alphadraw_[mkeep - 1] = alpha_;
}

UpdaterPackPayload NegbinThetaUpdaterAdapter::pack() {
    return {List::create(Named("alphadraw") = alphadraw_)};
}

AcceptanceMetrics NegbinThetaUpdaterAdapter::acceptance_metrics() const {
    AcceptanceMetrics am;
    am.acceptrbeta  = nacceptbeta_  / (R_ * nreg_ * 1.0) * 100;
    am.acceptralpha = nacceptalpha_ / (R_ * 1.0) * 100;
    return am;
}

double NegbinThetaUpdaterAdapter::current_loglike() const {
    return llnegbinpooled(regdata_, oldbetas_, alpha_);
}
