// Test harnesses for the tree kernels.
#include <RcppArmadillo.h>
#include <sstream>
#include "vecbart.h"
#include "vcbart.h"
#include "direct_hart_backend.h"
#include <memory>
#include "tree_cov.h"
#include "constant_cov.h"
#include "bart_pack.h"

using namespace Rcpp;
using namespace arma;

namespace {

// Structure key of a scalar tree, same format as vectree::key().
std::string tree_key(tree& t) {
    tree::npv nds;
    t.getnodes(nds);
    std::map<size_t, std::string> parts;
    for (tree::tree_p nd : nds) {
        if (!nd->getl()) continue;
        std::ostringstream os;
        os << nd->nid() << ":" << nd->getv() << ":" << nd->getc() << ";";
        parts[nd->nid()] = os.str();
    }
    std::string k;
    for (auto const& kv : parts) k += kv.second;
    return k;
}

std::vector<mat> slices(cube const& Q) {
    std::vector<mat> v(Q.n_slices);
    for (uword i = 0; i < Q.n_slices; i++) v[i] = Q.slice(i);
    return v;
}

List cutpoints(xinfo& xi) {
    List out(xi.size());
    for (size_t j = 0; j < xi.size(); j++)
        out[j] = NumericVector(xi[j].begin(), xi[j].end());
    return out;
}

}  // namespace

// Leaf marginal of a vector leaf with statistics (n, P, b) under N(0, C).
// [[Rcpp::export]]
double veclh_R(arma::mat const& P, arma::vec const& b, int n, arma::mat const& C) {
    vecsuff s(b.n_elem);
    s.n = n; s.P = P; s.b = b;
    return veclh(s, inv_sympd(C), log_det_sympd(C));
}

// Leaf conditional draws: ndraw x D.
// [[Rcpp::export]]
arma::mat vecdrawnodemu_R(arma::mat const& P, arma::vec const& b,
                          arma::mat const& C, int ndraw) {
    arn gen;
    vecsuff s(b.n_elem);
    s.n = 1; s.P = P; s.b = b;
    mat Cinv = inv_sympd(C);
    mat out(ndraw, b.n_elem);
    for (int k = 0; k < ndraw; k++) out.row(k) = vecdrawnodemu(s, Cinv, gen).t();
    return out;
}

// niter backfitting sweeps of a vecbart with m trees on Zt (p x n), y (n x D),
// Q (D x D x n): fits after every sweep, structure key of tree 1, leaf counts,
// split counts and probabilities, and the serialized treedraws.
// [[Rcpp::export]]
Rcpp::List vecbart_chain_R(arma::mat const& Zt, arma::mat const& y,
                           arma::cube const& Q, int m, arma::mat const& C,
                           double alpha, double beta, int nmin, int numcut,
                           int niter, bool sparse = false, int burn = 0) {
    size_t p = Zt.n_rows, n = Zt.n_cols, D = y.n_cols;
    if (y.n_rows != n || Q.n_slices != n || Q.n_rows != D || Q.n_cols != D)
        stop("vecbart_chain_R: y must be n x D and Q must be D x D x n.");
    std::vector<mat> Qv = slices(Q);
    mat Zc = Zt;                          // the ensemble keeps a pointer into it
    arn gen;

    vecbart f(m, D);
    f.setprior(alpha, beta);
    f.setleafguard(nmin);
    f.setdata(p, n, Zc.memptr(), numcut);
    f.setleafcov(C);
    f.setdart(0.5, 1.0, (double)p, false, sparse, 0.0, 1.0);

    std::stringstream tss;
    tss.precision(10);
    tss << niter << " " << m << " " << p << " " << D << std::endl;
    cube fits(n, D, niter);
    CharacterVector keys(niter);
    IntegerVector nleaves(niter);
    for (int it = 0; it < niter; it++) {
        if (sparse && it == burn + 1) f.startdart();
        f.draw(y, Qv, gen);
        fits.slice(it) = f.allfit();
        keys[it] = f.gettree(0).key();
        nleaves[it] = (int)f.nleaves();
        for (int j = 0; j < m; j++) tss << f.gettree(j);
    }
    IntegerVector nv(p);
    NumericVector pv(p);
    for (size_t j = 0; j < p; j++) {
        nv[j] = (int)f.getnv()[j]; pv[j] = f.getpv()[j];
    }
    List treedraws = List::create(Named("cutpoints") = cutpoints(f.getxinfo()),
                                  Named("trees")     = CharacterVector(tss.str()),
                                  Named("leaf_dim")  = (int)D);
    return List::create(Named("fits")       = fits,
                        Named("keys")       = keys,
                        Named("nleaves")    = nleaves,
                        Named("leaf_sumsq") = f.leaf_sumsq(),
                        Named("nv")         = nv,
                        Named("pv")         = pv,
                        Named("treedraws")  = treedraws);
}

// Varying-coefficient leaf score at statistics (P, b) under N(0, tau^2).
// [[Rcpp::export]]
double vclh_R(double P, double b, double tau) { return vclh(P, b, tau); }

// niter backfitting sweeps of a vcbart with m trees, y_i = a_i f(Z_i) + e_i,
// e_i ~ N(0, v_i): f(Z_i) after every sweep (n x niter), structure key of
// tree 1, leaf counts, split counts and probabilities, and the pwbart()
// treedraws.
// [[Rcpp::export]]
Rcpp::List vcbart_chain_R(arma::mat const& Zt, arma::vec const& y,
                          arma::vec const& a, arma::vec const& v, int m,
                          double tau, double alpha, double beta, int nmin,
                          int numcut, int niter, bool sparse = false,
                          int burn = 0) {
    size_t p = Zt.n_rows, n = Zt.n_cols;
    if (y.n_elem != n || a.n_elem != n || v.n_elem != n)
        stop("vcbart_chain_R: y, a, v must have one entry per column of Zt.");
    mat Zc = Zt;                          // the ensemble keeps pointers into these
    vec yc = y, ac = a, vc = v;
    arn gen;

    vcbart f(m);
    f.setprior(alpha, beta, tau);
    f.setleafguard(nmin);
    f.setdata(p, n, Zc.memptr(), yc.memptr(), numcut);
    f.setdart(0.5, 1.0, (double)p, false, sparse, 0.0, 1.0);

    std::stringstream tss;
    tss.precision(10);
    tss << niter << " " << m << " " << p << std::endl;
    mat fits(n, niter);
    CharacterVector keys(niter);
    IntegerVector nleaves(niter);
    for (int it = 0; it < niter; it++) {
        if (sparse && it == burn + 1) f.startdart();
        f.draw(ac.memptr(), vc.memptr(), gen);
        for (size_t i = 0; i < n; i++) fits(i, it) = f.f(i);
        keys[it] = tree_key(f.gettree(0));
        size_t L = 0;
        for (int j = 0; j < m; j++) L += f.gettree(j).nbots();
        nleaves[it] = (int)L;
        for (int j = 0; j < m; j++) tss << f.gettree(j);
    }
    IntegerVector nv(p);
    NumericVector pv(p);
    for (size_t j = 0; j < p; j++) {
        nv[j] = (int)f.getnv()[j]; pv[j] = f.getpv()[j];
    }
    List treedraws = List::create(Named("cutpoints") = cutpoints(f.getxinfo()),
                                  Named("trees")     = CharacterVector(tss.str()));
    return List::create(Named("fits")      = fits,
                        Named("keys")      = keys,
                        Named("nleaves")   = nleaves,
                        Named("nv")        = nv,
                        Named("pv")        = pv,
                        Named("treedraws") = treedraws);
}

// (L, D) at prescribed d (length D) and phi (D x D, strictly lower part):
// rootpi and L v through stump ensembles.
// [[Rcpp::export]]
Rcpp::List cov_eval_R(arma::vec const& d, arma::mat const& phi,
                      arma::vec const& v, bool diagonal) {
    size_t D = d.n_elem;
    mat Zt(1, 1, fill::zeros);
    vec y(1, fill::zeros);
    std::vector<varbart> dm(D);
    std::vector<std::vector<vcbart>> pm;
    for (size_t j = 0; j < D; j++) {
        dm[j] = varbart(1);
        dm[j].setdata(1, 1, Zt.memptr(), y.memptr(), 2);
        dm[j].setvarprior(3.0, 1.0, d(j));
    }
    if (!diagonal) {
        pm.resize(D);
        for (size_t j = 1; j < D; j++) {
            pm[j] = std::vector<vcbart>(j);
            for (size_t k = 0; k < j; k++) {
                pm[j][k] = vcbart(1);
                pm[j][k].gettree(0).settheta(phi(j, k));
                pm[j][k].setdata(1, 1, Zt.memptr(), y.memptr(), 2);
            }
        }
    }
    cov::cov_eval ev{dm, pm, 0};
    vec Lv = v;
    cov::apply_L(Lv, ev);
    return List::create(Named("rootpi") = cov::rootpi(ev), Named("Lv") = Lv);
}

namespace {

// Covariance strategies with pinnable blocks and readable state.
struct TreeCovProbe : TreeCov {
    using TreeCov::TreeCov;
    bool do_v = true, do_phi = true;
    void draw(mat const& E, rn& gen) override {
        refresh_residuals(E);
        if (do_v) draw_variance_trees(gen);
        if (do_phi) draw_phi_trees(gen);
        refresh_precision();
    }
    mat const& v() const { return vmat_; }
    std::vector<mat> const& Q() const { return Q_; }
    mat phi_at(int i) const {
        mat P(D_, D_, fill::zeros);
        if (full_cov_)
            for (int j = 1; j < D_; j++)
                for (int k = 0; k < j; k++) P(j, k) = phi_[j][k].f(i);
        return P;
    }
};

struct ConstantCovProbe : ConstantCov {
    using ConstantCov::ConstantCov;
    bool do_sigma = true;
    void draw(mat const& E, rn& gen) override {
        if (do_sigma) ConstantCov::draw(E, gen);
        else set_precision();
    }
    mat const& Sigma() const { return Sigma_; }
    mat const& rooti() const { return rooti_; }
    std::vector<mat> const& Q() const { return Q_; }
};

// Direct-model backend with pinnable mean blocks and readable state.
struct DirectHartProbe : DirectHartBackend {
    using DirectHartBackend::DirectHartBackend;
    void run(bool mu, bool f, bool tau) { sweep(mu, true, f, tau); }
    vec const& mu() const { return mu_; }
    vec const& tau2() const { return tau2_; }
    mat const& f() const { return f_.allfit(); }
    List state(std::vector<mat> const& Q) const {
        cube Qc(D_, D_, n_);
        mat pm(n_, D_);
        for (int i = 0; i < n_; i++) {
            Qc.slice(i) = Q[i];
            pm.row(i)   = prior_mean_[i].t();
        }
        vec ss = f_.leaf_sumsq();
        return List::create(Named("mu")   = NumericVector(mu_.begin(), mu_.end()),
                            Named("tau2") = NumericVector(tau2_.begin(), tau2_.end()),
                            Named("f") = f_.allfit(),
                            Named("prior_mean") = pm, Named("Q") = Qc,
                            Named("nleaves") = (int)f_.nleaves(),
                            Named("leaf_sumsq") = NumericVector(ss.begin(), ss.end()));
    }
};

}  // namespace

// Direct-model blocks on fixed theta: `warm` full sweeps, then `nsweep`
// sweeps with the blocks (mu, v, phi, f, tau) switched by `update`.  Returns
// per-sweep mu and tau^2, f and v at every unit, phi at unit 1, and the final
// state.
// [[Rcpp::export]]
Rcpp::List direct_hart_R(arma::mat const& theta, arma::mat const& Z,
                         arma::vec const& mubar, arma::mat const& Amu,
                         Rcpp::List const& bart_params,
                         Rcpp::List const& var_params,
                         Rcpp::List const& phi_params,
                         int warm, int nsweep,
                         Rcpp::LogicalVector const& update) {
    if (update.size() != 5) stop("direct_hart_R: update needs 5 entries.");
    int n = theta.n_rows, D = theta.n_cols;
    auto owned = std::make_unique<TreeCovProbe>(var_params, phi_params);
    TreeCovProbe* cov = owned.get();
    DirectHartProbe b(theta, Z, mubar, Amu, warm + nsweep, 1, bart_params,
                      std::move(owned));
    for (int r = 0; r < warm; r++) b.run(true, true, true);
    cov->do_v = update[1]; cov->do_phi = update[2];
    mat mu(nsweep, D), tau2(nsweep, D);
    cube f(n, D, nsweep), v(n, D, nsweep), phi(D, D, nsweep);
    for (int r = 0; r < nsweep; r++) {
        b.run(update[0], update[3], update[4]);
        mu.row(r)    = b.mu().t();
        tau2.row(r)  = b.tau2().t();
        f.slice(r)   = b.f();
        v.slice(r)   = cov->v();
        phi.slice(r) = cov->phi_at(0);
    }
    List st = b.state(cov->Q());
    st["v"] = cov->v();
    st["phi"] = cov->phi_at(0);
    return List::create(Named("mu") = mu, Named("tau2") = tau2,
                        Named("f") = f, Named("v") = v, Named("phi") = phi,
                        Named("state") = st);
}

// Constant-covariance direct blocks on fixed theta: `warm` full sweeps, then
// `nsweep` sweeps with the blocks (mu, Sigma, f, tau) switched by `update`.
// Returns per-sweep mu, tau^2, Sigma and f, and the final state.
// [[Rcpp::export]]
Rcpp::List direct_hart_constant_R(arma::mat const& theta, arma::mat const& Z,
                                  arma::vec const& mubar, arma::mat const& Amu,
                                  double nu, arma::mat const& V,
                                  Rcpp::List const& bart_params,
                                  int warm, int nsweep,
                                  Rcpp::LogicalVector const& update) {
    if (update.size() != 4) stop("direct_hart_constant_R: update needs 4 entries.");
    int n = theta.n_rows, D = theta.n_cols;
    auto owned = std::make_unique<ConstantCovProbe>(nu, V);
    ConstantCovProbe* cov = owned.get();
    DirectHartProbe b(theta, Z, mubar, Amu, warm + nsweep, 1, bart_params,
                      std::move(owned));
    for (int r = 0; r < warm; r++) b.run(true, true, true);
    cov->do_sigma = update[1];
    mat mu(nsweep, D), tau2(nsweep, D);
    cube f(n, D, nsweep), Sigma(D, D, nsweep);
    for (int r = 0; r < nsweep; r++) {
        b.run(update[0], update[2], update[3]);
        mu.row(r)      = b.mu().t();
        tau2.row(r)    = b.tau2().t();
        f.slice(r)     = b.f();
        Sigma.slice(r) = cov->Sigma();
    }
    List st = b.state(cov->Q());
    st["Sigma"] = cov->Sigma();
    st["rooti"] = cov->rooti();
    return List::create(Named("mu") = mu, Named("tau2") = tau2, Named("f") = f,
                        Named("Sigma") = Sigma, Named("state") = st);
}

namespace {

void state_prior(bart&) {}
void state_prior(varbart& f) { f.setvarprior(6.0, 1.0); }
void state_draw(bart& f, size_t n, arn& gen) { f.draw(0.3, gen); }
void state_draw(varbart& f, size_t n, arn& gen) { f.draw(gen); }
void state_draw(vcbart& f, size_t n, arn& gen) {
    vec a(n, fill::ones), v(n, fill::value(0.09));
    f.draw(a.memptr(), v.memptr(), gen);
}

template <class B>
List scalar_state(B& f, mat& Zt, vec& y, int niter, bool active,
                  Nullable<List> state) {
    size_t p = Zt.n_rows, n = Zt.n_cols;
    f.setprior(0.95, 2.0, 1.0);
    f.setleafguard(2);
    f.setdata(p, n, Zt.memptr(), y.memptr(), 20);
    state_prior(f);
    f.setdart(0.5, 1.0, (double)p, false, true, 0.0, 1.0);
    if (state.isNotNull()) restore_bart_state_(f, List(state));
    else if (active) f.startdart();
    arn gen;
    for (int it = 0; it < niter; it++) state_draw(f, n, gen);
    mat fitted(n, 1);
    for (size_t i = 0; i < n; i++) fitted(i, 0) = f.f(i);
    size_t leaves = 0;
    for (size_t j = 0; j < f.getm(); j++) leaves += f.gettree(j).nbots();
    return List::create(Named("snapshot") = pack_bart_state_(f),
                        Named("fit") = fitted, Named("nleaves") = (double)leaves);
}

}  // namespace

// Terminal state and fitted values on a freshly allocated ensemble.
// [[Rcpp::export]]
Rcpp::List ensemble_state_R(std::string kind, arma::mat Zt, arma::mat y,
                            int m, int niter, bool active,
                            Rcpp::Nullable<Rcpp::List> state = R_NilValue) {
    size_t p = Zt.n_rows, n = Zt.n_cols;
    if (p == 0 || n == 0 || y.n_rows != n || y.n_cols == 0 || m < 1 || niter < 0)
        stop("ensemble_state_R: invalid dimensions or iteration count");
    if (kind == "vecbart") {
        vecbart f(m, y.n_cols);
        f.setprior(0.95, 2.0);
        f.setleafguard(2);
        f.setdata(p, n, Zt.memptr(), 20);
        f.setleafcov(eye<mat>(y.n_cols, y.n_cols));
        f.setdart(0.5, 1.0, (double)p, false, true, 0.0, 1.0);
        if (state.isNotNull()) restore_bart_state_(f, List(state));
        else if (active) f.startdart();
        std::vector<mat> Q(n, eye<mat>(y.n_cols, y.n_cols) / 0.09);
        arn gen;
        for (int it = 0; it < niter; it++) f.draw(y, Q, gen);
        return List::create(Named("snapshot") = pack_bart_state_(f),
                            Named("fit") = f.allfit(),
                            Named("nleaves") = (double)f.nleaves(),
                            Named("leaf_sumsq") = f.leaf_sumsq());
    }
    if (y.n_cols != 1) stop("ensemble_state_R: scalar response required");
    vec yc = y.col(0);
    if (kind == "bart") {
        bart f(m);
        return scalar_state(f, Zt, yc, niter, active, state);
    }
    if (kind == "varbart") {
        varbart f(m);
        return scalar_state(f, Zt, yc, niter, active, state);
    }
    if (kind == "vcbart") {
        vcbart f(m);
        return scalar_state(f, Zt, yc, niter, active, state);
    }
    stop("ensemble_state_R: unknown ensemble");
}
