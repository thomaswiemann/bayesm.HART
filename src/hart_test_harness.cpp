// Test harnesses for the tree kernels.
#include <RcppArmadillo.h>
#include <sstream>
#include "vecbart.h"
#include "vcbart.h"
#include "direct_hart_backend.h"

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
// split counts and probabilities, and the pwvecbart() treedraws.
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

// Direct-model backend with pinnable blocks and readable state.
struct DirectHartProbe : DirectHartBackend {
    using DirectHartBackend::DirectHartBackend;
    void run(LogicalVector const& u) { sweep(u[0], u[1], u[2], u[3], u[4]); }
    vec const& mu() const { return mu_; }
    vec const& tau2() const { return tau2_; }
    mat const& f() const { return f_.allfit(); }
    mat const& v() const { return vmat_; }
    mat phi_at(int i) const {
        mat P(D_, D_, fill::zeros);
        if (full_cov_)
            for (int j = 1; j < D_; j++)
                for (int k = 0; k < j; k++) P(j, k) = phi_[j][k].f(i);
        return P;
    }
    List state() const {
        cube Q(D_, D_, n_);
        mat pm(n_, D_);
        for (int i = 0; i < n_; i++) {
            Q.slice(i) = Q_[i];
            pm.row(i)  = prior_mean_[i].t();
        }
        vec ss = f_.leaf_sumsq();
        return List::create(Named("mu")   = NumericVector(mu_.begin(), mu_.end()),
                            Named("tau2") = NumericVector(tau2_.begin(), tau2_.end()),
                            Named("f") = f_.allfit(), Named("v") = vmat_,
                            Named("phi") = phi_at(0), Named("Q") = Q,
                            Named("prior_mean") = pm,
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
    LogicalVector all(5, true);
    DirectHartProbe b(theta, Z, mubar, Amu, warm + nsweep, 1,
                      bart_params, var_params, phi_params);
    for (int r = 0; r < warm; r++) b.run(all);
    mat mu(nsweep, D), tau2(nsweep, D);
    cube f(n, D, nsweep), v(n, D, nsweep), phi(D, D, nsweep);
    for (int r = 0; r < nsweep; r++) {
        b.run(update);
        mu.row(r)    = b.mu().t();
        tau2.row(r)  = b.tau2().t();
        f.slice(r)   = b.f();
        v.slice(r)   = b.v();
        phi.slice(r) = b.phi_at(0);
    }
    return List::create(Named("mu") = mu, Named("tau2") = tau2,
                        Named("f") = f, Named("v") = v, Named("phi") = phi,
                        Named("state") = b.state());
}
