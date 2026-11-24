#include "unique_z_cache.h"

using namespace Rcpp;
using namespace arma;

void UniqueZCache::init(List const& params, int n, int nvar, int nkeep,
                        bool with_sigma, const char* who) {
    active = false;
    if (!params.containsElementNamed("z_index")) return;

    IntegerVector z_index = params["z_index"];
    if ((int)z_index.size() != n)
        stop("%s: z_index length must match the number of units.", who);
    int nu = params.containsElementNamed("n_unique")
        ? as<int>(params["n_unique"]) : 0;
    if (nu <= 0) nu = Rcpp::max(z_index);
    if (nu <= 0)
        stop("%s: n_unique must be positive when z_index is supplied.", who);

    std::vector<int> first(nu, -1);
    for (int i = 0; i < n; i++) {
        int u = z_index[i] - 1;  // 1-based from R
        if (u < 0 || u >= nu)
            stop("%s: z_index values must be in 1..n_unique.", who);
        if (first[u] < 0) first[u] = i;
    }
    first_idx.set_size(nu);
    for (int u = 0; u < nu; u++) {
        if (first[u] < 0)
            stop("%s: every unique index must appear at least once.", who);
        first_idx[u] = static_cast<uword>(first[u]);
    }

    delta_draws.zeros(nu, nvar, nkeep);
    if (with_sigma) sigma_draws_flat.zeros(nu, nvar * nvar, nkeep);
    active = true;
}

void UniqueZCache::store_delta(int mkeep, mat const& delta_Z) {
    if (!active) return;
    delta_draws.slice(mkeep - 1) = delta_Z.rows(first_idx);
}

void UniqueZCache::store_sigma(int mkeep, uword u, mat const& Sigma) {
    if (!active || sigma_draws_flat.n_elem == 0) return;
    sigma_draws_flat.slice(mkeep - 1).row(u) = trans(vectorise(Sigma));
}

void UniqueZCache::pack(List& out) const {
    if (!active) return;
    out["DeltaZ_unique_draws"] = delta_draws;
    if (sigma_draws_flat.n_elem > 0)
        out["SigmaZ_unique_draws_flat"] = sigma_draws_flat;
}
