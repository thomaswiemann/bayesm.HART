#include "bayesm.HART.h"
#include "unified_sampler_driver.h"
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <string>

void run_unified_sampler(
    int R, int keep, int nprint, int nreg,
    CovarianceBackend& backend,
    ThetaUpdater& updater,
    arma::mat& oldbetas,
    arma::cube& betadraw,
    arma::vec& loglike) {

    using clock_t = std::chrono::steady_clock;
    using ns_t = std::chrono::nanoseconds;

    const char* profile_env = std::getenv("BAYESM_HART_PROFILE");
    bool profile_enabled = false;
    if (profile_env != nullptr) {
        std::string profile_str(profile_env);
        profile_enabled = !profile_str.empty() && profile_str != "0";
    }

    std::uint64_t backend_ns_total = 0;
    std::uint64_t theta_loop_ns_total = 0;
    std::uint64_t post_iter_ns_total = 0;
    std::uint64_t store_ns_total = 0;

    if (nprint > 0) startMcmcTimer();

    for (int rep = 0; rep < R; rep++) {

        // Step A-D: backend draws covariance/mean-tree blocks (atomic).
        auto t_backend_start = profile_enabled ? clock_t::now() : clock_t::time_point{};
        backend.draw_iteration(rep);
        if (profile_enabled) {
            backend_ns_total += static_cast<std::uint64_t>(
                std::chrono::duration_cast<ns_t>(clock_t::now() - t_backend_start).count());
        }

        // Step E: per-unit theta update using backend-provided prior kernel.
        auto t_theta_start = profile_enabled ? clock_t::now() : clock_t::time_point{};
        for (int i = 0; i < nreg; i++) {
            arma::mat const& rootpi_i    = backend.precision_root(i);
            arma::vec const& prior_mean_i = backend.prior_mean(i);
            updater.update_unit(i, rootpi_i, prior_mean_i);
        }
        if (profile_enabled) {
            theta_loop_ns_total += static_cast<std::uint64_t>(
                std::chrono::duration_cast<ns_t>(clock_t::now() - t_theta_start).count());
        }

        // Per-iteration hook (e.g. NegBin alpha MH). No-op by default.
        auto t_post_start = profile_enabled ? clock_t::now() : clock_t::time_point{};
        updater.post_iteration();
        if (profile_enabled) {
            post_iter_ns_total += static_cast<std::uint64_t>(
                std::chrono::duration_cast<ns_t>(clock_t::now() - t_post_start).count());
        }

        if (nprint > 0 && (rep + 1) % nprint == 0) infoMcmcTimer(rep, R);

        if ((rep + 1) % keep == 0) {
            auto t_store_start = profile_enabled ? clock_t::now() : clock_t::time_point{};
            int mkeep = (rep + 1) / keep;
            betadraw.slice(mkeep - 1) = oldbetas;
            loglike[mkeep - 1] = updater.current_loglike();
            backend.store(mkeep);
            updater.store(mkeep);
            if (profile_enabled) {
                store_ns_total += static_cast<std::uint64_t>(
                    std::chrono::duration_cast<ns_t>(clock_t::now() - t_store_start).count());
            }
        }
    }

    if (nprint > 0) endMcmcTimer();

    if (profile_enabled && R > 0) {
        const double ns_to_ms = 1.0e-6;
        double backend_ms = backend_ns_total * ns_to_ms;
        double theta_ms = theta_loop_ns_total * ns_to_ms;
        double post_ms = post_iter_ns_total * ns_to_ms;
        double store_ms = store_ns_total * ns_to_ms;
        double measured_ms = backend_ms + theta_ms + post_ms + store_ms;

        Rcpp::Rcout << "[profile] unified_sampler draws=" << R << "\n";
        Rcpp::Rcout << "[profile] backend_ms_total=" << backend_ms
                    << " backend_ms_per_draw=" << (backend_ms / R) << "\n";
        Rcpp::Rcout << "[profile] theta_loop_ms_total=" << theta_ms
                    << " theta_loop_ms_per_draw=" << (theta_ms / R) << "\n";
        Rcpp::Rcout << "[profile] post_iter_ms_total=" << post_ms
                    << " post_iter_ms_per_draw=" << (post_ms / R) << "\n";
        Rcpp::Rcout << "[profile] store_ms_total=" << store_ms
                    << " store_ms_per_draw=" << (store_ms / R) << "\n";
        if (measured_ms > 0.0) {
            Rcpp::Rcout << "[profile] shares backend=" << (backend_ms / measured_ms)
                        << " theta=" << (theta_ms / measured_ms)
                        << " post=" << (post_ms / measured_ms)
                        << " store=" << (store_ms / measured_ms) << "\n";
        }
    }
}
