#include "vectree.h"
#include <set>
#include <sstream>

void vectree::birth(tree_p n, size_t v, size_t c,
                  arma::vec const& ml, arma::vec const& mr) {
    size_t id = n->nid();
    t_.birthp(n, v, c, 0.0, 0.0);
    leaf_.erase(id);
    leaf_[2 * id]     = ml;
    leaf_[2 * id + 1] = mr;
}

void vectree::death(tree_p n, arma::vec const& m) {
    size_t id = n->nid();
    t_.deathp(n, 0.0);
    leaf_.erase(2 * id);
    leaf_.erase(2 * id + 1);
    leaf_[id] = m;
}

void vectree::fit(xinfo& xi, size_t p, size_t n, double* x, arma::mat& fv) {
    fv.set_size(n, D_);
    for (size_t i = 0; i < n; i++) {
        tree_p bn = t_.bn(x + i * p, xi);
        fv.row(i) = leaf_.at(bn->nid()).t();
    }
}

std::string vectree::key() const {
    tree::cnpv nds;
    t_.getnodes(nds);
    std::set<size_t> ids;
    for (tree_cp n : nds) ids.insert(n->nid());
    std::map<size_t, std::string> parts;   // ordered by node id
    for (tree_cp n : nds) {
        size_t id = n->nid();
        if (!ids.count(2 * id)) continue;    // terminal
        std::ostringstream os;
        os << id << ":" << n->getv() << ":" << n->getc() << ";";
        parts[id] = os.str();
    }
    std::string key;
    for (auto const& kv : parts) key += kv.second;
    return key;
}

std::ostream& operator<<(std::ostream& os, vectree const& t) {
    tree::cnpv nds;
    t.t_.getnodes(nds);
    os << nds.size() << std::endl;
    for (vectree::tree_cp n : nds) {
        os << n->nid() << " " << n->getv() << " " << n->getc();
        auto it = t.leaf_.find(n->nid());
        for (size_t d = 0; d < t.D_; d++)
            os << " " << (it == t.leaf_.end() ? 0.0 : it->second(d));
        os << std::endl;
    }
    return os;
}

std::istream& operator>>(std::istream& is, vectree& t) {
    size_t nn;
    is >> nn;
    if (!is) return is;
    struct rec { size_t v, c; arma::vec m; };
    std::map<size_t, rec> recs;                // ordered by node id
    for (size_t i = 0; i < nn; i++) {
        size_t id, v, c;
        is >> id >> v >> c;
        arma::vec m(t.D_);
        for (size_t d = 0; d < t.D_; d++) is >> m(d);
        if (!is) return is;
        recs[id] = rec{v, c, m};
    }
    // Rebuild through the public tree API: parents precede children in node
    // id order, so every split target already exists when it is split.
    t.t_.tonull();
    t.leaf_.clear();
    for (auto const& kv : recs) {
        size_t id = kv.first;
        if (recs.count(2 * id)) t.t_.birth(id, kv.second.v, kv.second.c, 0.0, 0.0);
        else                    t.leaf_[id] = kv.second.m;
    }
    return is;
}
