#ifndef GUARD_vectree_h
#define GUARD_vectree_h

#include <map>
#include <string>
#include <RcppArmadillo.h>
#include "tree.h"

// A regression tree with vector-valued leaves.
//
// The split structure is an ordinary `tree` (split variable, cutpoint,
// parent/child links, node ids); its scalar leaf value is unused.  Leaf
// vectors live in a map keyed by node id, so a vectree is a plain value type:
// copying copies structure and leaves, and nothing points into caller-owned
// memory.  birth() and death() keep the map aligned with the terminal nodes.
//
// Stream format (one tree):
//   nnodes
//   nid v c m_1 ... m_D          one line per node; zeros for interior nodes
class vectree {
public:
    typedef tree::tree_p  tree_p;
    typedef tree::tree_cp tree_cp;

    vectree() : D_(0) {}
    explicit vectree(size_t D) : D_(D) { leaf_[1] = arma::vec(D, arma::fill::zeros); }

    size_t dim() const { return D_; }
    tree& topo() { return t_; }
    tree const& topo() const { return t_; }

    arma::vec const& leaf(tree_cp n) const { return leaf_.at(n->nid()); }
    void setleaf(tree_cp n, arma::vec const& m) { leaf_[n->nid()] = m; }
    size_t nleaves() const { return leaf_.size(); }
    std::map<size_t, arma::vec> const& leaves() const { return leaf_; }

    // Split terminal node n on (v, c) with child leaves ml, mr.
    void birth(tree_p n, size_t v, size_t c,
               arma::vec const& ml, arma::vec const& mr);
    // Collapse the children of nog node n into the single leaf m.
    void death(tree_p n, arma::vec const& m);

    // Evaluate at the p x n column stack x: fv(i, .) = leaf containing x_i.
    void fit(xinfo& xi, size_t p, size_t n, double* x, arma::mat& fv);

    // Structure key "nid:v:c;..." over interior nodes.
    std::string key() const;

    friend std::ostream& operator<<(std::ostream&, vectree const&);
    friend std::istream& operator>>(std::istream&, vectree&);

protected:
    tree t_;
    size_t D_;
    std::map<size_t, arma::vec> leaf_;
};

std::ostream& operator<<(std::ostream&, vectree const&);
std::istream& operator>>(std::istream&, vectree&);

#endif
