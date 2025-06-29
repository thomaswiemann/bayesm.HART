#ifndef GUARD_bart_logging_h
#define GUARD_bart_logging_h

#include <vector>
#include <string>

// Forward declarations  
class bart;
class tree;
typedef tree* tree_p;
typedef std::vector<double> xinfo_t; 
typedef std::vector<xinfo_t> xinfo;

namespace bart_logging {
	void log_draw_state(const bart& model, int tree_index, const char* stage,
		int debug_obs, double value);

	void log_tree_prediction(const bart& model, int tree_index, int obs_index,
		double theta, int node_id, double fpred);

	void log_prediction_details(const bart& model, int obs_index,
		const std::vector<double>& tree_preds,
		double total_pred);

	void log_tree_structure(const bart& model, int tree_index, const char* stage,
		const std::vector<tree*>& nodes);

	void log_node_details(const bart& model, int tree_index, int node_id, double theta);

	void log_tree_traversal(const bart& model, int tree_index, int obs_index,
		double* x, const xinfo& xi, tree_p node);

	void log_prediction_mismatch(const bart& model, int tree_index, int obs_index,
		double manual_pred, double auto_pred, const xinfo& xi);
}

#endif