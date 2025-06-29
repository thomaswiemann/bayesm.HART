#include "logging_bart.h"
#include "logging.h"
#include "bart.h"

namespace bart_logging {

	void init_loggers() {
		Logger::init("delta", "debug_delta.log");
	}

	void log_draw_state(const bart& model, int tree_index, const char* stage,
		int debug_obs, double value) {
		Logger::log("delta", std::string("\n") + stage + " tree " +
			std::to_string(tree_index) + " contribution:");
		Logger::log("delta", "  value[" + std::to_string(debug_obs) +
			"] = " + std::to_string(value));
	}

	void log_tree_prediction(const bart& model, int tree_index, int obs_index,
		double theta, int node_id, double fpred) {
		Logger::log("delta", "Tree " + std::to_string(tree_index) +
			" prediction details for obs " + std::to_string(obs_index) + ":");
		Logger::log("delta", "  Node ID: " + std::to_string(node_id));
		Logger::log("delta", "  Node theta: " + std::to_string(theta));
		Logger::log("delta", "  Predicted value: " + std::to_string(fpred));
	}

	void log_prediction_details(const bart& model, int obs_index,
		const std::vector<double>& tree_preds,
		double total_pred) {
		Logger::log("delta", "\nPrediction breakdown for obs " +
			std::to_string(obs_index) + ":");
		for (size_t i = 0; i < tree_preds.size(); i++) {
			Logger::log("delta", "  Tree " + std::to_string(i) +
				" contribution: " + std::to_string(tree_preds[i]));
		}
		Logger::log("delta", "  Total prediction: " + std::to_string(total_pred));
	}

	void log_tree_structure(const bart& model, int tree_index, const char* stage,
		const std::vector<tree*>& nodes) {
		Logger::log("delta", std::string("\n") + stage + " for tree " +
			std::to_string(tree_index) + ":");

		Logger::log("delta", "tree size: " + std::to_string(nodes.size()));

		for (size_t i = 0; i < nodes.size(); i++) {
			tree::tree_p node = nodes[i];
			size_t pid = 0;
			if (node->getp()) pid = node->getp()->nid();
			std::string pad(2 * node->depth(), ' ');
			std::string node_info = pad + "(id,parent): " + std::to_string(node->nid()) +
				", " + std::to_string(pid) +
				", (v,c): " + std::to_string(node->getv()) + ", " +
				std::to_string(node->getc()) +
				", theta: " + std::to_string(node->gettheta()) +
				", type: " + std::string(1, node->ntype()) +
				", depth: " + std::to_string(node->depth());
			Logger::log("delta", node_info);
		}

		Logger::log("delta", ""); // Empty line for readability
	}

	void log_node_details(const bart& model, int tree_index, int node_id, double theta) {
		Logger::log("delta", "Node details for tree " + std::to_string(tree_index) + ":");
		Logger::log("delta", "  Node ID: " + std::to_string(node_id));
		Logger::log("delta", "  Node theta: " + std::to_string(theta));
	}

	void log_tree_traversal(const bart& model, int tree_index, int obs_index,
		double* x, const xinfo& xi, tree_p node) {
		Logger::log("delta", "\nTraversing tree " + std::to_string(tree_index) +
			" for observation " + std::to_string(obs_index));

		// We can walk up from the node to reconstruct the path
		std::vector<tree_p> path;
		tree_p current = node;
		while (current) {
			path.push_back(current);
			current = current->getp();
		}

		// Now print the path in reverse (root to leaf)
		for (int i = path.size() - 1; i >= 0; i--) {
			tree_p n = path[i];
			if (n->getp()) { // Skip root node comparisons
				std::string msg = "  Node " + std::to_string(n->nid()) +
					": var " + std::to_string(n->getv()) +
					" = " + std::to_string(x[n->getv() + obs_index * model.getp()]) +
					" vs cut[" + std::to_string(n->getc()) +
					"] = " + std::to_string(xi[n->getv()][n->getc()]);
				Logger::log("delta", msg);
			}
		}

		Logger::log("delta", "  Landed at node " + std::to_string(node->nid()) +
			" with theta = " + std::to_string(node->gettheta()));
	}

	void log_prediction_mismatch(const bart& model, int tree_index, int obs_index,
		double manual_pred, double auto_pred, const xinfo& xi) {
		Logger::log("delta", "\nPrediction mismatch detected in tree " +
			std::to_string(tree_index) +
			" for observation " + std::to_string(obs_index));
		Logger::log("delta", "  Manual prediction: " + std::to_string(manual_pred));
		Logger::log("delta", "  Auto prediction: " + std::to_string(auto_pred));
	}

}