#include "logging.h"
#include <ctime>

std::map<std::string, std::ofstream> Logger::logFiles;
std::map<std::string, bool> Logger::initialized;

void Logger::ensureInitialized(const std::string& name) {
	if (!initialized[name]) {
		throw std::runtime_error("Logger '" + name + "' not initialized");
	}
}

void Logger::init(const std::string& name, const std::string& filename) {
	if (!initialized[name]) {
		logFiles[name].open(filename, std::ios::out | std::ios::app);
		initialized[name] = true;

		time_t now = time(0);
		logFiles[name] << "\n=== New Session Started at " << ctime(&now);
		logFiles[name] << "=== Logger: " << name << " ===\n\n";
	}
}

void Logger::log(const std::string& name, const std::string& message) {
	ensureInitialized(name);
	logFiles[name] << message << std::endl;
	logFiles[name].flush();
}

void Logger::logMatrix(const std::string& name, const std::string& label, const mat& matrix) {
	ensureInitialized(name);
	logFiles[name] << label << ":\n" << matrix;
	logFiles[name].flush();
}

void Logger::logVector(const std::string& name, const std::string& label, const vec& vector) {
	ensureInitialized(name);
	logFiles[name] << label << ":\n" << vector.t() << std::endl;
	logFiles[name].flush();
}

void Logger::logRowVector(const std::string& name, const std::string& label, const rowvec& vector) {
	ensureInitialized(name);
	logFiles[name] << label << ":\n" << vector << std::endl;
	logFiles[name].flush();
}

void Logger::logTrees(const std::string& name, const std::string& label, bart& model) {
	ensureInitialized(name);

	logFiles[name] << "\n=== " << label << " ===\n";

	size_t num_trees = model.getm();

	for (size_t j = 0; j < num_trees; j++) {
		logFiles[name] << "\nTree " << j << ":\n";

		// Get nodes for this tree
		tree::npv nodes;
		model.gettree(j).getnodes(nodes);

		// Log tree size first
		logFiles[name] << "tree size: " << model.gettree(j).treesize() << "\n";

		// Log node information
		for (size_t i = 0; i < nodes.size(); i++) {
			tree::tree_p node = nodes[i];

			// Get parent ID (0 if no parent)
			size_t pid = 0;
			if (node->getp()) pid = node->getp()->nid();

			// Calculate indent based on depth
			std::string pad(2 * node->depth(), ' ');

			logFiles[name] << pad
				<< "(id,parent): " << node->nid() << ", " << pid
				<< ", (v,c): " << node->getv() << ", " << node->getc()
				<< ", theta: " << node->gettheta()
				<< ", type: " << node->ntype()
				<< ", depth: " << node->depth()
				<< ", pointer: " << node << "\n";
		}
	}

	logFiles[name] << "\n----------------------------------------\n";
	logFiles[name].flush();
}

bool Logger::isInitialized(const std::string& name) {
	return initialized[name];
}

void Logger::close(const std::string& name) {
	if (initialized[name]) {
		time_t now = time(0);
		logFiles[name] << "\n=== Session Ended at " << ctime(&now);
		logFiles[name] << "=== Logger: " << name << " ===\n";
		logFiles[name].close();
		initialized[name] = false;
	}
}

void Logger::closeAll() {
	std::vector<std::string> toClose;
	for (const auto& kv : initialized) {
		if (kv.second) {
			toClose.push_back(kv.first);
		}
	}

	for (const auto& name : toClose) {
		close(name);
	}

	logFiles.clear();
	initialized.clear();
}