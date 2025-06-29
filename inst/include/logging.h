#ifndef __LOGGING_H__
#define __LOGGING_H__

#include <RcppArmadillo.h>
#include <fstream>
#include <string>
#include <map>
#include "bart.h"

using namespace arma;

class Logger {
private:
	static std::map<std::string, std::ofstream> logFiles;
	static std::map<std::string, bool> initialized;
	static void ensureInitialized(const std::string& name);

public:
	static void init(const std::string& name, const std::string& filename);
	static void log(const std::string& name, const std::string& message);

	template<typename T>
	static void log(const std::string& name, const std::string& label, const T& value) {
		if (initialized[name]) {
			logFiles[name] << label << ": " << value << std::endl;
			logFiles[name].flush();
		}
	}

	template<typename T>
	static void logAll(const std::string& label, const T& value) {
		for (const auto& kv : logFiles) {
			if (initialized[kv.first]) {
				log(kv.first, label, value);
			}
		}
	}

	static void logAll(const std::string& message) {
		for (const auto& kv : logFiles) {
			if (initialized[kv.first]) {
				log(kv.first, message);
			}
		}
	}

	static void logMatrix(const std::string& name, const std::string& label, const mat& matrix);
	static void logVector(const std::string& name, const std::string& label, const vec& vector);
	static void logRowVector(const std::string& name, const std::string& label, const rowvec& vector);
	static void logTrees(const std::string& name, const std::string& label, bart& model);
	static bool isInitialized(const std::string& name);
	static void close(const std::string& name);
	static void closeAll();
};

#endif