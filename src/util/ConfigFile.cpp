/*
 * ConfigFile.cpp
 *
 *  Created on: Jan 21, 2016
 *      Author: gjones
 */

#include "ConfigFile.h"

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <sys/types.h>
#include <sys/stat.h>

#include "Reporter.h"
#include "Util.h"

namespace GarethUtil {

using namespace std;
using namespace boost::filesystem;

/**
 * Return true if str is a directory path
 *
 * @param str
 * @return
 */
static bool existsAndIsDirectory(const string & str) {
	path path(str);

	return exists(path) && is_directory(path);
}

const string ConfigFile::configDirectory() {
	boost::optional<string> configDirEnv = getEnv("CONFIG_DIR");
	string configDir;
	if (configDirEnv) {
		configDir = configDirEnv.get();
		if (existsAndIsDirectory(configDir)) {
			return configDir;
		} else {
			auto message =
					(boost::format(
							"CONFIG_DIR set to %s, but that path does not exist or is not a directory")
							% configDir).str();
			REPORT(Reporter::FATAL) << message;
			throw runtime_error(message);
		}
	}

	for (const auto & path : { "../config", "/home/gjones/src/gape/config" }) {
		if (existsAndIsDirectory(path)) {
			REPORT(Reporter::NORMAL) << "Using " << path
					<< " as configuration directory";
			return path;
		}
	}

	auto message =
			"Unable to find configuration directory (set using environment variable CONFIG_DIR)";
	REPORT(Reporter::FATAL) << message;
	throw runtime_error(message);

}

string ConfigFile::configFilePath(const string &configFile) {
	path configDir = path(configDirectory());
	path file = path(configFile);
	file = configDir / file;
	if (!exists(file) || !is_regular_file(file)) {
		auto message =
				(boost::format(
						"Configuration file %s does not exist or is not a regualar file")
						% file).str();
		REPORT(Reporter::FATAL) << message;
		throw runtime_error(message);
	}
	return file.string();
}

ConfigFile::ConfigFile(const string & fileName_) :
		fileName(fileName_), fullPath(configFilePath(fileName_)) {
	in.open(fullPath);
	if (!in.is_open()) {
		auto message = (boost::format("Failed to open configuration file %s")
				% fullPath).str();
		REPORT(Reporter::FATAL) << message;
		throw runtime_error(message);
	}
}

bool ConfigFile::readConfigLine() {
	inLine = false;
	while (getline(in, currentLine)) {
		currentLine = trim(currentLine);
		if (currentLine.empty()) {
			continue;
		}
		if (startsWith(currentLine, "#")) {
			continue;
		}
		ss = istringstream(currentLine);
		inLine = true;
		return true;
	}
	return false;
}

const string & ConfigFile::getCurrentLine() {
	assert(inLine);
	return currentLine;
}

int ConfigFile::readInteger() {
	assert(inLine);
	int i;
	ss >> i;
	return i;
}

double ConfigFile::readDouble() {
	assert(inLine);
	double d;
	ss >> d;
	return d;
}

string ConfigFile::readString() {
	assert(inLine);
	string s;
	ss >> s;
	return s;
}

bool ConfigFile::readBoolean() {
	string s = readString();
	s = toLowerCase(s);
	if (s == "true" || s == "yes" || s == "1" || s == "on") {
		return true;
	}
	if (s == "false" || s == "no" || s == "0" || s == "off") {
		return false;
	}
	auto message =
			(boost::format("Unable to convert %s to boolean value") % s).str();
	REPORT(Reporter::FATAL) << message;
	throw runtime_error(message);
}

void ConfigFile::process(const string & fileName, function<void(ConfigFile &)>  processor) {
	ConfigFile configFile(fileName);
	while(configFile.readConfigLine()) {
		processor(configFile);
	}
}

} /* namespace GarethUtil */
