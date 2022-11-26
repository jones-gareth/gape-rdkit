/*
 * ConfigFile.h
 *
 *  Created on: Jan 21, 2016
 *      Author: gjones
 */

#ifndef SRC_UTIL_CONFIGFILE_H_
#define SRC_UTIL_CONFIGFILE_H_

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>

namespace GarethUtil {

using namespace std;

/**
 * Helper class to facilitate the reading of configuration files
 */
class ConfigFile {
public:
	ConfigFile(const string & fileName_) ;

	virtual ~ConfigFile() {
		in.close();
	}

	ConfigFile(const ConfigFile & rhs) = delete;
	ConfigFile & operator =(const ConfigFile & rhs) = delete;
	ConfigFile(ConfigFile && rhs) = delete;
	ConfigFile & operator =(ConfigFile && rhs) = delete;

	/**
	 * Reads the next line from the configuration file. Returns false if eof.
	 * @return
	 */
	bool readConfigLine();

	/**
	 * Reads the next integer from the current line
	 * @return
	 */
	int readInteger();

	/**
	 * Reads the next double from the current line
	 *
	 * @return
	 */
	double readDouble();

	/**
	 * Reads the next string from the current line
	 *
	 * @return
	 */
	string readString();

	/**
	 * Reads the next boolean from the current line
	 * @return
	 */
	bool readBoolean();

	/**
	 * Retruns the current line
	 *
	 * @return
	 */
	const string & getCurrentLine();

	/**
	 * A convenience function to facilitate the processing of configuration file
	 *
	 * @param fileName
	 * @param processor function to process each line
	 */
	static void process(const string & fileName, function<void(ConfigFile &)> processor);

	/**
	 * Return the path to the configuration directory
	 *
	 * @return
	 */
	static const string configDirectory();


	/**
	 * Returns the full native path of a configuration file
	 * @param configFile
	 * @return
	 */
	static string configFilePath(const string &configFile);



private:
	const string  & fileName;
	const string & fullPath;
	 ifstream in;
	 istringstream ss;
	 string currentLine;
	 bool inLine;


};
} /* namespace GarethUtil */

#endif /* SRC_UTIL_CONFIGFILE_H_ */
