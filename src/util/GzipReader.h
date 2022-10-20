/*
 * GzipReader.h
 *
 * A simple class that uses boost io streams to manage reading from gzipped files.
 * Don't let this get out of scope until you've finished all reading.
 *
 *  Created on: May 2, 2014
 *      Author: gjones
 */

#ifndef GZIPREADER_H_
#define GZIPREADER_H_

#include <iostream>
#include <fstream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

namespace GarethUtil {

using namespace std;

class GzipReader {
public:
	explicit GzipReader(string file);
	virtual ~GzipReader();

	/**
	 *
	 * @return stream to read from
	 */
	istream& getIn() {
		return in;
	}

private:
	boost::iostreams::filtering_istream in;
	ifstream inFileStream;
	string inFile;
};

} /* namespace GarethUtil */

#endif /* GZIPREADER_H_ */
