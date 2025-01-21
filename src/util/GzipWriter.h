/*
 * GzipWriter.h
 *
 * A simple class that uses boost to handle writing to gzipped files.
 * Don't let this object go out of scope until you have finished writing.
 *
 *  Created on: May 2, 2014
 *      Author: Gareth Jones
 *
 */


#ifndef GZIPWRITER_H_
#define GZIPWRITER_H_

#include <iostream>
#include <fstream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

namespace Gape {

using namespace std;

class GzipWriter {
public:
	explicit GzipWriter(const string &file);
	virtual ~GzipWriter();

	/**
	 * Do writing on this stream
	 * @return
	 */
	ostream& getOut()  {
		return out;
	}

private:
	boost::iostreams::filtering_ostream out;
	string outFile;
	boost::iostreams::file_descriptor_sink fileDescriptorSink;
	//ofstream outFileStream;
};

} /* namespace Gape */

#endif /* GZIPWRITER_H_ */
