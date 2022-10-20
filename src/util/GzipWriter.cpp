/*
 * GzipWriter.cpp
 *
 *  Created on: May 2, 2014
 *      Author: gjones
 */

#include "GzipWriter.h"

namespace GarethUtil {
using namespace std;

GzipWriter::GzipWriter(string file) :
		outFile(file), fileDescriptorSink(file) {
	// use the file_descriptor_sink instead of ofstream as that seems to hang output.
	//outFileStream.open(outFile, ios_base::out | ios_base::binary);
	out.push(boost::iostreams::gzip_compressor());
	out.push(fileDescriptorSink);
	//out.push(outFileStream);

}

GzipWriter::~GzipWriter() {
	//outFileStream.close();
	out.reset();
	// This has to be done after any pop and reset calls (otherwise you get an exception).
	// According  to the docs reset should close any device, this seems safe to do.
	fileDescriptorSink.close();
}

} /* namespace GarethUtil */
