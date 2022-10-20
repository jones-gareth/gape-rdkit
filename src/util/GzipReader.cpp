/*
 * GzipReader.cpp
 *
 *  Created on: May 2, 2014
 *      Author: gjones
 */

#include "GzipReader.h"

namespace GarethUtil {

using namespace std;

GzipReader::GzipReader(string file) : inFile(file) {
	in.push(boost::iostreams::gzip_decompressor());
	//in.push(boost::iostreams::file_descriptor_source(inFile));
	inFileStream.open(inFile, std::ios_base::in | std::ios_base::binary);
	if (!inFileStream.is_open()) {
		throw runtime_error("Failed to opene compressed file "+file);
	}
	in.push(inFileStream);

}

GzipReader::~GzipReader() {
	in.reset();
	inFileStream.close();
}

} /* namespace GarethUtil */
