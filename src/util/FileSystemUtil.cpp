
#include "FileSystemUtil.h"

namespace Gape {

using namespace std;

void writeObjectToFile(const string & fileName,
		function<void(ostream &)> func) {
	string test(fileName);
	toLowerCase(test);
	bool gzipped = endsWith(test, ".gz");
	if (gzipped) {
		GzipWriter gzipWriter(fileName);
		func(gzipWriter.getOut());
	} else {
		ofstream out(fileName);
		func(out);
		out.close();
	}
}

RunInDirectory::RunInDirectory(const string & directory_): directory(directory_) {
    if (!boost::filesystem::exists(directory)) {
        boost::filesystem::create_directory(directory);
    }
    startingDirectory = boost::filesystem::current_path();
    boost::filesystem::current_path(directory);
}

RunInDirectory::~RunInDirectory() {
    boost::filesystem::current_path(startingDirectory);
}



} // namespace