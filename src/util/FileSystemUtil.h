
#ifndef FS_UTIL_H_
#define FD_UTIL_H_

#include <string>
#include <boost/filesystem.hpp>

#include "GzipWriter.h"
#include "GzipReader.h"
#include "Util.h"

namespace GarethUtil {

/**
 * Write to a gzipped or normal file using the supplied function
 * @param fileName
 * @param func
 * @param object
 */
void writeObjectToFile(const string & fileName, function<void(ostream &)> func);

/**
 * Read an object from a gzipped or normal file using the supplied function.
 *
 * @param fileName
 * @param func
 * @return
 */
template<typename T>
T readObjectFromFile(const string & fileName, function<T(istream &)> func) {
    string test(fileName);
    toLowerCase(test);
    bool gzipped = endsWith(test, ".gz");

    ifstream in;
    if (gzipped) {
        GzipReader gzipReader(fileName);
        auto & in = gzipReader.getIn();
        return func(in);
    } else {
        ifstream in(fileName);
        if (!in.is_open()) {
            throw runtime_error("Failed to open fileName " + fileName);
        }
        T object = func(in);
        in.close();
        return object;
    }
}

/**
 * A RAII class to handle temporary change of directory. Will atttempt to make the directory if it doesn't exist.
 */
class RunInDirectory {
public:
    RunInDirectory(const string & directory);

    virtual ~RunInDirectory();

    RunInDirectory(const RunInDirectory & rhs) = delete;
    RunInDirectory & operator =(const RunInDirectory & rhs) = delete;
    RunInDirectory(RunInDirectory && rhs) = delete;
    RunInDirectory & operator =(RunInDirectory && rhs) = delete;

private:
    const string & directory;
    boost::filesystem::path startingDirectory;
};

} //namespace


#endif /* FS_UTIL_H_ */