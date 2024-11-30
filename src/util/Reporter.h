//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef LOGGER_H_
#define LOGGER_H_

#include <stdlib.h>
#include <iostream>
#include <string>
#include "Util.h"
#include "export.h"

// turn off logging for TRACE and DEBUG completely for optimized code,
// and TRACE for debug code
// (or whatever MIN_REPORTING_LEVEL is set to)

#ifndef MIN_REPORTING_LEVEL
#ifdef NDEBUG
#define MIN_REPORTING_LEVEL Gape::Reporter::ReportingLevel::DETAIL
#else
#define MIN_REPORTING_LEVEL Gape::Reporter::ReportingLevel::DEBUG
#endif
#endif

namespace Gape {

    using namespace std;

/*
 * Our logger class for generic reporting.
 *
 * See http://www.drdobbs.com/cpp/logging-in-c/201804215 for ideas
 */
    class Reporter {
    public:
        enum ReportingLevel {
            TRACE, DEBUG, DETAIL, NORMAL, INFO, WARN, FATAL
        };
    private:
        static ostream *fileStream;
        static ReportingLevel minReportingLevel;
        ReportingLevel reportingLevel;

        Reporter(const Reporter &);

        Reporter &operator=(const Reporter &);

        std::ostringstream os;
    public:
        Reporter() : reportingLevel(NORMAL) {};

        ~Reporter();

        static ostream &getFileStream() {
            return *fileStream;
        }


        static void setFileStream(ostream &fileStream_) {
            fileStream = &fileStream_;
        }

        static ReportingLevel getMinReportingLevel() {
            return minReportingLevel;
        }

        static void setMinReportingLevel(const ReportingLevel reportingLevel_) {
            minReportingLevel = reportingLevel_;
        }

        std::ostringstream &get(ReportingLevel level = NORMAL);

        static std::string levelToString(ReportingLevel level) {
            static const std::string labels[] = {"TRACE", "DEBUG", "DETAIL",
                                                 "NORMAL", "INFO", "WARNING", "FATAL"};
            return labels[level];
        }

        static bool isReportingAt(ReportingLevel level) {
            if (level < MIN_REPORTING_LEVEL) {
                return false;
            } else if (level < getMinReportingLevel()) {
                return false;
            }
            return true;
        }

    };
}

#define REPORT(level) \
if (level < MIN_REPORTING_LEVEL) ;\
else if (level < Gape::Reporter::getMinReportingLevel()) ; \
else Gape::Reporter().get(level) << __FILE__ << ":" << __LINE__ << " "


#endif /* LOGGER_H_ */
