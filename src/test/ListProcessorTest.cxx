/*
 * ListProcessorTest.cpp
 *
 *	Test file for testing list processing
 *
 *  Created on: Apr 29, 2014
 *      Author: Gareth Jones
 */


#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"

#include <boost/optional.hpp>
#include <atomic>
#include "../util/Reporter.h"
#include "../util/ListProcessor.h"


using namespace std;
using namespace Gape;

/**
 *  Test routines for list processor template class
 */
TEST_CASE("Check ListProcessor is working", "[ListProcessor]") {

    Reporter::setMinReportingLevel(Reporter::DETAIL);

    SECTION("Check Int List") {}
    {
        vector<int> values{-1, -2, -3, 1, 2, 3, 4, 5, 6};
        ListProcessor<int> list(values);
        REPORT(Reporter::DETAIL) << "list init " << list.toString();
        auto countFn = [](int c, int v) { return ++c; };
        auto count = list.reduce(0, countFn);
        CHECK(count == 9);
        auto filter = [](int v) { return v > 0; };
        list.filter(filter);
        REPORT(Reporter::DETAIL) << "list post filtering " << list.toString();
        CHECK(list.getValues().size() == static_cast<size_t>(6));

        list = ListProcessor<int>(values);
        auto sumFn = [](int c, int v) { return c + v; };
        auto sum = list.reduce(0, sumFn);
        CHECK(sum == 15);

        auto addFn = [](int v) { return ++v; };
        sum = list.map(addFn).reduce(0, sumFn);
        CHECK(sum == 24);
        REPORT(Reporter::DETAIL) << "list post addition " << list.toString();

        map<int, char> map = {{0,  'h'},
                              {-1, 'e'},
                              {-2, 'l'},
                              {2,  'l'},
                              {
                               3,  'o'},
                              {4,  ' '},
                              {5,  'w'},
                              {6,  'o'},
                              {7,  'r'}};
        auto mapFn = [map](int v) { return map.at(v); };
        auto charList = list.map<char>(mapFn);
        REPORT(Reporter::DETAIL) << "list post mapping " << charList.toString();
        auto acc = [](string s, char v) -> string { return s + v; };
        auto str = charList.reduce<string>("", acc);
        REPORT(Reporter::DETAIL) << "Accumulated string is " << str;
        CHECK(str == "hello wor");

        auto matcher = [](int v) -> bool { return v > 1 && v < 3; };
        auto find = list.findFirst(matcher);
        REPORT(Reporter::DETAIL) << "Found " << *find;
        CHECK(find.is_initialized());
        CHECK(2 == *find);

        auto matcher2 = [](int v) -> bool { return v > 20; };
        auto find2 = list.findFirst(matcher2);
        CHECK(!find2.is_initialized());
        CHECK(!find2);

        auto testSum = make_shared<int>(0);
        auto summer =
                [testSum](const int &v) -> void {
                    *testSum += v;
                    cout << "test " << *testSum << " value " << v << endl;
                };
        list.forEach(summer);
        REPORT(Reporter::DETAIL) << "test sum is " << *testSum;
        CHECK(24 == *testSum);

        CHECK(24 == list.sum());
        double average = list.average();
        REPORT(Reporter::DETAIL) << "Average is " << average;
        CHECK(equals(average, 24.0 / 9.0));

        vector<string> strValues{"hello", "world"};

        ListProcessor<string> strList(strValues);
        auto sortFunc = [](const string &s1, const string &s2) -> bool { return s2.compare(s1) < 0; };
        strList.sort(sortFunc);
        REPORT(Reporter::DETAIL) << " sorted str " << strList.toString();
        CHECK(0 == strList.getValues().at(0).compare("world"));
        CHECK(0 == strList.getValues().at(1).compare("hello"));

        // this next line should not compile
        // strList.average();


    }

    SECTION("Check double list") {}
    {
        ListProcessor<double> list({1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
        double sum = list.sum();
        double average = list.average();
        double min = list.min();
        double max = list.max();
        REPORT(Reporter::DETAIL) << "Average is " << average << " sum is " << sum << " min is " << min << " max is "
                                 << max;
        CHECK(equals(sum, 21.0));
        CHECK(equals(average, 21.0 / 6.0));
        CHECK(equals(min, 1.0));
        CHECK(equals(max, 6.0));
    }

}
