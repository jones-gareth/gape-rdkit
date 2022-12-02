/*
 * ListProcessorTest.cpp
 *
 *	Test file for testing list processing
 *
 *  Created on: Apr 29, 2014
 *      Author: gjones
 */

#include <gtest/gtest.h>

#include <boost/optional.hpp>
#include "../util/Reporter.h"
#include "../util/ListProcessor.h"
#include <atomic>

namespace {

using namespace std;
using namespace Gape;
using namespace ::testing;

/**
 *  Test routines for list processor template class
 */
class ListProcessorTest: public Test {

protected:
	ListProcessorTest() {
		Reporter::setMinReportingLevel(Reporter::DETAIL);
	}
};

TEST_F(ListProcessorTest, TestIntList) {
	vector<int> values { -1, -2, -3, 1, 2, 3, 4, 5, 6 };
	ListProcessor<int> list(values);
	REPORT(Reporter::DETAIL) << "list init " << list.toString();
	auto countFn = [] (int c, int v) {return ++c;};
	auto count = list.reduce(0, countFn);
	ASSERT_EQ(count, 9);
	auto filter = [] (int v) {return v > 0;};
	list.filter(filter);
	REPORT(Reporter::DETAIL) << "list post filtering " << list.toString();
	ASSERT_EQ(list.getValues().size(), static_cast<size_t>(6));

	list = ListProcessor<int>(values);
	auto sumFn = [] (int c, int v) {return c+v;};
	auto sum = list.reduce(0, sumFn);
	ASSERT_EQ(sum, 15);

	auto addFn = [] (int v) {return ++v;};
	sum = list.map(addFn).reduce(0, sumFn);
	ASSERT_EQ(sum, 24);
	REPORT(Reporter::DETAIL) << "list post addition " << list.toString();

	map<int, char> map = { { 0, 'h' }, { -1, 'e' }, { -2, 'l' }, { 2, 'l' }, {
			3, 'o' }, { 4, ' ' }, { 5, 'w' }, { 6, 'o' }, { 7, 'r' } };
	auto mapFn = [map] (int v) {return map.at(v);};
	auto charList = list.map<char>(mapFn);
	REPORT(Reporter::DETAIL) << "list post mapping " << charList.toString();
	auto acc = [] (string s, char v) -> string {return s+v;};
	auto str = charList.reduce<string>("", acc);
	REPORT(Reporter::DETAIL) << "Accumulated string is " << str;
	ASSERT_EQ(str, "hello wor");

	auto matcher = [](int v)->bool {return v >1 && v <3;};
	auto find = list.findFirst(matcher);
	REPORT(Reporter::DETAIL) << "Found " << *find;
	ASSERT_TRUE(find.is_initialized());
	ASSERT_EQ(2, *find);

	auto matcher2 = [](int v)->bool {return v >20;};
	auto find2 = list.findFirst(matcher2);
	ASSERT_FALSE(find2.is_initialized());
	ASSERT_FALSE(find2);

	auto testSum = make_shared<int>(0);
	auto summer =
			[testSum] (const int & v) -> void {*testSum += v; cout << "test "<<*testSum<< " value "<< v << endl;};
	list.forEach(summer);
	REPORT(Reporter::DETAIL) << "test sum is " << *testSum;
	ASSERT_EQ(24, *testSum);

	ASSERT_EQ(24, list.sum());
	double average = list.average();
	REPORT(Reporter::DETAIL) << "Average is " << average;
	ASSERT_TRUE(equals(average, 24.0 / 9.0));

	vector<string> strValues { "hello", "world" };

	ListProcessor<string> strList(strValues);
	auto sortFunc = [] ( const string & s1,  const string & s2) ->bool {return s2.compare(s1) < 0;};
	strList.sort(sortFunc);
	REPORT(Reporter::DETAIL) << " sorted str " << strList.toString();
	ASSERT_EQ(0, strList.getValues().at(0).compare("world"));
	ASSERT_EQ(0, strList.getValues().at(1).compare("hello"));

	// this next line should not compile
    // strList.average();


}

TEST_F(ListProcessorTest, TestDoubleList) {
	ListProcessor<double> list({1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
	double sum = list.sum();
	double average = list.average();
	double min = list.min();
	double max = list.max();
	REPORT(Reporter::DETAIL) << "Average is " << average <<" sum is " << sum << " min is "<< min<< " max is "<<max;
	ASSERT_TRUE(equals(sum, 21.0));
	ASSERT_TRUE(equals(average, 21.0 / 6.0));
	ASSERT_TRUE(equals(min, 1.0));
	ASSERT_TRUE(equals(max, 6.0));
}

}
