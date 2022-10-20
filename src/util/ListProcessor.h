/*
 * ListProcessor.h
 *
 *  Created on: Sep 3, 2015
 *      Author: gjones
 */

#ifndef SRC_UTIL_LISTPROCESSOR_H_
#define SRC_UTIL_LISTPROCESSOR_H_

#include <vector>
#include <algorithm>
#include <boost/optional.hpp>
#include <limits>

using namespace std;

namespace GarethUtil {

/**
 * A template to facilitate functional programming using lambda expressions on lists.
 *
 */
template<typename T>
class ListProcessor {
public:
	ListProcessor(const vector<T> & vec_) :
			values(vec_) {
	}

	virtual ~ListProcessor() {
	}

	ListProcessor(const ListProcessor & rhs) = default;
	ListProcessor & operator =(const ListProcessor & rhs) = default;
	ListProcessor(ListProcessor && rhs) = default;
	ListProcessor & operator =(ListProcessor && rhs) = default;

	/**
	 * Remove all items from the list that DON'T match the filter
	 * @param filter
	 * @return
	 */
	ListProcessor<T> & filter(const function<bool(const T &)> & filter) {
		filterList<T>(values, filter);
		return *this;
	}

	/**
	 * Find the first item in the list that matches the function.
	 *
	 * @param matcher
	 * @return
	 */
	boost::optional<T> findFirst(
			const function<bool(const T &)> & matcher) const {
		return findFirstInList<T>(values, matcher);
	}

	/**
	 * Replaces the current list values using the supplied transform function
	 *
	 * @param map
	 * @return
	 */
	ListProcessor<T> & map(const function<T(const T &)> & map) {
		std::transform(values.cbegin(), values.cend(), values.begin(), map);
		return *this;
	}

	/**
	 * Creates a new list of a different type using the supplied transform function
	 *
	 * @param map
	 * @return
	 */
	template<typename V>
	ListProcessor<V> map(const function<V(const T &)> & map) const {
		vector<V> mappedValues;
		mappedValues.reserve(values.size());
		std::transform(values.cbegin(), values.cend(),
				back_inserter(mappedValues), map);
		return ListProcessor<V>(mappedValues);
	}

	/**
	 * Reduces to a single value of the list type using the accumulator function.
	 *
	 * @param startingValue
	 * @param accumulator function takes accumulating value as first argument and current list value as second.
	 * @return
	 */
	T reduce(const T & startingValue,
			const function<T(const T &, const T &)> accumulator) const {
		T currentValue = startingValue;
		for (auto iter = values.cbegin(); iter != values.cend(); ++iter) {
			currentValue = accumulator(currentValue, *iter);
		}
		return currentValue;
	}

	/**
	 * Reduces to a single value of a different type using the accumulator function.
	 *
	 * @param startingValue
	 * @param accumulator function takes accumulating value as first argument and current list value as second.
	 * @return
	 */
	template<typename R>
	R reduce(const R & startingValue,
			const function<R(const R &, const T &)> accumulator) const {
		R currentValue = startingValue;
		for (auto iter = values.cbegin(); iter != values.cend(); ++iter) {
			currentValue = accumulator(currentValue, *iter);
			REPORT(Reporter::DEBUG) << "Current value " << currentValue << "iter"
					<< *iter;
		}
		return currentValue;
	}

	/**
	 * Runs the supplied function for all list values
	 *
	 * @param func
	 */
	ListProcessor<T> & forEach(const function<void(const T &)> & func) {
		for_each(values.cbegin(), values.cend(), func);
		return *this;
	}

	/**
	 * Runs the supplied function for all list values
	 *
	 * @param func
	 */
	const ListProcessor<T> & forEach(
			const function<void(const T &)> & func) const {
		for_each(values.cbegin(), values.cend(), func);
		return *this;
	}

	// unable to create a std::function signature for stl sort- this fails with compilation errors:
	//ListProcessor<T> & sort( function<bool(const  T &a ,  const T &b)> & func) {

	/**
	 * Sorts the list in place.
	 *
	 * @param func
	 * @return
	 */
	template<typename Binop>
	ListProcessor<T> & sort(Binop & func) {
		std::sort(values.begin(), values.end(), func);
		return *this;
	}

	const vector<T>& getValues() const {
		return values;
	}

	const string toString() const {
		return collectionToString(values);
	}

	const size_t count() const {
		return values.size();
	}

	// see
	// http://stackoverflow.com/questions/6972368/stdenable-if-to-conditionally-compile-a-member-function
	// for conditional compilation for arithmetic types only

	/**
	 * Determine the numeric sum- for lists of arithmetic types only
	 * @return
	 */
	template<class A = T>
	typename std::enable_if<std::is_arithmetic<A>::value, const T>::type sum() const {
		auto summer = [] (T s, T v) {return s + v;};
		return reduce(0, summer);
	}

	/**
	 * Determine the numeric average as a double- for lists of arithmetic types only
	 * @return
	 */
	template<class A = T>
	typename std::enable_if<std::is_arithmetic<A>::value, const double>::type average() const {
		T s = sum();
		double sum = static_cast<double>(s);
		double count = static_cast<double>(values.size());
		return sum / count;
	}

	/**
	 * Determine the maximum value- for lists of arithmetic types only
	 * @return
	 */
	template<class A = T>
	typename std::enable_if<std::is_arithmetic<A>::value, const T>::type max() const {
		auto check = [] (T max, T v) {return max > v ? max : v;};
		T min = numeric_limits<T>::is_signed ?
				-numeric_limits<T>::max() : numeric_limits<T>::min();
		return reduce(min, check);
	}

	/**
	 * Determine the minimum value- for lists of arithmetic types only
	 * @return
	 */
	template<class A = T>
	typename std::enable_if<std::is_arithmetic<A>::value, const T>::type min() const {
		auto check = [] (T min, T v) {return min > v ? v : min;};
		T max = numeric_limits<T>::max();
		return reduce(max, check);
	}

private:
	vector<T> values;
};
}

#endif /* SRC_UTIL_LISTPROCESSOR_H_ */
