/*
 * Array2D.h
 *
 * A simple class for storing and accessing a  2D array
 *
 * This uses a vector as the backing storage, but I've stopped using this
 * as vector<bool> does not work right in the STL
 *
 *  Created on: May 21, 2014
 *      Author: Gareth Jones
 *
 *
 */

#ifndef ARRAY2D_H_
#define ARRAY2D_H_

#include <array>
#include <cassert>
#include <memory>

namespace GarethUtil {

using namespace std;

template<typename T>
class Array2D {
public:

	Array2D(const Array2D & other) = delete;
	Array2D & operator =(const Array2D &) = delete;
	Array2D(Array2D &&) = delete;
	Array2D & operator =(Array2D &&) = delete;

	/**
	 * Create a two dimensional array
	 *
	 * @param r number of rows
	 * @param c number of columns
	 */
	Array2D(const size_t r, const size_t c) :
			nRows { r }, nColumns { c }, data(r * c) {
	}

	virtual ~Array2D() {
	}

	/**
	 * get a value from the array
	 *
	 * @param r
	 * @param c
	 * @return
	 */
	const T & get(const size_t r, const size_t c) const {
		assert(r < nRows);
		assert(c < nColumns);

		const T & val = data.at(r * nColumns + c);
		return val;
	}

	/**
	 * get a value from the array
	 *
	 * @param r
	 * @param c
	 * @return
	 */
	T  get(const size_t r, const size_t c) {
		assert(r < nRows);
		assert(c < nColumns);

		T  val = data.at(r * nColumns + c);
		return val;
	}

	/**
	 * set a value in the array
	 * @param r
	 * @param c
	 * @param value
	 */
	void set(const size_t r, const size_t c, const T & value) {
		assert(r < nRows);
		assert(c < nColumns);

		data.at(r * nColumns + c) = value;
	}

	/**
	 * Function call operator to allow easy getting and setting.
	 *
	 * @param r
	 * @param c
	 * @return
	 */
	T & operator()(const size_t r, const size_t c) {
		return data.at(r * nColumns + c);
	}

	/**
	 * Function call operator to allow easy getting
	 *
	 * @param r
	 * @param c
	 * @return
	 */
//	const T & operator()(const size_t r, const size_t c) const {
//		return data.at(r * nColumns + c);
//	}

	void copyTo(Array2D<T> & other) const {
		assert(other.getRows() == nRows);
		assert(other.getColumns() == nColumns);
		for (size_t i = 0; i < nRows; ++i) {
			for (size_t j = 0; j < nColumns; ++j) {
				other.set(i, j, get(i, j));
			}
		}
	}

	const size_t getColumns() const {
		return nColumns;
	}

	const size_t getRows() const {
		return nRows;
	}

	string toString() const {
		stringstream ss;
		ss << "[";
		for (size_t i = 0; i < nRows; ++i) {
			ss << "[";
			for (size_t j = 0; j < nColumns; ++j) {
				ss << get(i, j);
				if (j != nColumns - 1) {
					ss << ", ";
				}
			}
			ss << "]";
			if (i != nRows-1) {
				ss << endl;
			}
		}
		ss << "]";
		ss << endl;
		return ss.str();
	}
private:
	const size_t nRows, nColumns;
	//don't use vector for storage as vector doesn't work right for bool!
	vector<T> data;
};

} /* namespace GarethUtil */

#endif /* ARRAY2D_H_ */
