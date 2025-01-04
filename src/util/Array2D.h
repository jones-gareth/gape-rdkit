/*
 * Array2D.h
 *
 * A simple class for storing and accessing a  2D array
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
#include <sstream>

namespace Gape {

using namespace std;

template<typename T>
class Array2D {

public:

    Array2D(const Array2D & other) = delete;
    Array2D & operator =(const Array2D &) = delete;
    Array2D(Array2D &&) = delete;
    Array2D & operator =(Array2D &&) = delete;

    /**
     * Create a two dimensional array and initialize to default values
     *
     * @param r number of rows
     * @param c number of columns
     */
    Array2D(const size_t r, const size_t c) :
            nRows { r }, nColumns { c }, data(new T[r * c]) {
        initialize();
    }

    void initialize();

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

        return data[r * nColumns + c];
    }

    /**
     * get a value from the array
     *
     * @param r
     * @param c
     * @return
     */
    T & get(const size_t r, const size_t c) {
        assert(r < nRows);
        assert(c < nColumns);

        return data[r * nColumns + c];
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

        data[r * nColumns + c] = value;
        assert(get(r, c) == value);
    }

    /**
     * Function call operator to allow easy getting and setting.
     *
     * @param r
     * @param c
     * @return
     */
    T & operator()(const size_t r, const size_t c) {
        return data[r * nColumns + c];
    }

    /**
     * Function call operator to allow easy getting
     *
     * @param r
     * @param c
     * @return
     */
    const T & operator()(const size_t r, const size_t c) const {
        return data[r * nColumns + c];
    }

    friend ostream& operator<<(ostream& os, const Array2D<T> & arr) {
        os << arr.toString();
        return os;
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
            if (i != nRows - 1) {
                ss << endl;
            }
        }
        ss << "]";
        return ss.str();
    }

    void copyTo(Array2D<T> & other) const {
        assert(other.getRows() == nRows);
        assert(other.getColumns() == nColumns);
        for (size_t i = 0; i < nRows; ++i) {
            for (size_t j = 0; j < nColumns; ++j) {
                other(i, j) = operator()(i, j);
                assert(other.get(i, j) == get(i, j));
            }
        }
    }

    const size_t getColumns() const {
        return nColumns;
    }

    const size_t getRows() const {
        return nRows;
    }

    /**
     * Normalize the array by column vector.  Specialized for doubles only.
     */
    void normalizeColumns();

    /**
     * Determine the rms distance between rows.  Specialized for doubles only.
     * @param i
     * @param j
     * @return
     */
    double rowSqrDistance(int i, int j);

protected:
    /**
     * Constructor that only creates space, but does not set default values
     * @param r
     * @param c
     * @param dummy unused argument to differentiate with other constructor
     */
    Array2D(const size_t r, const size_t c, bool dummy) :
            nRows { r }, nColumns { c }, data(new T[r * c]) {
    }
    ;

    const size_t nRows, nColumns;
    //don't use vector for storage as vector doesn't work right for bool!
    const unique_ptr<T[]> data;
};

/**
 * Type trait to tell if a item is a unique pointer
 */
template<typename T>
struct isUniquePtr {
    static const bool value = false;
};

template<typename T>
struct isUniquePtr<unique_ptr<T>> {
    static const bool value = true;
};

/**
 * Class to initialize an Array2D to default values
 */
template<bool isUniquePtr>
struct Array2DInitializer {
    template<typename T>
    static void initialize(Array2D<T> & array) {
        T defaultValue { };
        for (auto i = 0ul; i < array.getRows(); i++) {
            for (auto j = 0ul; j < array.getColumns(); j++) {
                array.get(i, j) = defaultValue;
            }
        }
    }
};

/**
 * Class to initialize an Array2D of unique pointers to null pointers
 */
template<>
struct Array2DInitializer<true> {
    template<typename T>
    static void initialize(Array2D<T> & array) {
        for (auto i = 0ul; i < array.getRows(); i++) {
            for (auto j = 0ul; j < array.getColumns(); j++) {
                array.get(i, j) = nullptr;
            }
        }
    }
};

/**
 * Method that uses type traits to choose the correct initializer implementation
 */
template<typename T>
void Array2D<T>::initialize() {
    Array2DInitializer<isUniquePtr<T>::value>::initialize(*this);
}

/**
 * Normalize the array by column vector.  Specialized for doubles only.
 */
template<>
void Array2D<double>::normalizeColumns();

/**
 * Determine the rms distance between rows.  Specialized for doubles only.
 * @param i
 * @param j
 * @return
 */
template<>
double Array2D<double>::rowSqrDistance(const int i, const int j);

}


/* namespace Gape */

#endif /* ARRAY2D_H_ */
