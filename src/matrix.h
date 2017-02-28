/***************************************************************************
 *   Copyright (C) 2017 Jan Fostier (jan.fostier@ugent.be)                 *
 *   This file is part of BLStools                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cassert>

// ===========================================================================
// TEMPLATE PROTOTYPES
// ===========================================================================

template <class T>
class Vector;

template <class T>
class Matrix;

template <class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& M);

template <class T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& V);

template <class T>
std::istream& operator>>(std::istream& is, const Vector<T>& V);

// ===========================================================================
// MATRIX CLASS
// ===========================================================================

/**
 * Generic template matrix class with some specialized functions for
 * the complex deci type which use the BLAS and LAPACK
 */
template <class T>
class Matrix
{
public:
        int rows;       // number of rows
        int cols;       // number of columns
        T *data;        // actual storage for the elements

public:
        /**
         * Default constructor
         */
        Matrix();

        /**
         * Creates a nRows x nCols matrix
         * @param nRows Number of rows in the matrix
         * @param nCols Number of colums in the matrix
         */
        Matrix(int nRows, int nCols);

        /**
         * Creates a nRows x nCols matrix and initializes it
         * @param nRows Number of rows in the matrix
         * @param nCols Number of colums in the matrix
         * @param el Initializer object
         */
        Matrix(int nRows, int nCols, const T& el);

        /**
         * Copy constructor for matrix
         * @param M Matrix to copy
         */
        Matrix(const Matrix& M);

        /**
         * Reinterpret vector as a matrix (shallow copy)
         * @param V Original vector
         * @param nCols Number of columns
         */
        Matrix(const Vector<T>& V, int nCols);

        /**
         * Matrix destructor
         */
        ~Matrix();

        /**
         * Retrieve the number of rows in the matrix
         * @return Number of rows
         */
        int nRows() const { return abs(rows); };

        /**
         * Retrieve the number of columns in the matrix
         * @return Number of columns
         */
        int nCols() const { return abs(cols); };

        /**
         * Set the matrix dimensions without initializing memory
         * @param nRows Number of rows in the matrix
         * @param nCols Number of colums in the matrix
         */
        void setDimensions(int nRows, int nCols);

        /**
         * Allocate memory for the matrix
         */
        void allocateMemory();

        /**
         * Allocate memory for the matrix and initialize it
         * @param el Initializer object
         */
        void allocateMemory(const T& el);

        /**
         * Returns true if the data pointer has been initialized
         * @return True of false
         */
        bool isInitialized() const;

        /**
         * Overloaded parentheses to access/modify elements
         * @param row Row specification
         * @param col Column specification
         * @return Element at specified position
         */
        T& operator()(int row, int col) const;

        /**
         * Overloaded assignment operator
         * @param M Right hand side matrix
         * @return Reference to the lvalue matrix
         */
        Matrix& operator=(const Matrix& M);

        /**
         * Matrix-matrix multiplication
         * @param M Right hand side matrix
         * @return Product of the two matrices
         */
        const Matrix operator*(const Matrix &M);

        /**
         * Matrix-matrix addition
         * @param M Right hand side matrix
         * @return Sum of the two matrices
         */
        const Matrix operator+(const Matrix &M);

        /**
         * Matrix-matrix subtraction
         * @param M Right hand side matrix
         * @return Difference of the two matrices
         */
        const Matrix operator-(const Matrix &M);

        /**
         * Matrix-matrix addition
         * @param M Right hand side vector
         */
        void operator+=(const Matrix &M);

        /**
         * Scaling of a matrix
         * @param M Right hand size element
         */
        void operator*=(const T &M);

        /**
         * Matrix-vector multiplication
         * @param V Right hand side matrix
         * @return Product of the matrix-vector
         */
        const Vector<T> operator*(const Vector<T> &V);

        /**
         * Print matrix elements to the output stream
         * @param os Output stream to add to
         * @param M Matrix to print
         * @return Output stream with the matrix elements
         */
        friend std::ostream& (::operator<< <>)(std::ostream& os,
                                               const Matrix<T>& M);

        /**
         * Print the sequence imposed by the matrix
         * @param overlap Number of overlapping nucleotides
         */
        void printSequence(size_t overlap) const;

        /**
         * Fill the matrix with a certain element
         * @param el Element to fill the matrix with
         */
        void fill(const T& el);

        /**
         * Calculate the inverse of the matrix
         */
        void invert();

        /**
         * Perform a matrix-matrix multiplication A * submatrix(B)
         * @param A Left-hand m x k matrix
         * @param B Right-hand K x n matrix (K >= k)
         * @param rowB First row of B to consider
         * @param nColB Number of column of B to consider
         */
        void gemm(const Matrix &A, const Matrix &B, int rowB, int nColB);

        /** returns the 2-norm of the matrix
        */
        double norm() const;

        friend class Vector<T>;
};

template <class T>
Matrix<T>::Matrix() : rows(0), cols(0), data(NULL)
{

}

template <class T>
Matrix<T>::Matrix(int nRows, int nCols) : rows(nRows), cols(nCols)
{
        assert(rows > 0);
        assert(cols > 0);
        data = new T[rows*cols];
}

template <class T>
Matrix<T>::Matrix(int nRows, int nCols, const T& el) : rows(nRows), cols(nCols)
{
        assert(rows > 0);
        assert(cols > 0);
        data = new T[rows*cols];
        fill(el);
}

template <class T>
Matrix<T>::Matrix(const Matrix& M)
{
        assert(M.isInitialized());
        cols = M.nCols();
        rows = M.nRows();
        data = new T[rows*cols];
        for (int i = 0; i < rows*cols; i++)
                data[i] = M.data[i];
}

template <class T>
Matrix<T>::Matrix(const Vector<T>& V, int nCols)
{
        assert(V.isInitialized());
        assert(nCols > 0);
        data = V.data;
        rows = -V.size()/nRows;
        cols = -nCols;
}

template <class T>
Matrix<T>::~Matrix()
{
        if ((rows > 0) && (cols > 0))
               delete [] data;
}

template <class T>
void Matrix<T>::setDimensions(int nRows, int nCols)
{
        assert(nRows >= 0);
        assert(nCols >= 0);
        // delete data
        if ((rows > 0) && (cols > 0))
               delete [] data;
        data = NULL;
        rows = nRows;
        cols = nCols;
}

template <class T>
void Matrix<T>::allocateMemory()
{
        assert(data == NULL);
        if (rows*cols > 0)
                data = new T[rows*cols];
}

template <class T>
void Matrix<T>::allocateMemory(const T& el)
{
        allocateMemory();
        fill(el);
}

template <class T>
bool Matrix<T>::isInitialized() const
{
        if (rows*cols == 0) return true;
        return (data != NULL);
}

template <class T>
inline T& Matrix<T>::operator()(int row, int col) const
{
        assert(isInitialized());
        assert(row >= 0); assert(row < nRows());
        assert(col >= 0); assert(col < nCols());
        return data[col*nRows()+row];
}

template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix& M)
{
        assert(M.isInitialized());
        // check for self-assignment
        if (&M == this) return *this;
        // if the dimensions do not match: change them
        if ((nRows() != M.nRows()) || (nCols() != M.nCols())) {
                if ((rows > 0) && (cols > 0)) delete [] data;
                rows = M.nRows();
                cols = M.nCols();
                data = new T[rows*cols];
        }
        // copy element per element
        for (int i = 0; i < rows*cols; i++)
                data[i] = M.data[i];
        return *this;
}

template <class T>
const Matrix<T> Matrix<T>::operator*(const Matrix &M)
{
        assert(isInitialized());
        assert(M.isInitialized());
        assert(nCols() == M.nRows());
        Matrix result(nRows(), M.nCols(), 0);

        for (int i = 0; i < nRows(); i++) {
                for (int j = 0; j < M.nCols(); j++)
                        for(int k = 0; k < nCols(); k++)
                                result(i, j) += (*this)(i, k)*M(k, j);
        }
        return result;
}

template <class T>
const Matrix<T> Matrix<T>::operator+(const Matrix &M)
{
        assert(isInitialized());
        assert(M.isInitialized());
        assert(nCols() == M.nCols());
        assert(nRows() == M.nRows());
        Matrix result(nRows(), nCols());

        for (int i = 0; i < nRows()*nCols(); i++)
                result.data[i] = data[i] + M.data[i];
        return result;
}

template <class T>
void Matrix<T>::operator+=(const Matrix &M)
{
        assert(isInitialized());
        assert(M.isInitialized());
        assert(nCols() == M.nCols());
        assert(nRows() == M.nRows());
        for (int i = 0; i < nRows()*nCols(); i++)
                data[i] += M.data[i];
}

template <class T>
void Matrix<T>::operator*=(const T &M)
{
        assert(isInitialized());
        for (int i = 0; i < nRows()*nCols(); i++)
                data[i] *= M;
}

template <class T>
const Matrix<T> Matrix<T>::operator-(const Matrix &M)
{
        assert(isInitialized());
        assert(M.isInitialized());
        assert(nCols() == M.nCols());
        assert(nRows() == M.nRows());
        Matrix result(nRows(), nCols());

        for (int i = 0; i < nRows()*nCols(); i++)
                result.data[i] = data[i] - M.data[i];
        return result;
}

template <class T>
const Vector<T> Matrix<T>::operator*(const Vector<T> &V)
{
        assert(isInitialized());
        assert(nCols() == V.size());
        Vector<T> result(nRows(), 0);

        for (int i = 0; i < nRows(); i++) {
                for (int j = 0; j < nCols(); j++)
                        result(i) += (*this)(i, j)*V(j);
        }
        return result;
}

template <class T>
void Matrix<T>::fill(const T& el)
{
        assert(isInitialized());
        for (int i = 0; i < nRows()*nCols(); i++)
                data[i] = el;
}

// ===========================================================================
// VECTOR CLASS
// ===========================================================================

/**
 * Generic template vector class with some specialized functions for
 * the complex deci type which use the BLAS and LAPACK.
 */
template <class T>
class Vector
{
public:
        int sz;        // number of elements
        T *data;       // actual storage for the elements

public:
        /**
         * Default constructor
         */
        Vector();

        /**
         * Create a vector of a certain size
         * @param size Number of elements in the vector
         */
        Vector(int size);

        /**
         * Create a vector of a certain size and initialize it
         * @param size Number of elements in the vector
         * @param el Initializer object
         */
        Vector(int size, const T& el);

        /**
         * Copy constructor for vector
         * @param V Vector to copy
         */
        Vector(const Vector& V);

        /**
         * Get a subvector (shallow copy)
         * @param V Original vector
         * @param first First element
         * @param last Last element
         */
        Vector(const Vector& V, int first, int last);

        /**
         * Create a vector wrapped from T pointer
         * @param Data Pointer to the data
         * @param size Number of elements
         */
        Vector(T *data, int size);

        /**
         * Vector destructor
         */
        ~Vector();

        /**
         * Return the size of the vector
         * @return Size of the vector
         */
        int size() const { return abs(sz); };

        /**
         * Set the size of the vector without initializing memory
         * @param size Number of elements in the vector
         */
        void setDimensions(int size);

        /**
         * Allocate memory for the vector
         */
        void allocateMemory();

        /**
         * Allocate memory for the vector and initialize it
         * @param el Initializer object
         */
        void allocateMemory(const T& el);

        /**
         * Returns true if the data pointer has been initialized
         * @return True of false
         */
        bool isInitialized() const;

        /**
          * Overloaded parentheses to access/modify elements
          * @param row Row specification
          * @return Element at specified position
          */
        T& operator()(int row) const;

        /**
          * Overloaded brackets to access/modify elements
          * @param row Row specification
          * @return Element at specified position
          */
        T& operator[](int row) const;

        /**
         * Overloaded assignment operator
         * @param V Right hand side vector
         * @return Reference to the lvalue vector
         */
        Vector& operator=(const Vector& V);

        /**
         * Vector-vector addition
         * @param V Right hand side vector
         * @return Sum of the two vectors
         */
        const Vector operator+(const Vector &V);

        /**
         * Vector-vector subtraction
         * @param V Right hand side vector
         * @return Subtraction of the two vectors
         */
        const Vector operator-(const Vector &V);

        /**
         * Scaling of a vector
         * @param V Right hand size element
         */
        void operator*=(const T &V);

        /**
         * Print vector elements to the output stream
         * @param os Output stream to add to
         * @param V Vector to print
         * @return Output stream with the vector elements
         */
        friend std::ostream& (::operator<< <>)(std::ostream& os,
                                               const Vector& V);

        /**
         * Load vector elements to from output stream
         * @param is Input stream to add to
         * @param V Vector to load
         * @return Input stream the vector elements removed
         */
        friend std::istream& (::operator>> <>)(std::istream& is,
                                               const Vector& V);

        /**
         * Fill the vector with a certain element
         * @param el Element to fill the vector with
         */
        void fill(const T& el);

        /**
         * Copy vector X to vector Y
         * @param X Right-hand vector
         */
        void copy(const Vector &X);

        friend class Matrix<T>;
};

template <class T>
Vector<T>::Vector() : sz(0), data(NULL)
{

}

template <class T>
Vector<T>::Vector(int size) : sz(size)
{
        assert(sz >= 0);
        if (sz == 0) return;
        data = new T[sz];
}

template <class T>
Vector<T>::Vector(int size, const T& el) : sz(size)
{
        assert(sz >= 0);
        if (sz == 0) return;
        data = new T[sz];
        fill(el);
}

template <class T>
Vector<T>::Vector(const Vector& V)
{
        assert(V.isInitialized());
        sz = V.size();
        data = new T[sz];
        for (int i = 0; i < sz; i++)
                data[i] = V.data[i];
}

template <class T>
Vector<T>::Vector(const Vector& V, int first, int last)
{
        assert(V.isInitialized());
        assert(first <= last);
        assert(first >= 0);
        assert(last = V.size());
        data = V.data + first;
        sz = -(last-first+1);
}

template <class T>
Vector<T>::Vector(T *data_, int size)
{
        assert(size > 0);
        if (size == 0) return;
        data = data_;
        sz = -size;
}

template <class T>
Vector<T>::~Vector()
{
        if (sz > 0)
                delete [] data;
}

template <class T>
void Vector<T>::setDimensions(int size)
{
        // delete data
        if (sz > 0) {
                delete [] data;
                data = NULL;
        }

        assert(size >= 0);
        sz = size;
}

template <class T>
void Vector<T>::allocateMemory()
{
        assert(data == NULL);
        if (sz > 0)
                data = new T[sz];
}

template <class T>
void Vector<T>::allocateMemory(const T& el)
{
        allocateMemory();
        fill(el);
}

template <class T>
bool Vector<T>::isInitialized() const
{
        if (sz == 0) return true;
        return (data != NULL);
}

template <class T>
inline T& Vector<T>::operator()(int row) const
{
        assert(isInitialized());
        assert(row >= 0); assert(row < size());
        return data[row];
}

template <class T>
inline T& Vector<T>::operator[](int row) const
{
        assert(isInitialized());
        assert(row >= 0); assert(row < size());
        return data[row];
}

template <class T>
Vector<T>& Vector<T>::operator=(const Vector& V)
{
        assert(V.isInitialized());
        // check for self-assignment
        if (&V == this) return *this;
        // if the dimension does not match: change it
        if (size() != V.size()) {
                if (sz > 0) delete [] data;
                sz = V.size();
                data = new T[sz];
        }

        for (int i = 0; i < sz; i++)
                data[i] = V.data[i];
        return *this;
}

template <class T>
const Vector<T> Vector<T>::operator+(const Vector &V)
{
        assert(isInitialized());
        assert(V.isInitialized());
        assert(size() == V.size());
        Vector result(size());

        for (int i = 0; i < size(); i++)
                result.data[i] = data[i] + V.data[i];
        return result;
}

template <class T>
const Vector<T> Vector<T>::operator-(const Vector &V)
{
        assert(isInitialized());
        assert(V.isInitialized());
        assert(size() == V.size());
        Vector result(size());

        for (int i = 0; i < size(); i++)
                result.data[i] = data[i] - V.data[i];
        return result;
}

template <class T>
void Vector<T>::operator*=(const T &E)
{
        assert(isInitialized());
        for (int i = 0; i < size(); i++)
                data[i] *= E;
}

template <class T>
void Vector<T>::fill(const T& el)
{
        assert(isInitialized());
        for (int i = 0; i < size(); i++)
                data[i] = el;
}

#endif
