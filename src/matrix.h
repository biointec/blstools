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
#include <vector>

// ===========================================================================
// TEMPLATE PROTOTYPES
// ===========================================================================

template <class T>
class Matrix;

template <class T>
class SubMatrix;

template <class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& M);

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
         * Perform a matrix-matrix multiplication A * submatrix(B)
         * @param A Left-hand m x k matrix
         * @param B Right-hand K x n matrix (K >= k)
         * @param rowB First row of B to consider
         * @param nColB Number of column of B to consider
         */
        void gemm(const Matrix &A, const Matrix &B, int rowB, int nColB);

        /**
         * Perform a matrix-matrix multiplication A * submatrix(B)
         * @param A Left-hand m x k matrix
         * @param B Right-hand K x n matrix (K >= k)
         */
        void gemm(const Matrix &A, const SubMatrix<T> &B);

        /**
         * Perform a matrix-matrix multiplication A * submatrix(B)
         * @param A Left-hand m x k matrix
         * @param B Right-hand K x n matrix (K >= k)
         */
        void gemm(const Matrix &A, const SubMatrix<T> &B,
                  const std::vector<std::pair<size_t, size_t> >& matBlocksA);

        /** returns the 2-norm of the matrix
        */
        double norm() const;

        friend class SubMatrix<T>;
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
Matrix<T>::~Matrix()
{
        if (data != NULL)
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
void Matrix<T>::fill(const T& el)
{
        assert(isInitialized());
        for (int i = 0; i < nRows()*nCols(); i++)
                data[i] = el;
}

// ===========================================================================
// SUBMATRIX CLASS
// ===========================================================================

template <class T>
class SubMatrix
{
public:
        const Matrix<T>& M;     // parent matrix
        int firstRow;           // first row index
        int rows;               // number of rows
        int firstCol;           // first colum index
        int cols;               // number of columns

public:
        /**
         * Default constructor
         */
        SubMatrix(const Matrix<T>& M, int firstRow, int nRows, int firstCol,
                  int nCols) : M(M), firstRow(firstRow), rows(nRows),
                  firstCol(firstCol), cols(nCols) {}

        /**
         * Retrieve the number of rows in the matrix
         * @return Number of rows
         */
        int nRows() const { return abs(rows); }

        /**
         * Retrieve the number of columns in the matrix
         * @return Number of columns
         */
        int nCols() const { return abs(cols); }

        /**
         * Retrieve the first row of the submatrix
         * @return First rows
         */
        int getFirstRow() const { return firstRow; }

        /**
         * Retrieve the first columns of the submatrix
         * @return First column
         */
        int getFirstCol() const { return firstCol; }

        /**
         * Retrieve the leading dimensions
         * @return Leading dimensions
         */
        int getLD() const { return M.nRows(); };

        /**
         * Get the data pointer to the first element of the submatrix
         *
         */
        T* getDataPtr() const {
                return M.data + firstCol * M.nRows() + firstRow;
        }
};

#endif
