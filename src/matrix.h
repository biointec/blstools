/***************************************************************************
 *   Copyright (C) 2017-2018 Jan Fostier (jan.fostier@ugent.be)            *
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
// MATRIX CLASS (Column major storage)
// ===========================================================================

/**
 * Generic template matrix class with some specialized functions for float types
 */
template <class T>
class Matrix
{
public:
        size_t rows;    // number of rows
        size_t cols;    // number of columns
        T *data;        // actual storage for the elements

public:
        /**
         * Default constructor
         */
        Matrix() : rows(0), cols(0), data(NULL) {}

        /**
         * Create a nRows x nCols matrix
         * @param nRows Number of rows in the matrix
         * @param nCols Number of colums in the matrix
         */
        Matrix(size_t nRows, size_t nCols) : rows(nRows), cols(nCols) {
                assert(rows > 0);
                assert(cols > 0);
                data = new T[rows*cols];
        }

        /**
         * Create a nRows x nCols matrix and initialize it
         * @param nRows Number of rows in the matrix
         * @param nCols Number of colums in the matrix
         * @param el Initializer object
         */
        Matrix(size_t nRows, size_t nCols, const T& el) : Matrix(nRows, nCols) {
                fill(el);
        }

        /**
         * Forbid copy constructor
         * @param M Matrix to copy
         */
        Matrix(const Matrix& M) = delete;

        /**
         * Matrix destructor
         */
        ~Matrix() {
                if (data != NULL)
                        delete [] data;
        }

        /**
         * Resize the matrix and initialize it
         * @param nRows Number of rows in the matrix
         * @param nCols Number of colums in the matrix
         * @param el Initializer object
         */
        void resize(size_t nRows, size_t nCols, const T& el) {
                assert(nRows > 0);
                assert(nCols > 0);

                // delete data
                if (data != NULL)
                        delete [] data;

                // create memory for new matrix
                rows = nRows;
                cols = nCols;
                data = new T[rows*cols];
                fill(el);
        }

        /**
         * Fill the matrix with a certain element
         * @param el Element to fill the matrix with
         */
        void fill(const T& el) {
                for (size_t i = 0; i < rows*cols; i++)
                        data[i] = el;
        }

        /**
         * Retrieve the number of rows in the matrix
         * @return Number of rows
         */
        size_t nRows() const {
                return rows;
        }

        /**
         * Retrieve the number of columns in the matrix
         * @return Number of columns
         */
        size_t nCols() const {
                return cols;
        }

        /**
         * Overloaded parentheses to access/modify elements
         * @param row Row specification
         * @param col Column specification
         * @return Element at specified position
         */
        T& operator()(size_t row, size_t col) const {
                return data[col*rows+row];
        }

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
         * Perform a matrix-matrix multiplication submatrix(A) * B
         * @param A Left-hand m x K matrix (K >= k)
         * @param B Right-hand k x n matrix
         */
        void gemm(const SubMatrix<T>& A, const Matrix& B);

        /**
         * Perform a matrix-matrix multiplication submatrix(A) * B
         * @param A Left-hand m x k matrix
         * @param B Right-hand K x n matrix (K >= k)
         * @param matBlocksB Matrix blocks of non-zero elements of B
         */
        void gemm(const SubMatrix<T>& A, const Matrix& B,
                  const std::vector<std::pair<size_t, size_t> >& matBlocksB);

        friend class SubMatrix<T>;
};

// ===========================================================================
// SUBMATRIX CLASS
// ===========================================================================

template <class T>
class SubMatrix
{
public:
        const Matrix<T>& M;     // parent matrix
        size_t firstRow;        // first row index
        size_t rows;            // number of rows
        size_t firstCol;        // first colum index
        size_t cols;            // number of columns

public:
        /**
         * Default constructor
         * @param M Const-reference to parent matrix
         * @param firstRow First row of submatrix
         * @param nRows Number of rows in submatrix
         * @param firstCol First column of submatrix
         * @param nCols Number of columns in submatrix
         */
        SubMatrix(const Matrix<T>& M, size_t firstRow, size_t nRows,
                  size_t firstCol, size_t nCols) : M(M), firstRow(firstRow),
                  rows(nRows), firstCol(firstCol), cols(nCols) {}

        /**
         * Retrieve the number of rows in the matrix
         * @return Number of rows
         */
        size_t nRows() const {
                return rows;
        }

        /**
         * Retrieve the number of columns in the matrix
         * @return Number of columns
         */
        size_t nCols() const {
                return cols;
        }

        /**
         * Retrieve the first row of the submatrix
         * @return First rows
         */
        size_t getFirstRow() const {
                return firstRow;
        }

        /**
         * Retrieve the first columns of the submatrix
         * @return First column
         */
        size_t getFirstCol() const {
                return firstCol;
        }

        /**
         * Retrieve the leading dimensions
         * @return Leading dimensions
         */
        size_t getLD() const {
                return M.nRows();
        }

        /**
         * Get the data pointer to the first element of the submatrix
         * @return Pointer to the first element of the submatrix
         */
        T* getDataPtr() const {
                return M.data + firstCol * M.nRows() + firstRow;
        }
};

#endif
