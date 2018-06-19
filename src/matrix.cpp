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

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <config.h>

#include "matrix.h"

using namespace std;

// ============================================================================
// BLAS SINGLE/DOUBLE PRECISION FUNCTION PROTOTYPES
// ============================================================================

#define sgemm_f77 F77_FUNC (sgemm, SGEMM)

// general matrix-matrix multiplication
extern "C" void sgemm_f77(const char* transA, const char* transB,
                          const int* m, const int* n, const int* k,
                          const float* alpha,
                          const float* A, const int* LDA,
                          const float* B, const int* LDB,
                          const float* beta,
                          float *C, const int* LDC);

// ===========================================================================
// MATRIX CLASS
// ===========================================================================

template<>
void Matrix<float>::printSequence(size_t overlap) const
{
        cout << "Print matrix: " << nRows() << " x " << nCols() << endl;
        int K = nCols() / 4 - overlap;
        cout << "K: " << K << ", overlap: " << overlap << endl;

        for (size_t r = 0; r < nRows(); r++) {
                for (size_t c = 0; c < nCols(); c += 4) {
                        if ((*this)(r, c+0) == 1)
                                cout << "A\t";
                        if ((*this)(r, c+1) == 1)
                                cout << "C\t";
                        if ((*this)(r, c+2) == 1)
                                cout << "G\t";
                        if ((*this)(r, c+3) == 1)
                                cout << "T\t";
                }
        }

        cout << endl;
}

/**
 * Perform a matrix-matrix multiplication submatrix(A) * B (float specialized)
 * @param A Left-hand m x k submatrix
 * @param B Right-hand k x n matrix
 */
template <>
void Matrix<float>::gemm(const SubMatrix<float>& A, const Matrix& B)
{
        const float zero = 0.0f;
        const float one = 1.0f;

        int m = A.nRows(), n = B.nCols(), k = B.nRows();
        int LDA = A.nRows(), LDB = B.nRows(), LDC = nRows();
        sgemm_f77("N", "N", &m, &n, &k, &one, A.getDataPtr(), &LDA,
                  B.data, &LDB, &zero, data, &LDC);
}

template <>
void Matrix<float>::gemm(const SubMatrix<float>& A, const Matrix& B,
                         const vector<pair<size_t, size_t> >& matBlocksB)
{
        const float zero = 0.0f;
        const float one = 1.0f;

        int LDA = A.nRows(), LDB = B.nRows(), LDC = nRows();
        int m = A.nRows();

        for (size_t i = 0; i < matBlocksB.size(); i++) {
                size_t colStart = (i == 0) ? 0 : matBlocksB[i-1].first;
                size_t colEnd = matBlocksB[i].first;

                int n = colEnd-colStart;
                int k = matBlocksB[i].second;

                sgemm_f77("N", "N", &m, &n, &k, &one, A.getDataPtr(), &LDA,
                          B.data + colStart*LDB, &LDB, &zero,
                          data + colStart*LDC, &LDC);
        }
}

/**
 * Specialized print matrix elements to the output stream for floats
 * @param os Output stream to add to
 * @param M Matrix to print
 * @return Output stream with the matrix elements
 */
template <>
std::ostream& operator<<(std::ostream& os, const Matrix<float>& M)
{
        std::streamsize oldPrec = os.precision();
        os.precision(2);
        for (size_t r = 0; r < M.nRows(); r++) {
                for (size_t c = 0; c < (M.nCols()-1); c++)
                        os << M.data[c*M.nRows()+r] << "\t";
                os << M.data[(M.nCols()-1)*M.nRows()+r] << endl;
        }
        cout.precision(oldPrec);
        return os;
}
