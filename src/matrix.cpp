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

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <complex>

#include "matrix.h"
#include "blas.h"

using namespace std;

// ===========================================================================
// MATRIX CLASS
// ===========================================================================

/**
 * Specialized gemm function for deci types with BLAS
 * @param A Left-hand n x p matrix
 * @param B Right-hand p x m matrix
 */
template <>
void Matrix<deci>::gemm(const Matrix &A, const Matrix &B, int rowB, int nColB)
{
        assert(isInitialized());
        assert(A.isInitialized());
        assert(B.isInitialized());

        const deci zero = 0;
        const deci one = 1;

        int m = nRows(), n = nColB;
        int LDA = A.nRows(), k = A.nCols(), LDB = B.nRows();
        gemm_f77("N", "N", &m, &n, &k, &one, A.data, &LDA,
                 B.data + rowB, &LDB, &zero, data, &m);
}

/**
 * Specialized gemm function for deci types with BLAS
 * @param A Left-hand n x p matrix
 * @param B Right-hand p x m matrix
 */
template <>
void Matrix<deci>::gemm(const SubMatrix<deci> &A, const Matrix &B)
{
        const deci zero = 0;
        const deci one = 1;

        int m = A.nRows(), n = B.nCols();
        int LDA = A.nRows(), k = B.nRows(), LDB = B.nRows();
        gemm_f77("N", "N", &m, &n, &k, &one, A.getDataPtr(), &LDA,
                 B.data, &LDB, &zero, data, &m);
}

template <>
void Matrix<deci>::gemm(const SubMatrix<deci>& A, const Matrix& B,
                        const vector<pair<size_t, size_t> >& matBlocksA)
{
      /*  const deci zero = 0;
        const deci one = 1;

        int LDA = A.nRows(), LDB = B.nRows(), LDC = nRows();
        int n = B.nCols();

        for (size_t i = 0; i < matBlocksA.size(); i++) {
                size_t start = (i == 0) ? 0 : matBlocksA[i-1].first;
                size_t end = matBlocksA[i].first;

                int m = end-start;
                int k = matBlocksA[i].second;

                gemm_f77("N", "N", &m, &n, &k, &one, A.data + start, &LDA,
                         B.getDataPtr(), &LDB, &zero, data + start, &LDC);
        }*/
}

template <>
double Matrix<dcomplex>::norm() const
{
	double result = 0.0;
	for (int cr=0;cr<nRows();cr++)
		for (int cc=0;cc<nCols();cc++)
			result += abs((*this)(cr,cc)) * abs((*this)(cr,cc));
	return sqrt(result);
}

/**
 * Specialized print matrix elements to the output stream for deci complex
 * @param os Output stream to add to
 * @param M Matrix to print
 * @return Output stream with the matrix elements
 */
template <>
std::ostream& operator<<(std::ostream& os, const Matrix<cplx> &M)
{
        int oldPrec = os.precision();
        os.precision(15);
        for (int r = 0; r < M.nRows(); r++) {
                for (int c = 0; c < (M.nCols()-1); c++) {
                        cplx el = M.data[c*M.nRows()+r];
                        os << real(el) << "\t" << imag(el) << "\t";
                }
                cplx el = M.data[(M.nCols()-1)*M.nRows()+r];
                os << real(el) << "\t" << imag(el) << endl;
        }
        cout.precision(oldPrec);
        return os;
}

/**
 * Specialized print matrix elements to the output stream for decis
 * @param os Output stream to add to
 * @param M Matrix to print
 * @return Output stream with the matrix elements
 */
template <>
std::ostream& operator<<(std::ostream& os, const Matrix<deci> &M)
{
        int oldPrec = os.precision();
        os.precision(15);
        for (int r = 0; r < M.nRows(); r++) {
                for (int c = 0; c < (M.nCols()-1); c++)
                        os << M.data[c*M.nRows()+r] << "\t";
                os << M.data[(M.nCols()-1)*M.nRows()+r] << endl;
        }
        cout.precision(oldPrec);
        return os;
}

template<>
void Matrix<deci>::printSequence(size_t overlap) const
{
        cout << "Print matrix: " << nRows() << " x " << nCols() << endl;
        int K = nRows() / 4 - overlap;
        cout << "K: " << K << ", overlap: " << overlap << endl;

        for (int r = 0; r < nRows(); r += 4) {
                if (r/4 == K)
                        cout << "------------------------\n";
                cout << r/4+1 << "\t";
                for (int c = 0; c < nCols(); c++) {
                        if ((*this)(r+0, c) == 1)
                                cout << "A\t";
                        if ((*this)(r+1, c) == 1)
                                cout << "C\t";
                        if ((*this)(r+2, c) == 1)
                                cout << "G\t";
                        if ((*this)(r+3, c) == 1)
                                cout << "T\t";
                }

                cout << endl;
        }

}
