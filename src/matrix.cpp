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
        size_t K = nRows() / 4 - overlap;
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

// ===========================================================================
// VECTOR CLASS
// ===========================================================================

/**
 * Specialized copy function for deci complex types with BLAS
 * @param X Right-hand side vector
 */
template <>
void Vector<cplx>::copy(const Vector &X)
{
        // check dimensions
        assert(X.size() == size());

        // cblas_zcopy(size, X.data, 1, data, 1);
        memcpy(data, X.data, size()*sizeof(cplx));
}

/**
 * Print vector elements to the output stream
 * @param os Output stream to add to
 * @param V Vector to print
 * @return Output stream with the vector elements
 */
template <>
std::ostream& operator<<(std::ostream& os, const Vector<cplx>& V)
{
        int oldPrec = os.precision();
        os.precision(15);
        for (int r = 0; r < V.size(); r++) {
                cplx el = V.data[r];
                os << real(el) << "\t" << imag(el) << endl;
        }
        os.precision(oldPrec);
        return os;
}

template <>
std::istream& operator>>(std::istream& is, const Vector<cplx>& V)
{
        for (int r = 0; r < V.size(); r++) {
                double real, imag;
                is >> real >> imag;
                V.data[r] = cplx(real, imag);
        }
        return is;
}

/**
 * Print vector elements to the output stream
 * @param os Output stream to add to
 * @param V Vector to print
 * @return Output stream with the vector elements
 */
template <>
std::ostream& operator<<(std::ostream& os, const Vector<deci>& V)
{
        int oldPrec = os.precision();
        os.precision(15);
        for (int r = 0; r < V.size(); r++)
                os << V.data[r] << endl;
        os.precision(oldPrec);
        return os;
}