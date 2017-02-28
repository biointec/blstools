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

#ifndef BLASMATHPACK_H
#define BLASMATHPACK_H

#include "base.h"

// ============================================================================
// BLAS SINGLE/DOUBLE PRECISION FUNCTION PROTOTYPES
// ============================================================================

#ifdef USE_FLOAT
        #define gemm_f77 sgemm_f77
#else
        #define gemm_f77 dgemm_f77
#endif

#define sgemm_f77 F77_FUNC (sgemm, SGEMM)

// general matrix-matrix multiplication
extern "C" void sgemm_f77(const char* transA, const char* transB, const int* M,
                          const int* N, const int* K, const float* alpha,
                          const float* A, const int* lda, const float *B,
                          const int *ldb, const float *beta, float *C,
                          const int*ldc);

#endif
