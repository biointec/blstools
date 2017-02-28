/***************************************************************************
 *   Copyright (C) 2006, 2007, 2008, 2009 by Jan Fostier                   *
 *   jan.fostier@intec.ugent.be                                            *
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

#ifndef PIM_H
#define PIM_H

#include "base.h"

// ============================================================================
// PIM SINGLE/DOUBLE PRECISION FUNCTION PROTOTYPES
// ============================================================================

#ifdef USE_FLOAT
        #define pimsetpar_f77 pimssetpar_f77
        #define init_f77 cinit_f77
        #define pimtfqmr2_f77 pimctfqmr2_f77
        #define pimbicgstab_f77 pimcbicgstab_f77
        #define pimrgmres_f77 pimcrgmres_f77
#else
        #define pimsetpar_f77 pimdsetpar_f77
        #define init_f77 zinit_f77
        #define pimtfqmr2_f77 pimztfqmr2_f77
        #define pimbicgstab_f77 pimzbicgstab_f77
        #define pimrgmres_f77 pimzrgmres_f77
#endif

// ============================================================================
// PIM DOUBLE PRECISION LIBRARY FUNCTION PROTOTYPES
// ============================================================================

#define pimdsetpar_f77 F77_FUNC (pimdsetpar, PIMDSETPAR)
#define zinit_f77 F77_FUNC (zinit, ZINIT)
#define pimztfqmr2_f77 F77_FUNC (pimztfqmr2, PIMZTFQMR2)
#define pimzbicgstab_f77 F77_FUNC (pimzbicgstab, PIMZBICGSTAB)
#define pimzrgmres_f77 F77_FUNC (pimzrgmres, PIMZRGMRES)

// set parameters for iterative solving
extern "C" void pimdsetpar_f77(int *ipar, double *dpar, int *lda, int *n,
                               int *blksz, int *loclen, int *basisdim,
                               int *nprocs, int *procid, int *precontype,
                               int *stoptype, int *maxit, double *epsilon);
// initialises a vector of length n with value alpha
extern "C" void zinit_f77(int *n, dcomplex *alpha, dcomplex *x, int *incx);
// iterative solving using tfqmr
extern "C" void pimztfqmr2_f77(dcomplex *X, const dcomplex *B, dcomplex *WRK,
                              int *IPAR, double *DPAR,
                              void (*matvec_)(dcomplex*, dcomplex*, int*),
                              void (*preconl_)(dcomplex*, dcomplex*, int *),
                              void (*preconr_)(dcomplex*, dcomplex*, int *),
                              void (*pzsum_)(int *, dcomplex *, int*),
                              void (*pdznrm)(double*, int *, const dcomplex*, int *),
                              void (*progress)(int *, double*, double *),
                              int*);
// iterative solving using bi-cgstab
extern "C" void pimzbicgstab_f77(dcomplex *X, const dcomplex *B, dcomplex *WRK,
                                 int *IPAR, double *DPAR,
                                 void (*matvec_)(dcomplex *, dcomplex *, int *),
                                 void (*preconl)(dcomplex *, dcomplex *, int *),
                                 void (*preconr)(dcomplex *, dcomplex *, int *),
                                 void (*pzsum_)(int *, dcomplex *, int *),
                                 void (*pdznrm)(double *, int *, const dcomplex *,
                                 int *), void (*progress)(int *, int *, double *,
                                 dcomplex *, dcomplex *, dcomplex *));
// iterative solving using restarted-gmres
extern "C" void pimzrgmres_f77(dcomplex *X, const dcomplex *B, dcomplex *WRK,
                               int *IPAR, double *DPAR,
                               void (* matvec_)(dcomplex *, dcomplex *, int *),
                               void (*preconl_)(dcomplex *, dcomplex *, int *),
                               void (*preconr_)(dcomplex *, dcomplex *, int *),
                               void (*pzsum_)(int *, dcomplex *, int *),
                               void (*pdznrm)(double*, int *, const dcomplex *,
                               int *), void (*progress)(int *, int *, double *,
                               dcomplex *, dcomplex *, dcomplex *));

// ============================================================================
// PIM SINGLE PRECISION LIBRARY FUNCTION PROTOTYPES
// ============================================================================

#define pimssetpar_f77 F77_FUNC (pimssetpar, PIMSSETPAR)
#define cinit_f77 F77_FUNC (cinit, CINIT)
#define pimctfqmr2_f77 F77_FUNC (pimctfqmr2, PIMCTFQMR2)
#define pimcbicgstab_f77 F77_FUNC (pimcbicgstab, PIMCBICGSTAB)
#define pimcrgmres_f77 F77_FUNC (pimcrgmres, PIMCRGMRES)

// set parameters for iterative solving
extern "C" void pimssetpar_f77(int *ipar, float *dpar, int *lda, int *n,
                               int *blksz, int *loclen, int *basisdim,
                               int *nprocs, int *procid, int *precontype,
                               int *stoptype, int *maxit, float *epsilon);
// initialises a vector of length n with value alpha
extern "C" void cinit_f77(int *n, scomplex *alpha, scomplex *x, int *incx);
// iterative solving using tfqmr
extern "C" void pimctfqmr2_f77(scomplex *X, const scomplex *B, scomplex *WRK,
                              int *IPAR, float *DPAR,
                              void (*matvec_)(scomplex *, scomplex *, int *),
                              void (*preconl_)(scomplex *, scomplex *, int *),
                              void (*preconr_)(scomplex *, scomplex *, int *),
                              void (*pzsum_)(int *, scomplex *, int *),
                              void (*pdznrm)(float*, int *, const scomplex *, int *),
                              void (*progress)(int *, float *, float *),
                              int*);
// iterative solving using bi-cgstab
extern "C" void pimcbicgstab_f77(scomplex *X, const scomplex *B, scomplex *WRK,
                                 int *IPAR, float *DPAR,
                                 void (*matvec_)(scomplex *, scomplex *, int *),
                                 void (*preconl)(scomplex *, scomplex *, int *),
                                 void (*preconr)(scomplex *, scomplex *, int *),
                                 void (*pzsum_)(int *, scomplex *, int *),
                                 void (*pdznrm)(float *, int *, const scomplex *,
                                 int *), void (*progress)(int *, int *, float *,
                                 scomplex *, scomplex *, scomplex *));
// iterative solving using restarted-gmres
extern "C" void pimcrgmres_f77(scomplex *X, const scomplex *B, scomplex *WRK,
                               int *IPAR, float *DPAR,
                               void (* matvec_)(scomplex *, scomplex *, int *),
                               void (*preconl_)(scomplex *, scomplex *, int *),
                               void (*preconr_)(scomplex *, scomplex *, int *),
                               void (*pzsum_)(int *, scomplex *, int *),
                               void (*pdznrm)(float*, int *, const scomplex *,
                               int *), void (*progress)(int *, int *, float *,
                               scomplex *, scomplex *, scomplex *));

#endif
