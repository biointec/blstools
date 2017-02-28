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

#ifndef AMOS_H
#define AMOS_H

// ============================================================================
// AMOS SINGLE/DOUBLE PRECISION FUNCTION PROTOTYPES
// ============================================================================

#ifdef USE_FLOAT
        #define besj_f77 cbesj_f77
        #define besy_f77 cbesy_f77
        #define besh_f77 cbesh_f77
#else
        #define besj_f77 zbesj_f77
        #define besy_f77 zbesy_f77
        #define besh_f77 zbesh_f77
#endif

// ============================================================================
// AMOS DOUBLE PRECISION LIBRARY FUNCTION PROTOTYPES
// ============================================================================

#define zbesj_f77 F77_FUNC (zbesj, ZBESJ)
#define zbesy_f77 F77_FUNC (zbesy, ZBESY)
#define zbesh_f77 F77_FUNC (zbesh, ZBESH)

extern "C" void zbesh_f77(double *zr, double *zi, double *fnu, int *kode,
                          int *m, int *n, double *cyr, double *cyi, int *nz,
                          int *ierr);
extern "C" void zbesj_f77(double *zr, double *zi, double *fnu, int *kode,
                          int *NN, double *cyr, double *cyi, int *nz,
                          int *ierr);
extern "C" void zbesy_f77(double *zr, double *zi, double *fnu, int *kode,
                          int *NN, double *cyr, double *cyi, int *nz,
                          double *work1, double *work2, int *ierr);
extern "C" void zbesk_f77(double *zr, double *zi, double *fnu, int *kode,
                          int *n, double *cyr, double *cyi, int *nz, int *ierr);

// ============================================================================
// AMOS SINGLE PRECISION LIBRARY FUNCTION PROTOTYPES
// ============================================================================

#define cbesj_f77 F77_FUNC (cbesj, CBESJ)
#define cbesy_f77 F77_FUNC (cbesy, CBESY)
#define cbesh_f77 F77_FUNC (cbesh, CBESH)

extern "C" void cbesh_f77(float *zr, float *zi, float *fnu, int *kode,
                          int *m, int *n, float *cyr, float *cyi, int *nz,
                          int *ierr);
extern "C" void cbesj_f77(float *zr, float *zi, float *fnu, int *kode,
                          int *NN, float *cyr, float *cyi, int *nz,
                          int *ierr);
extern "C" void cbesy_f77(float *zr, float *zi, float *fnu, int *kode,
                          int *NN, float *cyr, float *cyi, int *nz,
                          float *work1, float *work2, int *ierr);
extern "C" void cbesk_f77(float *zr, float *zi, float *fnu, int *kode,
                          int *n, float *cyr, float *cyi, int *nz, int *ierr);

#endif
