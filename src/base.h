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

#ifndef BASEMATHPACK_H
#define BASEMATHPACK_H

#ifdef HAVE_CONFIG_H
        #include <config.h>
#endif

#define USE_FLOAT

#include <complex>

#ifdef USE_FLOAT
        typedef float deci;
        #define ZERO_THRESHOLD (1e-5)
#else
        typedef double deci;
        #define ZERO_THRESHOLD (1e-12)
#endif

typedef std::complex<double> dcomplex;
typedef std::complex<float> scomplex;
typedef std::complex<deci> cplx;

#endif
