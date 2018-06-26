/***************************************************************************
 *   Copyright (C) 2017-2018 Jan Fostier (jan.fostier@ugent.be)            *
 *   This file is part of Blamm                                            *
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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <cstdlib>

// ============================================================================
// SETTINGS CLASS
// ============================================================================

class Settings {
public:
        size_t matrix_S_w;                      // no of rows in matrix S
        size_t matrix_S_h;                      // no of cols in matrix S
        size_t matrix_P_tile_min_size;          // minimum size of matrix P tile
        size_t matrix_P_tile_min_area;          // minimum area of matrix P tile
        double matrix_P_tile_min_zero_frac;     // minimum zero fraction in tile

public:
        /**
         * Default constructor
         */
        Settings();
};

#endif