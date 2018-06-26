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

#include "settings.h"

#include <string>
#include <fstream>
#include <iostream>
#include <regex>

using namespace std;

// ============================================================================
// SETTINGS CLASS
// ============================================================================

Settings::Settings() : matrix_S_w(250), matrix_S_h(1000),
        matrix_P_tile_min_size(32), matrix_P_tile_min_area(128*128),
        matrix_P_tile_min_zero_frac(0.1)
{
        ifstream ifs("settings.cnf");
        if (!ifs) {
                cout << "File settings.cnf not found, using default values" << endl;
                return;
        }

        cout << "Loading configuration from file settings.cnf" << endl;
        while (ifs) {
                string line;
                getline(ifs, line);

                // remove leading spaces
                line = regex_replace(line, regex("^ +"), "");

                // skip empty lines
                if (line.empty())
                        continue;

                // skip comment lines
                if (line[0] == '#')
                        continue;

                istringstream iss(line);
                string key;
                iss >> key;

                if (key == "MATRIX_S_W") {
                        iss >> matrix_S_w;
                } else if (key == "MATRIX_S_H") {
                        iss >> matrix_S_h;
                } else if (key == "MATRIX_P_TILE_MIN_SIZE") {
                        iss >> matrix_P_tile_min_size;
                } else if (key == "MATRIX_P_TILE_MIN_AREA") {
                        iss >> matrix_P_tile_min_area;
                } else if (key == "MATRIX_P_TILE_MIN_ZERO_FRAC") {
                        iss >> matrix_P_tile_min_zero_frac;
                } else {
                        cerr << "WARNING: settings.cnf contains unknown key: " << key << endl;
                }
        }
}
