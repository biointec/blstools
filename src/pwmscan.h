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

#ifndef PWMSCAN_H
#define PWMSCAN_H

#include <array>
#include <vector>

#include "matrix.h"
#include "motif.h"

// ============================================================================
// PWMSCAN
// ============================================================================

class PWMScan
{
private:
        std::string outputFilename;
        bool listProvided;
        float threshold;

        /**
         * Print module instructions
         */
        void printUsage() const;

public:
        /**
         * Constructor (run scan module)
         * @param argc Command line argument count
         * @param argv Command line argument values
         */
        PWMScan(int argc, char **argv);

        void loadFasta(const std::string& filename, std::string& sequence);

        void countFrequencies(const std::string& sequence,
                              std::array<float, 5>& frequencies);
};

#endif
