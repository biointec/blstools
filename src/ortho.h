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

#ifndef ORTHO_H
#define ORTHO_H

#include <string>
#include <map>
#include <set>
#include <vector>

class OrthoGroup
{
private:
        std::set<std::string> genes;

public:
        OrthoGroup() {}

        void insert(const std::string& genename) {
                genes.insert(genename);
        }

        size_t size() const {
                return genes.size();
        }
};

class Ortho
{
private:
        /**
         * Print module instructions
         */
        void printUsage() const;

        std::map<std::string, OrthoGroup> orthoGroups;
        std::multimap<std::string, std::string> seq2ortho;
        std::vector<std::string> seqIndex;
        std::vector<std::string> motifIndex;

public:
        /**
         * Constructor (run Ortho module)
         * @param argc Command line argument count
         * @param argv Command line argument values
         */
        Ortho(int argc, char **argv);
};

#endif
