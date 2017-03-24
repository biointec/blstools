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

#ifndef MOTIF_H
#define MOTIF_H

#include <vector>
#include <array>

#include "matrix.h"

// ============================================================================
// MOTIF
// ============================================================================

class Motif {
private:
        std::string name;
        std::vector<std::array<float, 4> > motif;

public:
        Motif(const std::string& name_) : name(name_) {};

        /**
         * Add a character to the current motif
         * @param c character to add
         */
        void addCharacter(const std::array<float, 4>& c) {
                motif.push_back(c);
        }

        void posFreq2PWM(const std::array<float, 4>& bgFreq);

        /**
         * Get the length of a motif
         * @return The length of a motif
         */
        size_t size() const {
                return motif.size();
        }

        std::array<float,4> operator[](size_t i) const {
                return motif[i];
        }

        std::string getName() const {
                return name;
        }

        /**
         * Reverse complement the motif
         */
        void revCompl();

        /**
         * Compute the PWM score given a pattern
         * @param pattern Pattern to compare against
         */
        float getScore(const std::string& pattern) const;

        /**
         * operator<< overloading
         * @param os Output stream
         * @param m Motif
         */
        friend std::ostream& operator<< (std::ostream& os, const Motif& m);
};

// ============================================================================
// MOTIF CONTAINER
// ============================================================================

class MotifContainer {
private:
        std::vector<Motif> motifs;

public:
        /**
         * Create a motif file from disk
         * @param filename File name
         * @param bfFreq Background frequencies
         */
        MotifContainer(const std::string& filename,
                       const std::array<float, 5>& bgFreq);

        void generateMatrix(Matrix<float>& M);

        void addReverseCompl();

        /**
         * Get the number of motifs
         */
        size_t size() const {
                return motifs.size();
        }

        Motif operator[](size_t index) const {
                return motifs[index];
        }

        size_t getMaxMotifLen() const;
};

// ============================================================================
// MOTIF OCCURRENCES
// ============================================================================

class MotifOccurrence {
private:
        size_t motifID;
        size_t sequenceID;
        size_t sequencePos;
        float score;

public:
        /**
         * Default constructor
         */
        MotifOccurrence() {};

        /**
         * Constructor
         *
         */
        MotifOccurrence(size_t motifID_, size_t sequenceID_,
                        size_t sequencePos_, float score_);

        size_t getMotifID() const {
                return motifID;
        }

        size_t getSequenceID() const {
                return sequenceID;
        }

        size_t getSequencePos() const {
                return sequencePos;
        }

        float getScore() const {
                return score;
        }

        /**
         * operator<< overloading
         * @param os Output stream
         * @param pt PhylogeneticTree
         */
        friend std::ostream& operator<< (std::ostream& os, const MotifOccurrence& m);
};

#endif
