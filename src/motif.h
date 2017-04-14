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
        std::string name;               // name of the motif
        std::vector<std::array<float, 4> > motif;               // PWM

public:
        Motif(const std::string& name_) :
                name(name_) {};

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

        /**
         * Set the name of the motif
         * @param name_ The name of the motif
         */
        void setName(const std::string& name_) {
                name = name_;
        }

        /**
         * Get the name of the motif
         * @return The name of the motif
         */
        std::string getName() const {
                return name;
        }

        /**
         * Reverse complement the motif
         */
        void revCompl();

        /**
         * Get the maximum score for this motif
         * @return The maximum score for this motif
         */
        float getMaxScore() const;

        /**
         * Get the minimum score for this motif
         * @return The minimum score for this motif
         */
        float getMinScore() const;

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
        std::vector<float> threshold;

public:
        /**
         * Create a motif file from disk
         * @param filename Filename of the motif input file
         * @param bgFreq Background frequencies
         */
        MotifContainer(const std::string& filename,
                       const std::array<float, 5>& bgFreq);

        /**
         * @param P Pre-allocated matrix P (output)
         * @param row2MotifID Index to relate a row of P to a motif ID
         * @param revCompl Flag to indicate to also add reverse complements
         */
        void generateMatrix(Matrix<float>& P, std::vector<size_t>& row2MotifID,
                            bool revCompl);

        /**
         * Set a global absolute score threshold for each motif
         * @param absThreshold Absolute score threshold
         */
        void setAbsThreshold(float absThreshold) {
                threshold = std::vector<float>(motifs.size(), absThreshold);
        }

        /**
         * Set a relative score threshold [0..1] for each motif
         * @param relThreshold Relative score threshold
         */
        void setRelThreshold(float relThreshold);

        /**
         * Write the motif names
         * @param filename File name of the motif file
         */
        void writeMotifNames(const std::string& filename);

        /**
         * Write the possum file
         * @param filename File name of the possum file
         */
        void writePossumFile(const std::string& filename);

        /**
         * Get the number of motifs
         */
        size_t size() const {
                return motifs.size();
        }

        const Motif& operator[](size_t index) const {
                return motifs[index];
        }

        float getThreshold(size_t index) const {
                return threshold[index];
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
         * Constructor with arguments
         * @param motifID_ The motif identifier
         * @param sequenceID_ The sequence identifier
         * @param sequencePos_ The sequence position identifier
         * @param score_ The motif score
         */
        MotifOccurrence(size_t motifID_, size_t sequenceID_,
                        size_t sequencePos_, float score_);

        /**
         * Get the motif identifier
         * @return The motif identifier
         */
        size_t getMotifID() const {
                return motifID;
        }

        /**
         * Get the sequence identifier
         * @return The sequence identifier
         */
        size_t getSequenceID() const {
                return sequenceID;
        }

        /**
         * Get the sequence position
         * @return The sequence position
         */
        size_t getSequencePos() const {
                return sequencePos;
        }

        /**
         * Get the motif score
         * @return The motif score
         */
        float getScore() const {
                return score;
        }

        /**
         * operator<< overloading
         * @param os Output stream
         * @param m Motif occurrence
         */
        friend std::ostream& operator<< (std::ostream& os, const MotifOccurrence& m);
};

#endif
