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
#include <atomic>

#include "matrix.h"

// ============================================================================
// SCORE HISTOGRAM
// ============================================================================

class ScoreHistogram
{
private:
        std::atomic<size_t>* counts;          // score counts
        float minScore;
        float maxScore;
        float width;
        size_t numBins;

public:
        /**
         * Default constructor
         */
        ScoreHistogram() : counts(NULL), minScore(0.0), maxScore(0.0),
                width(0.0), numBins(0) {};

        /**
         * Default constructor
         * @param minScore Minimum score for the histogram
         * @param maxScore Maximum score for the histogram
         * @param numBins Number of bins in the histogram
         */
        ScoreHistogram(float minScore, float maxScore, size_t numBins) :
                minScore(minScore), maxScore(maxScore), numBins(numBins)
        {
                width = (maxScore - minScore) / (float)numBins;
                counts = new std::atomic<size_t>[numBins];
                for (int i = 0; i < numBins; i++)
                        counts[i] = 0;
        }

        /**
         * Copy constructor
         * @param S Histogram to copy
         */
        ScoreHistogram(const ScoreHistogram& S) {
                minScore = S.minScore;
                maxScore = S.maxScore;
                width = S.width;
                numBins = S.numBins;

                counts = new std::atomic<size_t>[numBins];
                for (int i = 0; i < numBins; i++)
                        counts[i].store(S.counts[i].load());
        }

        /**
         * Destructor
         */
        ~ScoreHistogram() {
                if (counts != NULL)
                        delete [] counts;
        }

        void addObservation(float score) {
                int binIdx = int((score - minScore) / width);
                binIdx = std::max<int>(0, binIdx);
                binIdx = std::min<int>(numBins-1, binIdx);

                counts[binIdx]++;
        }

        void print() const {
                for (int i = 0; i < numBins; i++) {
                      //  std::cout << i*width + minScore << "\t";
                       // std::cout << (i+1)*width + minScore << "\t";
                        std::cout << counts[i] << std::endl;
                }
        }

        /**
         * Compute the average of the histogram
         */
        float getAverage() const;

        /**
         * Compute the score cutoff corresponding to a certain pvalue
         * @param pvalue p-value
         * @return The score cutoff
         */
        float getScoreCutoff(float pvalue) const;

        /**
         * Write a GNUplot file containing the histogram
         * @param dir Output directory (must end by '/')
         * @param baseFilename Base filename (.dat and .gnu will be added)
         * @param label Histogram label string
         */
        void writeGNUPlotFile(const std::string& dir,
                              const std::string& baseFilename,
                              const std::string& label) const;

        /**
         * Load a histogram from disk
         * @param dir Output directory (must end by '/')
         * @param baseFilename Base filename (.dat will be added)
         */
        void loadHistogram(const std::string& dir,
                           const std::string& baseFilename);
};


// ============================================================================
// MOTIF
// ============================================================================

class Motif {
private:
        std::string name;                               // name of the motif
        size_t motifID;                                 // motif identifier
        std::vector<std::array<size_t, 4> > PFM;        // position frequency matrix
        std::vector<std::array<float, 4> > PWM;         // position weight matrix
        float threshold;                                // occurrence threshold cutoff
        bool revComp;                                   // is this a rev-compl motif

public:
        /**
         * Default constructor
         */
        Motif() : threshold(0.0), revComp(false) {}

        /**
         * Default constructor
         * @param name Name of the motif
         */
        Motif(const std::string& name, size_t motifID) : name(name),
                motifID(motifID), threshold(0.0), revComp(false) {};

        /**
         * Get the motif identifier
         * return The motif identifier
         */
        size_t getID() const {
                return motifID;
        }

        /**
         * Is the motif a reverse complement?
         * @return True or false
         */
        bool isRevCompl() const {
                return revComp;
        }

        /**
         * Add a character to the positoin frequency matrix (PFM)
         * @param c counts of the character to add
         */
        void addCharacter(const std::array<size_t, 4>& c) {
                PFM.push_back(c);
        }

        /**
         * Convert a PFM into a PWM given background nucleotide counts
         * @param bgCounts Background nucleotide counts (ACTG)
         * @param pseudoCounts Speudo counts for the *motif* PFM
         */
        void PFM2PWM(const std::array<size_t, 4>& bgCounts,
                     size_t pseudoCounts = 0);

        /**
         * Get the length of a motif
         * @return The length of a motif
         */
        size_t size() const {
                return PFM.size();
        }

        /**
         * Operator[] overloading to access the PWM
         * @param i position index
         */
        const std::array<float,4>& operator[](size_t i) const {
                return PWM[i];
        }

        /**
         * Get the name of the motif
         * @return The name of the motif
         */
        const std::string& getName() const {
                return name;
        }

        /**
         * Get the base name of the motif (with __XX for random permutations)
         * @return The name of the motif
         */
        std::string getBaseName() const {
                std::string retval = name;
                ssize_t i = name.size()-1;
                for ( ; i >= 0; i--)
                        if (!isdigit(name[i]))
                                break;
                if (i < 1)
                        return retval;
                if ((name[i-1] != '_') || (name[i] != '_'))
                        return retval;
                return retval.substr(0, i-1);
        }

        /**
         * Check whether a motif is a permutation (__XX at the end)
         * @return true or false
         */
        bool isPermutation() const {
                ssize_t i = name.size()-1;
                for ( ; i >= 0; i--)
                        if (!isdigit(name[i]))
                                break;
                if (i < 1)
                        return false;
                return ((name[i-1] == '_') && (name[i] == '_'));
        }

        /**
         * Reverse complement the motif (both PFM and PWM)
         */
        void revCompl();

        /**
         * Write the MOODS file
         * @param filename Filename
         */
        void writeMOODSFile(const std::string& filename) const;

        /**
         * Get the maximum PWM score for this motif
         * @return The maximum PWM score for this motif
         */
        float getMaxScore() const;

        /**
         * Get the minimum PWM score for this motif
         * @return The minimum PWM score for this motif
         */
        float getMinScore() const;

        /**
         * Compute the PWM score given a pattern
         * @param pattern Pattern to compare against
         */
        float getScore(const std::string& pattern) const;

        /**
         * Set the motif cutoff threshold
         * @param threshold_ Motif score threshold
         */
        void setThreshold(float threshold_) {
                threshold = threshold_;
        }

        /**
         * Get the motif cutoff threshold
         * @return Motif score threshold
         */
        float getThreshold() const {
                return threshold;
        }

        bool operator<(const Motif& rhs) const {
                return size() < rhs.size();
        }

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

        Matrix<float> P;
        std::vector<std::pair<size_t, size_t> > matBlock;

        std::vector<size_t> row2MotifID;

public:
        /**
         * Create a motif file from disk
         * @param filename Filename of the motif input file
         * @param loadPermutations Also load the random permutations
         */
        MotifContainer(const std::string& filename, bool loadPermutations);

        /**
         * Load cluster-buster motifs from disk
         * @param filename Filename of the motif input file
         * @param motifs Vector to store motifs (output)
         */
        void loadCBMotifs(const std::string& filename,
                          std::vector<Motif>& motifs);

        /**
         * Load Jaspar motifs from disk
         * @param filename Filename of the motif input file
         * @param motifs Vector to store motifs (output)
         */
        void loadJasparMotifs(const std::string& filename,
                              std::vector<Motif>& motifs);

        /**
         * Add the reverse complement motifs to the container
         */
        void addReverseComplements();

        /**
         * Generate a motif (pattern) matrix P
         */
        void generateMatrix();

        /**
         * Get a const-reference to matrix P
         * @return A const-reference to matrix P
         */
        const Matrix<float>& getMatrix() const {
                return P;
        }

        /**
         * Get a const-reference to a matrix block
         * @return A const-reference to a matrix block
         */
        const std::vector<std::pair<size_t, size_t> >& getMatrixBlock() const {
                return matBlock;
        }

        /**
         * Get the motif ID contained in row i of the matrix
         * @return The motif ID at row i of the matrix
         */
        size_t getMotifIDAtRow(size_t i) const {
                return row2MotifID[i];
        }

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
         * Write the MOODS file
         */
        void writeMOODSFiles();

        /**
         * Get the number of motifs
         */
        size_t size() const {
                return motifs.size();
        }

        const Motif& operator[](size_t index) const {
                return motifs[index];
        }

        size_t getMaxMotifLen() const;

        /**
         * Return an iterator pointing to the first motif in the container
         * @return An iterator to the beginning of the motif container
         */
        std::vector<Motif>::const_iterator begin() const {
                return motifs.begin();
        }

        /**
         * Return an iterator pointing past the final motif in the container
         * @return An iterator to the end of the motif container
         */
        std::vector<Motif>::const_iterator end() const {
                return motifs.end();
        }

        /**
         * Return an iterator pointing to the first motif in the container
         * @return An iterator to the beginning of the motif container
         */
        std::vector<Motif>::iterator begin() {
                return motifs.begin();
        }

        /**
         * Return an iterator pointing past the final motif in the container
         * @return An iterator to the end of the motif container
         */
        std::vector<Motif>::iterator end() {
                return motifs.end();
        }
};

// ============================================================================
// MOTIF OCCURRENCES
// ============================================================================

class MotifOccurrence {
private:
        size_t motifID;
        size_t sequenceID;
        size_t sequencePos;
        char strand;
        float score;

public:
        /**
         * Default constructor
         */
        MotifOccurrence() {};

        /**
         * Constructor with arguments
         * @param motifID The motif identifier
         * @param sequenceID The sequence identifier
         * @param sequencePos The sequence position identifier
         * @param strand + or - strand
         * @param score The motif score
         */
        MotifOccurrence(size_t motifID, size_t sequenceID, size_t sequencePos,
                        char strand, float score) :
                motifID(motifID), sequenceID(sequenceID),
                sequencePos(sequencePos), strand(strand), score(score) {}

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
         * Get the sequence strand
         * @return The sequence strand
         */
        char getStrand() const {
                return strand;
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
