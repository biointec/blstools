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
#include "sequence.h"

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class SpeciesContainer;

// ============================================================================
// PWMSCAN
// ============================================================================

class PWMScan
{
private:
        bool simpleMode;
	bool cudaMode;

        std::string outputFilename;

        bool absThSpecified;
        float absThreshold;

        bool relThSpecified;
        float relThreshold;

        bool pvalueSpecified;
        float pvalue;

        std::string histdir;    // histogram directory

        size_t numThreads;

        bool revCompl;

        size_t pseudocounts;


        size_t m;               // number of motifs
        size_t n;               // choose freely
        size_t k;               // maximum motif length
        size_t overlap;         // overlap
        size_t K;               // choose freely

        std::mutex myMutex;
        size_t totMatches;

        std::condition_variable cvEmpty;        // condition variable to signal buffer is empty
        std::condition_variable cvFull;         // condition variable to signal buffer is full

        bool active;                    // flag (thread is active or not)

        std::vector<MotifOccurrence> buffer;

        /**
         * Print module instructions
         */
        void printUsage() const;

        /**
         * Given a result matrix, extract the PWM occurrences
         * @param R Result matrix R
         * @param offset Offset used to generate R = P * sub(S)
         * @param sm Sequence matrix
         * @param motifContainer Motif container
         * @param motifOcc Vector containing all motif occurrences
         */
        void extractOccurrences(const Matrix<float>& R, size_t offset,
                                SeqMatrix& sm,
                                const MotifContainer& motifContainer,
                                std::vector<MotifOccurrence>& motifOcc);

        /**
         * Thread function that does the actual scanning (BLAS algorithm)
         * @param speciesID Species identifier
         * @param motifContainer Motif MotifContainer
         * @param seqBatch Fasta sequence batch
         */
        void scanThreadBLAS(size_t speciesID, const MotifContainer& motifContainer,
                        FastaBatch& seqBatch);

        /**
         * Thread function that does the actual scanning (naive algorithm)
         * @param speciesID Species identifier
         * @param motifContainer Motif MotifContainer
         * @param seqBatch Fasta sequence batch
         */
        void scanThreadNaive(size_t speciesID, const MotifContainer& motifContainer,
                             FastaBatch& seqBatch);

        /**
         * Scan sequences for PWM occurrences using BLAS
         * @param motifContainer Motif MotifContainer
         * @param seqBatch Fasta sequence batch
         */
        void scanPWMBLAS(size_t speciesID, const MotifContainer& motifContainer,
                         FastaBatch& seqBatch);

        /**
         * Scan sequences for PWM occurrences using the naive algorithm
         * @param motifContainer Motif MotifContainer
         * @param seqBatch Fasta sequence batch
         */
        void scanPWMNaive(size_t speciesID, const MotifContainer& motifContainer,
                          FastaBatch& seqBatch);

#ifdef HAVE_CUDA
        /**
         * Thread function that does the actual scanning (CUBLAS algorithm)
         * @param devID Device ID
         * @param speciesID Species identifier
         * @param motifContainer Motif MotifContainer
         * @param seqBatch Fasta sequence batch
         */
        void scanThreadCUBLAS(int devID, size_t speciesID,
                              const MotifContainer& motifContainer,
                              FastaBatch& seqBatch);

        /**
         * Scan sequences for PWM occurrences using CUBLAS
         * @param motifContainer Motif MotifContainer
         * @param seqBatch Fasta sequence batch
         */
        void scanPWMCUBLAS(size_t speciesID, const MotifContainer& motifContainer,
                           FastaBatch& seqBatch);
#endif
        /**
         * Commit some occurrences onto the output thread
         * @param chunk A number of occurrences
         */
        void commitOccurrences(const std::vector<MotifOccurrence>& chunk);

        /**
         * Entry point for the output thread
         * @param filename File name of the output file
         * @param sc Species container
         */
        void outputThread(const std::string& filename,
                          const SpeciesContainer& sc,
                          const MotifContainer& mc);

public:
        /**
         * Constructor (run scan module)
         * @param argc Command line argument count
         * @param argv Command line argument values
         */
        PWMScan(int argc, char **argv);
};

#endif
