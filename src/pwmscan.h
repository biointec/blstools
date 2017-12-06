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
         * @param os Output stream to write hits to
         */
        void scanThreadBLAS(size_t speciesID, const MotifContainer& motifContainer,
                        FastaBatch& seqBatch, std::ostream& os);

        /**
         * Thread function that does the actual scanning (naive algorithm)
         * @param speciesID Species identifier
         * @param motifContainer Motif MotifContainer
         * @param seqBatch Fasta sequence batch
         * @param os Output stream to write hits to
         */
        void scanThreadNaive(size_t speciesID, const MotifContainer& motifContainer,
                             FastaBatch& seqBatch, std::ostream& os);

        /**
         * Scan sequences for PWM occurrences using BLAS
         * @param motifContainer Motif MotifContainer
         * @param seqBatch Fasta sequence batch
         * @param os Output stream to write hits to
         */
        void scanPWMBLAS(size_t speciesID, const MotifContainer& motifContainer,
                         FastaBatch& seqBatch, std::ostream& os);

        /**
         * Scan sequences for PWM occurrences using the naive algorithm
         * @param motifContainer Motif MotifContainer
         * @param seqBatch Fasta sequence batch
         * @param os Output stream to write hits to
         */
        void scanPWMNaive(size_t speciesID, const MotifContainer& motifContainer,
                          FastaBatch& seqBatch, std::ostream& os);

#ifdef HAVE_CUDA
        /**
         * Thread function that does the actual scanning (CUBLAS algorithm)
         * @param devID Device ID
         * @param speciesID Species identifier
         * @param motifContainer Motif MotifContainer
         * @param seqBatch Fasta sequence batch
         * @param os Output stream to write hits to
         */
        void scanThreadCUBLAS(int devID, size_t speciesID,
                              const MotifContainer& motifContainer,
                              FastaBatch& seqBatch, std::ostream& os);

        /**
         * Scan sequences for PWM occurrences using CUBLAS
         * @param motifContainer Motif MotifContainer
         * @param seqBatch Fasta sequence batch
         * @param os Output stream to write hits to
         */
        void scanPWMCUBLAS(size_t speciesID, const MotifContainer& motifContainer,
                           FastaBatch& seqBatch, std::ostream& os);
#endif

public:
        /**
         * Constructor (run scan module)
         * @param argc Command line argument count
         * @param argv Command line argument values
         */
        PWMScan(int argc, char **argv);
};

#endif
