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

#ifndef PWMSCAN_H
#define PWMSCAN_H

#include <array>
#include <vector>

#include "matrix.h"
#include "motif.h"
#include "sequence.h"
#include "settings.h"
#include "species.h"

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
        Settings settings;              // settings object
        MotifContainer motifContainer;  // motif container
        SpeciesContainer speciesContainer;      // species container

        bool simpleMode;                // naive algorithm (non-BLAS)
        bool cudaMode;                  // CUDA version

        std::mutex myMutex;             // mutex object
        size_t totMatches;              // total number of matches
        std::string outputFilename;     // occurrences filename
        std::ofstream os;               // output stream
        std::vector<MotifOccurrence> buffer;    // temporary buffer

        bool absThSpecified;            // absolute PWM threshold specified
        float absThreshold;             // asbolute PWM threshold

        bool relThSpecified;            // relative PWM threshold specified
        float relThreshold;             // relative PWM threshold

        bool pvalueSpecified;           // p-value threshold specified
        float pvalue;                   // p-value threshold

        std::string histdir;            // histogram directory

        size_t numThreads;              // number of threads used
        bool revCompl;                  // scan both directions (fwd/bwd)

        /**
         * Print module instructions
         */
        void printUsage() const;

        /**
         * Write the occurrences to disk
         * @param occurrences Occurrences to write (const-ref)
         */
        void writeOccToDisk(const std::vector<MotifOccurrence>& occurrences);

        /**
         * Write the occurrences to disk
         * @param occurrences Occurrences to write (COPY)
         */
        void writeOccToDiskCopy(const std::vector<MotifOccurrence> occurrences) {
                writeOccToDisk(occurrences);
        }

        /**
         * Given a result matrix, extract the PWM occurrences
         * @param R Result matrix R
         * @param offset Offset used to generate R = sub(S) * P
         * @param speciesID Species identifier
         * @param sm Sequence matrix
         * @param motifOcc Vector containing all motif occurrences (output)
         */
        void extractOccurrences(const Matrix& R, size_t offset,
                                size_t speciesID, SeqMatrix& sm,
                                std::vector<MotifOccurrence>& motifOcc);

        /**
         * Given a compacted occurrence stream, extract the PWM occurrences
         * @param LDR Leading dimension of matrix R
         * @param offset_v Vector with offsets
         * @param occIdx Occurrence index (in R)
         * @param occScore Occurrence PWM occScore
         * @param speciesID Species identifier
         * @param sm Sequence matrix
         * @param motifOcc Vector containing all motif occurrences (output)
         */
        void extractOccurrences2(int LDR, const std::map<int, int>& offset_v,
                                 int *occIdx, float *occScore, size_t speciesID,
                                 SeqMatrix& sm, std::vector<MotifOccurrence>& motifOcc);

        /**
         * Thread function that does the actual scanning (BLAS algorithm)
         * @param speciesID Species identifier
         * @param seqBatch Fasta sequence batch
         */
        void scanThreadBLAS(size_t speciesID,
                            const MotifContainer& motifContainer,
                            FastaBatch& seqBatch);

        /**
         * Thread function that does the actual scanning (naive algorithm)
         * @param speciesID Species identifier
         * @param seqBatch Fasta sequence batch
         */
        void scanThreadNaive(size_t speciesID,
                             FastaBatch& seqBatch);

        /**
         * Scan sequences for PWM occurrences using BLAS
         * @param seqBatch Fasta sequence batch
         */
        void scanPWMBLAS(size_t speciesID,
                         FastaBatch& seqBatch);

        /**
         * Scan sequences for PWM occurrences using the naive algorithm
         * @param seqBatch Fasta sequence batch
         */
        void scanPWMNaive(size_t speciesID,
                          FastaBatch& seqBatch);

#ifdef HAVE_CUDA
        /**
         * Thread function that does the actual scanning (CUBLAS algorithm)
         * @param devID Device ID
         * @param speciesID Species identifier
         * @param seqBatch Fasta sequence batch
         */
        void scanThreadCUBLAS(int devID, size_t speciesID,
                              FastaBatch& seqBatch);

        /**
         * Scan sequences for PWM occurrences using CUBLAS
         * @param seqBatch Fasta sequence batch
         */
        void scanPWMCUBLAS(size_t speciesID,
                           FastaBatch& seqBatch);
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
