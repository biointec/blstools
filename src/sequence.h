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

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <mutex>

#include "matrix.h"
#include "motif.h"

// ============================================================================
// FUNCTION PROTOTYPE
// ============================================================================

class FastaBatch;

// ============================================================================
// SEQUENCE POSITION
// ============================================================================

class SeqPos {

private:
        size_t seqIdx;  // sequence index
        size_t seqPos;  // position in the sequence

public:
        /**
         * Default constructor
         */
        SeqPos() : seqIdx(0), seqPos(0) {}

        /**
         * Default constructor
         * @param seqIdx_ Sequence index
         * @param seqPos_ Sequence position
         */
        SeqPos(size_t seqIdx_, size_t seqPos_) :
                seqIdx(seqIdx_), seqPos(seqPos_) {}

        /**
         * Get the sequence index
         * @return The sequence index
         */
        size_t getSeqIndex() const {
                return seqIdx;
        }

        /**
         * Get the sequence position
         * @return The sequence position
         */
        size_t getSeqPos() const {
                return seqPos;
        }
};

// ============================================================================
// SEQUENCE BLOCK
// ============================================================================

class SeqBlock {

private:
        std::string block;              // block of one or more sequences
        std::map<size_t, SeqPos> block2seq;     // block to sequence map

public:
        /**
         * Given a block position, get the corresponding sequence position
         * @param blockPos Block position
         * @return The sequence position
         */
        size_t getSeqPos(size_t blockPos) const {
                auto it = block2seq.upper_bound(blockPos);
                assert(it != block2seq.begin());
                --it;
                return it->second.getSeqPos() + blockPos - it->first;
        }

        /**
         * Given a block position, get the remaining length within the sequence
         * @param blockPos Block position
         * @return The remaining length within this sequence
         */
        size_t getRemainingSeqLen(size_t blockPos) const {
                auto it = block2seq.upper_bound(blockPos);
                assert(it != block2seq.begin());
                size_t nextBlockPos = (it == block2seq.end()) ? block.size() : it->first;
                return nextBlockPos - blockPos;
        }

        /**
         * Given a block position, get the corresponding sequence index
         * @param blockPos Block position
         * @return The sequence index
         */
        size_t getSeqIdx(size_t blockPos) const {
                auto it = block2seq.upper_bound(blockPos);
                assert(it != block2seq.begin());
                --it;
                return it->second.getSeqIndex();
        }

        /**
         * Append a new sequence to the block
         * @param data Sequence content
         * @param sp Sequence position
         */
        void append(const std::string& data, const SeqPos& sp) {
                if ((block2seq.empty()) ||
                    (block2seq.rbegin()->second.getSeqIndex() != sp.getSeqIndex()))
                        block2seq[block.size()] = sp;
                block.append(data);
        }

        /**
         * Clear the sequence block
         */
        void clear() {
                block.clear();
                block2seq.clear();
        }

        /**
         * Is the sequence block empty()
         * @return True if empty
         */
        bool empty() const {
                return block.empty();
        }

        /**
         * Get the size of the sequence block
         * @return The size of the sequence block
         */
        size_t size() const {
                return block.size();
        }

        /**
         * Get a specific character
         * @param blockPos Block position
         */
        char operator[](size_t blockPos) {
                return block[blockPos];
        }

        /**
         * Copy the suffix of this block into a new block
         * @param newBlock New block (output)
         * @param startPos Start position of the suffix
         */
        void getSuffixBlock(SeqBlock& newBlock, size_t startPos);

        /**
         * operator<< overloading
         * @param os Output stream
         * @param pt PhylogeneticTree
         */
        friend std::ostream& operator<< (std::ostream& os, const SeqBlock& sb);
};

// ============================================================================
// FASTA BATCH HANDLER
// ============================================================================

class FastaBatch {

private:
        /**
         * Close the current input file stream and move to the next one
         * @return False if no more files are left, true otherwise
         */
        bool moveToNextFile();

        // following variables are modified as side effect of moveToNextFile()
        std::ifstream _ifs;                     // current input file stream
        size_t _currFileIdx;                    // current file index

        /**
         * Get the next line of sequence content
         * @param line String where the line is stored
         * @param seqIdx Sequence index from which line is retreived
         * @param seqPos Position of the line in the current sequence
         * @return False if no more lines are left, true otherwise
         */
        bool getNextLine(std::string& line, size_t& seqIdx, size_t& seqPos);

        // following variables are modified as side effect of getNextLine()
        std::vector<std::string> _seqNames;     // sequence names
        size_t _currSeqLen;                     // current sequence length
        size_t _totSeqLen;                      // total sequence length

        /**
         * Get a raw block of (appended) sequence data
         * @param block Block of bulk sequence data
         * @param maxSize Maximum size of the block
         * @return True if the block is non-empty
         */
        bool getNextBlock(SeqBlock& block, size_t maxSize);

        // following variables are modified as side effect of getNextBlock()
        std::string _actSeq;                    // active sequence
        size_t _actSeqIdx;                      // active sequence index
        size_t _actSeqPos;                      // active sequence position
        size_t _actPos;                         // position within active sequence

        // following variables are modified as side effect of getNextOverlappingBlock()
        SeqBlock _nextBlock;

        int char2Idx[256];
        std::vector<std::string> seqFiles;
        std::array<float, 5> nucleotideFreq;    // nucleotide frequency

        std::mutex m;

public:
        /**
         * Default constructor
         */
        FastaBatch() : _currFileIdx(0), _currSeqLen(0), _totSeqLen(0),
                       _actSeqIdx(0), _actSeqPos(0), _actPos(0) {}

        /**
         * Calculate the background frequencies
         */
        void calcBGFrequencies();

        /**
         * Write the sequence names
         * @param filename File name of the output file
         */
        void writeSeqNames(const std::string& filename);

        /**
         * Add a single sequence file to process
         * @param filename File name
         */
        void addFile(const std::string& filename);

        /**
         * Add a list of sequence files to process
         * @param listname File name of the list
         */
        void addList(const std::string& listname);

        /**
         * Reset the fasta batch handler
         */
        void reset() {
                _currFileIdx = _currSeqLen = _totSeqLen = 0;
                _actSeqIdx = _actSeqPos = _actPos = 0;
                if (_ifs.is_open())
                        _ifs.close();
                _seqNames.clear();
                _actSeq.clear();
                _nextBlock.clear();
        }

        /**
         * Get the background nucleotide frequencies
         * @return The background nucleotide frequencies
         */
        std::array<float, 5> getNucleotideFreq() const {
                return nucleotideFreq;
        }

        /**
         * Given a sequence index, get the sequence name
         * @param seqIdx Sequence index
         * @return The sequence name
         */
        std::string getSeqName(size_t seqIdx) const {
                return _seqNames[seqIdx];
        }

        /**
         * Get the number of sequence files in this batch
         * @return The number of sequence files in this batch
         */
        size_t size() const {
                return seqFiles.size();
        }

        /**
         * Get the total number of sequences in all files
         * @return The total number of sequences
         */
        size_t getNumSequences() const {
                return _seqNames.size();
        }

        /**
         * Get the total sequence length
         * @return The total sequence length
         */
        size_t getTotalSeqLength() const {
                return _totSeqLen;
        }

        /**
         * Get a raw block of (appended) sequence data
         * @param block Block of bulk sequence data
         * @param maxSize Maximum size of the block
         * @param overlap Overlap with the previous block
         * @return True if the block is non-empty
         */
        bool getNextOverlappingBlock(SeqBlock& block, size_t maxSize, size_t overlap);
};

// ============================================================================
// SEQUENCE MATRIX
// ============================================================================

class SeqMatrix {

private:
        size_t K;                       // number of sequence rows
        size_t overlap;                 // number of overlap rows
        size_t numCol;                  // number of columns
        size_t numOccCol;               // number of occupied columns

        Matrix<float> S;                // actual matrix
        SeqBlock block;                 // sequence block

        void extractOccurrences(const Matrix<float>& R,
                                std::vector<MotifOccurrence>& motifOcc,
                                size_t offset, const MotifContainer& motifs);

public:
        /**
         * Default constructor
         * @param K_ Number of sequence rows
         * @param overlap_ Number of overlap rows
         * @param numCol_ Number of sequence columns
         */
        SeqMatrix(size_t K_, size_t overlap_, size_t numCol_) :
                K(K_), overlap(overlap_), numCol(numCol_), numOccCol(0) {
                        S = Matrix<float>(4*(K+overlap), numCol);
                }

        /**
         * Fill the next sequence matrix
         * @param bf Fasta batch reader
         * @return True if more data is present
         */
        bool getNextSeqMatrix(FastaBatch& bf);

        /**
         * Find the motif findOccurrences
         * @param P Pattern matrix
         * @param motifOcc Motif occurrences (output)
         * @param motifs Motif container
         */
        void findOccurrences(const Matrix<float>& P,
                             std::vector<MotifOccurrence>& motifOcc,
                             const MotifContainer& motifs);

        /**
         * operator<< overloading
         * @param os Output stream
         * @param pt PhylogeneticTree
         */
        friend std::ostream& operator<< (std::ostream& os, const SeqMatrix& sm);
};

#endif
