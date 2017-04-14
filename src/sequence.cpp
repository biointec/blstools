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

#include <stdexcept>
#include <fstream>
#include <array>
#include <iostream>
#include <sstream>

#include "sequence.h"

using namespace std;

// ============================================================================
// SEQUENCE BLOCK
// ============================================================================

void SeqBlock::getSuffixBlock(SeqBlock& sb, size_t startPos)
{
        // if the current block is too small, get out
        if (startPos > size())
                return;

        sb.clear();
        sb.block = block.substr(startPos);

        auto it = block2seq.upper_bound(startPos);
        assert(it != block2seq.begin());
        --it;

        // copy the primary sequence marker
        sb.block2seq[0] = SeqPos(it->second.getSeqIndex(),
                                 it->second.getSeqPos() + startPos - it->first);

        // copy the remainder of sequence markers (if any)
        for (it++ ; it != block2seq.end(); it++)
                sb.block2seq[it->first-startPos] = it->second;
}

std::ostream& operator<< (std::ostream& os, const SeqBlock& sb)
{
        os << sb.block << endl;
        for (auto el : sb.block2seq)
                os << el.first << ": (" << el.second.getSeqIndex()
                   << ", " << el.second.getSeqPos() << ")" << endl;
        return os;
}

// ============================================================================
// FASTA BATCH HANDLER
// ============================================================================

bool FastaBatch::moveToNextFile()
{
        // close input stream if necessary
        if (_ifs.is_open())
                _ifs.close();

        // return false if no files are left
        if (_currFileIdx >= seqFiles.size())
                return false;

        // open new sequence file
        _ifs.open(seqFiles[_currFileIdx]);
        if (!_ifs)
                throw runtime_error("Could not open file: " + seqFiles[_currFileIdx]);

        _currFileIdx++;
        return true;
}

bool FastaBatch::getNextLine(std::string& line, size_t& seqIdx, size_t& seqPos)
{
        while (true) {
                if (!getline(_ifs, line))
                        if (!moveToNextFile())
                                return false;

                if (line.empty())
                        continue;

                if (line.front() == '>') {
                        _currSeqLen = 0;
                        _seqNames.push_back(line.substr(1));
                        continue;
                }

                assert(!_seqNames.empty());

                seqIdx = _seqNames.size() - 1;
                seqPos = _currSeqLen;

                _totSeqLen += line.size();
                _currSeqLen += line.size();
                return true;
        }

        return true;
}

void FastaBatch::addFile(const string& filename)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + filename);

        seqFiles.push_back(filename);
}

void FastaBatch::addList(const std::string& listname)
{
        ifstream ifs(listname.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + listname);

        // read all files from the list and check for their existance
        while (ifs) {
                string filename;
                getline(ifs, filename);
                if (filename.empty())
                        continue;

                ifstream fn(filename);
                if (!fn)
                        throw runtime_error("Could not open fasta file: " + filename);

                seqFiles.push_back(filename);
        }
}

void FastaBatch::calcBGFrequencies()
{
        for (size_t i = 0; i < 256; i++)
                char2Idx[i] = 4;
        char2Idx[(unsigned short)'A'] = 0;
        char2Idx[(unsigned short)'a'] = 0;
        char2Idx[(unsigned short)'C'] = 1;
        char2Idx[(unsigned short)'c'] = 1;
        char2Idx[(unsigned short)'G'] = 2;
        char2Idx[(unsigned short)'g'] = 2;
        char2Idx[(unsigned short)'T'] = 3;
        char2Idx[(unsigned short)'t'] = 3;

        std::array<size_t, 5> nuclCount;
        nuclCount.fill(0);

        string line; size_t dummy;
        while (getNextLine(line, dummy, dummy)) {
                for (auto c : line)
                        nuclCount[char2Idx[(unsigned short)c]]++;
        }

        for (size_t i = 0; i < 5; i++)
                nucleotideFreq[i] = (float)nuclCount[i] / (float)_totSeqLen;
}

void FastaBatch::writeSeqNames(const std::string& filename)
{
        ofstream ofs(filename.c_str());

        for (auto it : _seqNames)
                ofs << it << "\n";

        ofs.close();
}

bool FastaBatch::getNextBlock(SeqBlock& block, size_t maxSize)
{
        // first copy from the active sequence (sequence data previously read)
        size_t thisSize = min<size_t>(maxSize - block.size(), _actSeq.size() - _actPos);
        block.append(_actSeq.substr(_actPos, thisSize),
                     SeqPos(_actSeqIdx, _actSeqPos + _actPos));
        _actPos += thisSize;

        // read additional data from disk
        while (block.size() < maxSize) {
                _actPos = 0;
                if (!getNextLine(_actSeq, _actSeqIdx, _actSeqPos))
                        return !block.empty();

                thisSize = min<size_t>(maxSize - block.size(), _actSeq.size());
                block.append(_actSeq.substr(0, thisSize),
                             SeqPos(_actSeqIdx, _actSeqPos));
                _actPos = thisSize;
        }

        return !block.empty();
}

bool FastaBatch::getNextOverlappingBlock(SeqBlock& block, size_t maxSize,
                                         size_t overlap)
{
        lock_guard<mutex> lock(m);

        block = _nextBlock;

        bool retVal = getNextBlock(block, maxSize);
        block.getSuffixBlock(_nextBlock, block.size() - overlap);

        return retVal;
}

// ============================================================================
// SEQUENCE MATRIX
// ============================================================================

bool SeqMatrix::getNextSeqMatrix(FastaBatch& bf)
{
        // compute the maximum block size
        size_t maxBlockSize = K * numCol + overlap;

        // get the next block
        bf.getNextOverlappingBlock(block, maxBlockSize, overlap);

        // compute the number of non-zero columns in the sequence matrix
        numOccCol = (block.size() - overlap + K - 1) / K;

        // fill the sequence matrix
        S.fill(0.0f);
        for (size_t j = 0; j < numOccCol; j++) {
                for (size_t i = 0; i < (K + overlap); i++) {
                        if (block[j*K+i] == 'A')
                                S(4*i+0, j) = 1;
                        if (block[j*K+i] == 'C')
                                S(4*i+1, j) = 1;
                        if (block[j*K+i] == 'G')
                                S(4*i+2, j) = 1;
                        if (block[j*K+i] == 'T')
                                S(4*i+3, j) = 1;
                }
        }

        return numOccCol > 0;
}

void SeqMatrix::extractOccurrences(const Matrix<float>& R, const vector<size_t>& row2motifID,
                                   std::vector<MotifOccurrence>& motifOcc,
                                   size_t offset, const MotifContainer& motifs)
{
        for (size_t j = 0; j < numOccCol; j++) {
                for (size_t i = 0; i < R.nRows(); i++) {
                        float thisScore = R(i,j);
                        size_t motifID = row2motifID[i];

                        if (thisScore < motifs.getThreshold(motifID))
                                continue;

                        // at this point an occurrence is found
                        size_t seqID = block.getSeqIdx(j*K+offset);
                        size_t seqPos = block.getSeqPos(j*K+offset);
                        size_t remSeqLen = block.getRemainingSeqLen(j*K+offset);

                        if (motifs[motifID].size() > remSeqLen)
                                continue;

                        motifOcc.push_back(MotifOccurrence(motifID, seqID, seqPos, thisScore));
                }
        }
}

void SeqMatrix::findOccurrences(const Matrix<float>& P,
                                const vector<size_t>& row2motifID,
                                vector<MotifOccurrence>& motifOcc,
                                const MotifContainer& motifs)
{
        motifOcc.clear();

        Matrix<float> R(P.nRows(), numOccCol);
        for (int offset = 0; offset < K; offset++) {
                R.gemm(P, S, 4*offset, numOccCol);
                extractOccurrences(R, row2motifID, motifOcc, offset, motifs);
        }
}

std::ostream& operator<< (std::ostream& os, const SeqMatrix& sm)
{
        os << "Matrix dims: " << sm.S.nRows() << " x " << sm.S.nCols()
           << ", of which " << sm.numOccCol << " are occupied\n";
        os << sm.block << "\n";
        sm.S.printSequence(sm.overlap);
        return os;
}
