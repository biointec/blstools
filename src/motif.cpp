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

#include <fstream>
#include <sstream>
#include <cmath>

#include <iostream>
#include <algorithm>
#include <numeric>

#include "motif.h"

using namespace std;

// ============================================================================
// SCORE HISTOGRAM
// ============================================================================

float ScoreHistogram::getAverage() const
{
        double sum = 0.0, total = 0.0;
        for (size_t i = 0; i < numBins; i++) {
                sum += ((0.5 + i) * width + minScore) * counts[i];
                total += counts[i];
        }

        return sum / total;
}

float ScoreHistogram::getScoreCutoff(float pvalue) const
{
        // count the total number of observations in the histogram
        double totObs = 0.0;
        for (size_t i = 0; i < numBins; i++)
                totObs += counts[i];

        // compute the fraction of best observations
        double bestObs = pvalue * totObs;

        double curr = bestObs;
        for (ssize_t i = numBins-1; i >= 0; i--) {
                if (counts[i] < curr)
                        curr -= counts[i];
                else {
                        double frac = curr / counts[i];
                        float cutoffi = frac * i + (1.0-frac) * (i+1);
                        return cutoffi * width + minScore;
                }
        }

        return maxScore;
}

void ScoreHistogram::writeGNUPlotFile(const string& dir,
                                      const string& baseFilename,
                                      const string& label) const
{
        string filename = dir + baseFilename + ".dat";
        ofstream ofs(filename.c_str());
        if (!ofs)
                throw runtime_error("Error: cannot write to file " + filename);

        ofs << numBins << "\t" << minScore << "\t" << maxScore << "\n";
        for (size_t i = 0; i < numBins; i++)
                ofs << (0.5 + i) * width + minScore << "\t" << counts[i] << "\n";
        ofs.close();

        size_t maxy = 0;
        for (size_t i = 0; i < numBins; i++)
                maxy = max<size_t>(maxy, counts[i]);
        maxy *= 1.1;

        filename = dir + baseFilename + ".gnu";
        ofs.open(filename.c_str());
        if (!ofs)
                throw runtime_error("Error: cannot write to file " + filename);

        ofs << "set output \"" << baseFilename << ".ps\"\n";
        ofs << "set key autotitle columnhead\n";
        ofs << "set terminal postscript landscape\n";
        ofs << "set terminal postscript noenhanced\n";
        ofs << "set xrange [" << minScore << ":" << maxScore << "]\n";
        ofs << "set yrange [" << 0 << ":" << maxy << "]\n";
        ofs << "set xlabel \'PWM score\'" << endl;
        ofs << "set ylabel \'count\'" << endl;
        ofs << "plot \"" << baseFilename << ".dat\" using 1:2 title \'"
            << label << "\' with boxes\n";

        ofs.close();
}

void ScoreHistogram::loadHistogram(const std::string& dir,
                                   const std::string& baseFilename)
{
        string filename = dir + baseFilename + ".dat";
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Error: cannot read file " + filename + ". "
                                    "Did you run the hist module?");

        ifs >> numBins >> minScore >> maxScore;
        width = (maxScore - minScore) / (float)numBins;
        if (counts != NULL)
                delete [] counts;
        counts = new std::atomic<size_t>[numBins];

        for (size_t i = 0; i < numBins; i++) {
                float bin; size_t val;
                ifs >> bin >> val;
                counts[i] = val;
        }

        if (!ifs)
                throw runtime_error("Unexpected end-of-file reached");
}

// ============================================================================
// MOTIF
// ============================================================================

size_t char2idx(char c)
{
        if (c == 'A' || c == 'a')
                return 0;
        if (c == 'C' || c == 'c')
                return 1;
        if (c == 'G' || c == 'g')
                return 2;
        if (c == 'T' || c == 't')
                return 3;
        return 4;
}

void Motif::PFM2PWM(const std::array<size_t, 4>& bgCounts, size_t pseudoCounts)
{
        // compute the background probability for ACGT
        size_t bgTotCounts = accumulate(bgCounts.begin(), bgCounts.end(), 0ull);
        bgTotCounts += 4 * pseudoCounts;
        array<float, 4> bgProb;
        for (size_t i = 0; i < 4; i++)
                bgProb[i] = static_cast<double>(bgCounts[i] + pseudoCounts) / bgTotCounts;

        // if the motif is a reverse-complementary motif, also complement the bgProb
        if (revComp) {
                swap<float>(bgProb[0], bgProb[3]);
                swap<float>(bgProb[1], bgProb[2]);
        }

        // compute the PWM
        PWM.resize(PFM.size());
        for (size_t i = 0; i < PFM.size(); i++) {
                size_t totCounts = accumulate(PFM[i].begin(), PFM[i].end(), 0ull);
                totCounts += 4 * pseudoCounts;

                for (size_t j = 0; j < 4; j++) {
                        // compute the PPM
                        PWM[i][j] = static_cast<double>(PFM[i][j] + pseudoCounts) / totCounts;
                        // and convert to PWM
                        PWM[i][j] = log2(PWM[i][j] / bgProb[j]);
                }
        }
}

float Motif::getScore(const std::string& pattern) const
{
        // make sure the pattern and motif have the same size
        assert (pattern.size() == size());

        float score = 0.0f;
        for (size_t i = 0; i < pattern.size(); i++) {
                size_t j = char2idx(pattern[i]);
                if (j < 4)      // if ACTG character
                        score += PWM[i][j];
        }

        return score;
}

float Motif::getMaxScore() const
{
        float maxScore = 0.0f;

        for (auto& pos : PWM) {
                float maxAC = max<float>(pos[0], pos[1]);
                float maxGT = max<float>(pos[2], pos[3]);
                maxScore += max<float>(maxAC, maxGT);
        }

        return maxScore;
}

float Motif::getMinScore() const
{
        float minScore = 0.0f;

        for (auto& pos : PWM) {
                float minAC = min<float>(pos[0], pos[1]);
                float minGT = min<float>(pos[2], pos[3]);
                minScore += min<float>(minAC, minGT);
        }

        return minScore;
}

void Motif::revCompl()
{
        // reverse complement the position frequency matrix
        reverse(PFM.begin(), PFM.end());

        for (size_t i = 0; i < PFM.size(); i++) {
                array<size_t, 4> copy = PFM[i];
                PFM[i] = array<size_t, 4>{copy[3], copy[2], copy[1], copy[0]};
        }

        // reverse complement the position weight matrix
        reverse(PWM.begin(), PWM.end());

        for (size_t i = 0; i < PWM.size(); i++) {
                array<float, 4> copy = PWM[i];
                PWM[i] = array<float, 4>{copy[3], copy[2], copy[1], copy[0]};
        }

        revComp = !revComp;
}

void Motif::writeMOODSFile(const std::string& filename) const
{
        ofstream ofs(filename.c_str());

        for (size_t i = 0; i < 4; i++) {
                ofs << PFM[0][i];
                for (size_t j = 1; j < PFM.size(); j++)
                        ofs << "\t" << PFM[j][i];
                ofs << "\n";
        }
}

ostream& operator<< (ostream& os, const Motif& m)
{
        os << m.name << "\n";
        for (auto pos : m.PWM)
                os << pos[0] << " " << pos[1] << " " << pos[2] << " " << pos[3] << "\n";
        return os;
}

// ============================================================================
// MOTIF CONTAINER
// ============================================================================

MotifContainer::MotifContainer(const std::string& filename, bool loadPermutations)
{
        vector<Motif> allMotifs;

        if ((filename.size() > 7) && (filename.substr(filename.size() - 7) == ".jaspar"))
                loadJasparMotifs(filename, allMotifs);
        else
                loadCBMotifs(filename, allMotifs);

        // copy the temporary allmotifs structure to motifs
        for (const auto& motif : allMotifs) {
                if (loadPermutations || !motif.isPermutation())
                        motifs.push_back(motif);
        }
}

void MotifContainer::loadCBMotifs(const std::string& filename,
                                  vector<Motif>& motifs)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + filename);

        size_t motifID = 0;
        while (ifs.good()) {
                string temp;
                getline(ifs, temp);
                if (temp.empty())
                        continue;
                if (temp.front() == '>') {      // add a new motif
                        motifs.push_back(Motif(temp.substr(1), motifID++));
                        continue;
                }
                istringstream iss(temp);
                size_t A, C, G, T;
                iss >> A >> C >> G >> T;

                motifs.back().addCharacter({A, C, G, T});
        }

        sort(motifs.begin(), motifs.end());
}

void MotifContainer::loadJasparMotifs(const std::string& filename,
                                      vector<Motif>& motifs)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + filename);

        size_t motifID = 0;
        while (ifs.good()) {
                string motifName, temp;
                ifs >> motifName;
                if (!motifName.empty())
                        motifName = motifName.substr(1);
                getline(ifs, temp);
                if (!ifs)
                        break;

                vector<vector<size_t> > freq(4);
                for (size_t i = 0; i < 4; i++) {
                        getline(ifs, temp);

                        istringstream iss(temp);
                        iss >> temp;
                        iss >> temp;

                        while (iss) {
                                size_t count;
                                iss >> count;
                                if (!iss)
                                        break;

                                freq[i].push_back(count);
                        }
                }

                motifs.push_back(Motif(motifName, motifID++));
                for (size_t i = 0; i < freq[0].size(); i++)
                        motifs.back().addCharacter({freq[0][i], freq[1][i], freq[2][i], freq[3][i]});
        }

        sort(motifs.begin(), motifs.end());
}

void MotifContainer::addReverseComplements()
{
        vector<Motif> copy = motifs;
        motifs.clear();

        for (Motif& m : copy) {
                motifs.push_back(m);
                m.revCompl();
                motifs.push_back(m);
        }
}

TilePair MotifContainer::findBestSplit(const MatrixTile& input)
{
        size_t bestZeroArea = 0;
        size_t bestJ = input.colStart + 1;

        for (size_t j = input.colStart + 1; j < input.colEnd; j++) {
                size_t numCols = j - input.colStart;
                size_t numRows = input.rowEnd - 4*motifs[j-1].size();
                size_t zeroArea = numCols * numRows;

                if (zeroArea > bestZeroArea) {
                        bestZeroArea = zeroArea;
                        bestJ = j;
                }
        }

        MatrixTile left(input.rowStart, input.colStart, 4*motifs[bestJ-1].size(), bestJ);
        MatrixTile right(input.rowStart, bestJ, input.rowEnd, input.colEnd);

        return make_pair(left, right);
}

bool MotifContainer::keepSplit(const TilePair& tilePair)
{
        const size_t minSize = 32;
        const size_t minArea = 128*128;
        const double minZeroFrac = 0.10;

        const MatrixTile& f = tilePair.first;
        const MatrixTile& s = tilePair.second;

        size_t origArea = (s.rowEnd-f.rowStart) * (s.colEnd-f.colStart);
        size_t area1 = (f.rowEnd-f.rowStart) * (f.colEnd-f.colStart);
        size_t area2 = (s.rowEnd-s.rowStart) * (s.colEnd-s.colStart);

        // check the minimum matrix dimensions
        if (f.colEnd - f.colStart < minSize)
                return false;
        if (f.rowEnd - f.rowStart < minSize)
                return false;
        if (s.colEnd - s.colStart < minSize)
                return false;
        if (s.rowEnd - s.rowStart < minSize)
                return false;

        // check the minimum area
        if (area1 < minArea)
                return false;
        if (area2 < minArea)
                return false;

        // check the relative area gain
        if (double(area1 + area2)/double(origArea) > (1.0-minZeroFrac))
                return false;

        return true;
}

void MotifContainer::generateMatrix()
{
        // allocate memory for matrix P
        size_t m = motifs.size();
        size_t k = 4 * getMaxMotifLen();

        P.resize(k, m, 0.0f);

        cout << "Matrix P has dimensions: " << k << " x " << m << endl;

        // fill the matrix
        for (size_t i = 0; i < motifs.size(); i++) {
                // fill in the forward motif
                const Motif& fwd = motifs[i];
                for (size_t j = 0; j < fwd.size(); j++)
                        for (size_t o = 0; o < 4; o++)
                                P(4*j+o, i) = fwd[j][o];
                col2MotifID.push_back(i);
        }

        matrixTiles.push_back(MatrixTile(0, 0, P.nRows(), P.nCols()));
        while (true) {
                vector<MatrixTile> newTiles;
                bool didSomething = false;

                for (const auto& it : matrixTiles) {
                        auto res = findBestSplit(it);

                        if (keepSplit(res)) {
                                // keep the split
                                didSomething = true;
                                newTiles.push_back(res.first);
                                newTiles.push_back(res.second);
                        } else {
                                // keep the original
                                newTiles.push_back(it);
                        }
                }

                matrixTiles = newTiles;

                if (!didSomething)
                        break;
        }

        auto oldPrecision = cout.precision();
        cout.precision(2);

        // compute the zero fraction
        size_t zeroElements = 0;
        for (size_t i = 0; i < P.nCols(); i++)
                zeroElements += P.nRows() - 4 * motifs[i].size();

        double zeroFrac = (double)zeroElements / (double)(P.nRows() * P.nCols());
        cout << "Matrix P initially contains " << 100.0*zeroFrac << "% zeros\n";
        cout << "Matrix P has been partitioned into " << matrixTiles.size() << " tile(s):\n";

        for (const auto& it : matrixTiles)
                cout << "\t" << it << "\n";

        for (const auto& it : matrixTiles)
                zeroElements -= (P.nRows() - it.rowEnd) * (it.colEnd - it.colStart);
        zeroFrac = (double)zeroElements / (double)(P.nRows() * P.nCols());
        cout << "Tiled matrix P contains " << 100.0*zeroFrac << "% zeros\n";

        cout.precision(oldPrecision);

        /*ofstream ofs("hist.dat");
        for (size_t i = 0; i < height.size(); i++)
                ofs << i+1 << "\t" << height[i] << "\n";*/
}

void MotifContainer::writeMotifNames(const std::string& filename)
{
        ofstream ofs(filename.c_str());

        for (auto m : motifs)
                ofs << m.getName() << "\n";

        ofs.close();
}

void MotifContainer::writePossumFile(const std::string& filename)
{
        ofstream ofs(filename.c_str());
        for (auto& m : motifs) {
                ofs << "BEGIN GROUP" << endl;
                ofs << "BEGIN FLOAT" << endl;
                ofs << "ID " << m.getName() << endl;
                ofs << "AC " << "dummy" << endl;
                ofs << "DE " << "dummy description" << endl;
                ofs << "AP DNA" << endl;
                ofs << "LE " << m.size() << endl;
                for (size_t i = 0; i < m.size(); i++)
                        ofs << "MA " << m[i][0] << " " << m[i][1] << " "
                            << m[i][2] << " " << m[i][3] << endl;
                ofs << "END" << endl;
                ofs << "END" << endl;
        }
}

void MotifContainer::writeMOODSFiles()
{
        for (const auto& it : motifs)
                it.writeMOODSFile(it.getName() + ".pfm");
}

size_t MotifContainer::getMaxMotifLen() const
{
        size_t maxMotifLen = 0;
        for (const auto& it : motifs)
                maxMotifLen = max<size_t>(maxMotifLen, it.size());
        return maxMotifLen;
}

// ============================================================================
// MOTIF OCCURRENCES
// ============================================================================

ostream& operator<< (ostream& os, const MotifOccurrence& m)
{
        os << m.getMotifID() << "\t" << m.getSpeciesID() << "\t"
           << m.getSequenceID() << "\t" << m.getSequencePos() << "\t"
           << m.getStrand() << "\t" << m.getScore();
        return os;
}
