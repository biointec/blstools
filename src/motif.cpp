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

        for (int i = 0; i < numBins; i++) {
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
        size_t bgTotCounts = accumulate(bgCounts.begin(), bgCounts.end(), 0);
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
                size_t totCounts = accumulate(PFM[i].begin(), PFM[i].end(), 0);
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

        float score = 0.0;
        for (size_t i = 0; i < pattern.size(); i++) {
                size_t j = char2idx(pattern[i]);
                if (j < 4)      // if ACTG character
                        score += PWM[i][j];
        }

        return score;
}

float Motif::getMaxScore() const
{
        float maxScore = 0.0;

        for (auto& pos : PWM) {
                float maxAC = max<float>(pos[0], pos[1]);
                float maxGT = max<float>(pos[2], pos[3]);
                maxScore += max<float>(maxAC, maxGT);
        }

        return maxScore;
}

float Motif::getMinScore() const
{
        float minScore = 0.0;

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
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + filename);

        vector<Motif> allMotifs; size_t motifID = 0;
        while (ifs.good()) {
                string temp;
                getline(ifs, temp);
                if (temp.empty())
                        continue;
                if (temp.front() == '>') {      // add a new motif
                        allMotifs.push_back(Motif(temp.substr(1), motifID++));
                        continue;
                }
                istringstream iss(temp);
                size_t A, C, G, T;
                iss >> A >> C >> G >> T;

                allMotifs.back().addCharacter({A, C, G, T});
        }

        // copy the temporary allmotifs structure to motifs
        for (const auto& motif : allMotifs) {
                if (loadPermutations || !motif.isPermutation())
                        motifs.push_back(motif);
        }

        ifs.close();
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

void MotifContainer::generateMatrix()
{
        // allocate memory for matrix P
        size_t m = motifs.size();
        size_t k = 4 * getMaxMotifLen();

        P.setDimensions(m, k);
        P.allocateMemory(0);

        // fill the matrix
        for (size_t i = 0; i < motifs.size(); i++) {
                // fill in the forward motif
                const Motif& fwd = motifs[i];
                for (size_t j = 0; j < fwd.size(); j++)
                        for (size_t o = 0; o < 4; o++)
                                P(i, 4*j+o) = fwd[j][o];
                row2MotifID.push_back(i);
        }
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
        os << m.getMotifID() << "\t" << m.getSequenceID() << "\t"
           << m.getSequencePos() << "\t"  << m.getStrand()
           << "\t" << m.getScore();
        return os;
}
