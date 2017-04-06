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

#include "motif.h"

using namespace std;

// ============================================================================
// MOTIF
// ============================================================================

int char2idx(char c)
{
        if (c == 'A' || c == 'a')
                return 0;
        if (c == 'C' || c == 'c')
                return 1;
        if (c == 'G' || c == 'g')
                return 2;
        return 3;
}

void Motif::posFreq2PWM(const array<float, 4>& bgFreq)
{
        for (auto& pos : motif) {
                float totFreq = 0;
                for (size_t i = 0; i < 4; i++)
                        totFreq += pos[i];

                for (size_t i = 0; i < 4; i++)
                        pos[i] = log2((pos[i] + 1) / (totFreq + 4) / bgFreq[i]);
        }
}

float Motif::getScore(const std::string& pattern) const
{
        // make sure the pattern and motif have the same size
        assert (pattern.size() == size());

        float score = 0.0;
        for (size_t i = 0; i < pattern.size(); i++) {
                //cout << pattern[i] << " " << motif[i][char2idx(pattern[i])] << endl;
                score += motif[i][char2idx(pattern[i])];
        }

        return score;
}

float Motif::getMaxScore() const
{
        float maxScore = 0.0;

        for (auto& pos : motif) {
                float maxAC = max<float>(pos[0], pos[1]);
                float maxGT = max<float>(pos[2], pos[3]);
                maxScore += max<float>(maxAC, maxGT);
        }

        return maxScore;
}

float Motif::getMinScore() const
{
        float minScore = 0.0;

        for (auto& pos : motif) {
                float minAC = min<float>(pos[0], pos[1]);
                float minGT = min<float>(pos[2], pos[3]);
                minScore += min<float>(minAC, minGT);
        }

        return minScore;
}


void Motif::revCompl()
{
        name.append("_RC");
        reverse(motif.begin(), motif.end());

        for (size_t i = 0; i < motif.size(); i++) {
                array<float, 4> copy = motif[i];
                motif[i] = array<float, 4>{copy[3], copy[2], copy[1], copy[0]};
        }

}

ostream& operator<< (ostream& os, const Motif& m)
{
        os << m.name << "\n";
        for (auto pos : m.motif)
                os << pos[0] << " " << pos[1] << " " << pos[2] << " " << pos[3] << "\n";
        return os;
}

// ============================================================================
// MOTIF CONTAINER
// ============================================================================

MotifContainer::MotifContainer(const std::string& filename,
                               const std::array<float, 5>& bgFreq)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + filename);

        while (ifs.good()) {
                string temp;
                getline(ifs, temp);
                if (temp.empty())
                        continue;
                if (temp.front() == '>') {
                        motifs.push_back(Motif(temp.substr(1)));
                        continue;
                }
                istringstream iss(temp);
                int a, b, c, d;
                iss >> a >> b >> c >> d;

                motifs.back().addCharacter({(float)a, (float)b, (float)c, (float)d});
        }

        ifs.close();

        for (auto& m : motifs)
                m.posFreq2PWM({0.25, 0.25, 0.25, 0.25}/*{bgFreq[0], bgFreq[1], bgFreq[2], bgFreq[3]}*/);

        /*ofstream ofs("jaspar2.possum");
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
        }*/
}

void MotifContainer::generateMatrix(Matrix<float>& M)
{
        M.fill(0);
        size_t numMotifs = motifs.size();

        for (size_t i = 0; i < numMotifs; i++)
                for (size_t j = 0; j < motifs[i].size(); j++)
                        for (size_t o = 0; o < 4; o++)
                                M(i, 4*j+o) = motifs[i][j][o];
}

void MotifContainer::addReverseCompl()
{
        vector<Motif> copy = motifs;
        motifs.insert(motifs.end(), copy.begin(), copy.end());

        for (size_t i = motifs.size() / 2; i < motifs.size(); i++)
                motifs[i].revCompl();
}

void MotifContainer::setRelThreshold(float relThreshold)
{
        threshold.clear();
        threshold.reserve(motifs.size());

        for (const auto& m : motifs) {
                float maxScore = m.getMaxScore();
                float minScore = m.getMinScore();

                threshold.push_back(relThreshold * (maxScore - minScore) + minScore);

                //cout << "Motif " << m.getName() << "[min: " << minScore
                //     << ", max: " << maxScore << ", th: " << threshold.back() << "]" << endl;
        }
}

void MotifContainer::writeMotifNames(const std::string& filename)
{
        ofstream ofs(filename.c_str());

        for (auto m : motifs)
                ofs << m.getName() << "\n";

        ofs.close();
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

MotifOccurrence::MotifOccurrence(size_t motifID_, size_t sequenceID_,
                                 size_t sequencePos_, float score_) :
                                 motifID(motifID_), sequenceID(sequenceID_),
                                 sequencePos(sequencePos_), score(score_)
{
}

ostream& operator<< (ostream& os, const MotifOccurrence& m)
{
        //os << "Seq ID: " << m.sequenceID << ", seq pos: " << m.sequencePos
        //   << ", motif ID: " << m.getMotifID() << ", score:" << m.score;
        os << m.getMotifID() << "\t" << m.getSequenceID() << "\t"
           << m.getSequencePos() << "\t" << "\t" << m.getScore();
        return os;
}
