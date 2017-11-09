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

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "ortho.h"
#include "phylotree.h"
#include "species.h"
#include "motif.h"

using namespace std;

void OrthoCount::print()
{
        cout << motifName << endl;
        for (auto it : motifCounts)
                cout << it << " ";
        cout << "\n";
        vector<float> median = getRandomMedian();
        for (auto it : median)
                cout << it << " ";
        cout << "\n";
        vector<float> C = computeCScores();
        for (auto it : C)
                cout << it << " ";
        cout << "\n";
}

// ============================================================================
// ORTHOCONTAINER
// ============================================================================

void OrthoContainer::load(const string& filename)
{
        // read the input file and create the orthoGroups
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + filename);

        while (true) {
                string refSpecies, refSeqName, orthoSpecies, orthoSeqName;
                ifs >> refSpecies >> refSeqName >> orthoSpecies >> orthoSeqName;
                if (!ifs)
                        break;

                orthoGroups[refSeqName].insert(refSpecies, refSeqName);
                orthoGroups[refSeqName].insert(orthoSpecies, orthoSeqName);
        }

        // create a sequence 2 ortho index
        for (auto it = orthoGroups.begin(); it != orthoGroups.end(); it++) {
                const string& orthoName = it->first;
                const OrthoGroup& og = it->second;
                for (auto it2 = og.seqBegin(); it2 != og.seqEnd(); it2++)
                        seq2ortho.insert(make_pair(*it2, orthoName));
        }
}

// ============================================================================
// ORTHO MODULE
// ============================================================================

std::ostream& operator<< (std::ostream& os, const OrthoGroup& og)
{
        for (auto it : og.sequences)
                os << it << "\n";
        return os;
}

void Ortho::printUsage() const
{
        cout << "Usage: blstools ortho [options] motifs.input sequences.input orthogroups.input phylotree.input occurrences.input\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help message\n";

        cout << " [file_options]\n";
        cout << "  -o\t--output\tfilename for output BLS values [default = stdout]\n\n";

        cout << " File \"phylotree.input\" should contain the phylogenetic tree in Newick format, e.g.:\n";
        cout << "  \"(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);\"\n\n";

        cout << " File \"leafs.input\" should contain a list of space or tab separated\n";
        cout << "  leaf names for which the BLS should be computed (one line per entry)\n\n";

        cout << " Example:\n";
        cout << "  blstools bls -o output.bls tree.newick leafs.input\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

Ortho::Ortho(int argc, char ** argv)
{
        // check for sufficient arguments
        if (argc < 7) {
                printUsage();
                exit(EXIT_FAILURE);
        }

        // process optional arguments
        for (int i = 2; i < argc-5; i++) {
                string arg(argv[i]);

                if ((arg == "-h") || (arg == "--help")) {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else {
                        printUsage();
                        exit(EXIT_FAILURE);
                }
        }

        string motifFilename = argv[argc-5];
        string manifestFilename = argv[argc-4];
        string orthoFilename = argv[argc-3];
        string phyloFilename = argv[argc-2];
        string occFilename = argv[argc-1];

        // A) Load the motifs
        MotifContainer motifContainer(motifFilename, true);
        cout << "Loaded " << motifContainer.size() << " motifs from disk\n";

        // B) Load the manifest file
        string dictFilename = manifestFilename + ".dict";

        SpeciesContainer speciesContainer;
        speciesContainer.load(dictFilename);
        cout << "Loaded dictionaire with " << speciesContainer.size() << " species\n";

        // C) Read the ortho groups
        OrthoContainer orthoContainer;
        orthoContainer.load(orthoFilename);
        cout << "Loaded " << orthoContainer.size() << " orthology groups\n";

         // D) process the phylogenetic tree
        ifstream ifs(phyloFilename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + phyloFilename);
        string newickStr;
        getline(ifs, newickStr);
        ifs.close();

        PhylogeneticTree pt(newickStr);
        pt.normalizeBranchLength();

        // check whether the species in the phylotree match with the ones in the manifest file
        set<string> speciesNames = pt.getAllNames();
        for (auto it : speciesContainer)
                if (speciesNames.find(it.getName()) == speciesNames.end())
                        throw runtime_error("ERROR: Species name \"" + it.getName() + "\" does not occur in " + phyloFilename);

        cout << "Loaded phylogenetic tree" << endl;

        // E) Read the occurrence file
        ofstream ofs("BLS.txt");
        ofstream ofsCounts("motifCounts.txt");

        ifs.open(occFilename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + occFilename);

        size_t motifID, speciesID, nextMotifID, seqID;
        map<string, set<string> > orthoSpecComb;
        map<string, set<string> > orthoGeneComb;

        ifs >> motifID;
        size_t numBLSIntv = 10;

        while (ifs) {
                string temp;
                // motifID speciesID seqID seqPos strand score
                ifs >> speciesID >> seqID >> temp >> temp >> temp;
                ifs >> nextMotifID;

                if (speciesID > speciesContainer.size())
                        throw runtime_error("ERROR: File " + occFilename + " contains a speciesID not present in file " + manifestFilename);

                if (motifID > motifContainer.size())
                        throw runtime_error("ERROR: File " + occFilename + " contains a motifID not present in file " + motifFilename);

                if (nextMotifID < motifID)
                        throw runtime_error("ERROR: File " + occFilename + " is not sorted. Sort this file prior to running the ortho module");

                string species = speciesContainer.getSpecies(speciesID).getName();
                string seq = speciesContainer.getSpecies(speciesID).getSeqName(seqID);

                // for all ortho groups that contain sequence "seq"
                auto range = orthoContainer.equal_range(seq);
                for (auto it = range.first; it != range.second; it++) {
                        orthoSpecComb[it->second].insert(species);
                        orthoGeneComb[it->second].insert(seq);
                }

                // we have more hits from the same motif: continue
                if ((ifs) && (nextMotifID == motifID))
                        continue;

                vector<size_t> counts(numBLSIntv, 0);
                for (auto it : orthoSpecComb) {
                        const OrthoGroup& og = orthoContainer.getOrthoGroup(it.first);

                        float BLS = pt.getBLS(it.second);
                        // float maxBLS = pt.getBLS(og.getSpecies());
                        //if (maxBLS > 0)
                         //       BLS /= maxBLS;

                        ofs << motifID << "\t" << motifContainer[motifID].getName() << "\t" << it.first << "\t" << BLS;

                        for (auto gene : orthoGeneComb[it.first])
                                ofs << "\t" << gene;
                        ofs << endl;

                        float intWidth = 1.0 / numBLSIntv;
                        int end = floor(BLS / intWidth) + 1;
                        if (end > numBLSIntv)
                                end = numBLSIntv;

                        for (size_t i = 0; i < end; i++)
                                counts[i]++;
                }

                cout << motifID << "\t" << motifContainer[motifID].getName();
                ofsCounts << motifID << "\t" << motifContainer[motifID].getName();
                for (int i = 0; i < numBLSIntv; i++) {
                        cout << "\t" << counts[i];
                        ofsCounts << "\t" << counts[i];
                }
                cout << "\n";
                ofsCounts << "\n";

                orthoSpecComb.clear();
                orthoGeneComb.clear();

                motifID = nextMotifID;
        }

        ifs.close();
        ofs.close();
        ofsCounts.close();

        ofs.open("motifCScores.txt");
        ifs.open("motifCounts.txt");
        if (!ifs)
                throw runtime_error("Could not open file: counts.txt");

        OrthoCount orthoCount("");
        while (ifs) {
                size_t motifID;
                vector<size_t> counts(numBLSIntv, 0);
                string motifName;

                ifs >> motifID >> motifName;
                for (size_t i = 0; i < numBLSIntv; i++)
                        ifs >> counts[i];

                // if we encounter a new motif
                if (!motifContainer[motifID].isPermutation()) {
                        if (orthoCount.size() > 0) {
                                ofs << motifName;
                                vector<float> C = orthoCount.computeCScores();
                                for (auto it : C)
                                        ofs << "\t" << it;
                                ofs << "\n";
                        }

                        orthoCount = OrthoCount(motifContainer[motifID].getName());
                        orthoCount.setMotifCounts(counts);
                } else {        // it's a random motif
                        orthoCount.addRandomCounts(counts);
                }
        }

        if (orthoCount.size() > 0) {
                ofs << orthoCount.getName();
                vector<float> C = orthoCount.computeCScores();
                for (auto it : C)
                        ofs << "\t" << it;
                ofs << "\n";
        }

        ofs.close();
}
