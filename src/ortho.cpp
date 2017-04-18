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

#include "ortho.h"
#include "phylotree.h"

using namespace std;

void Ortho::printUsage() const
{
        cout << "Usage: blstools ortho [options] occurrences.input sequences.idx motifs.idx orthogroups.input phylotree.newick\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help message\n\n";
        cout << "  -p\t--print\t\toutput tree to the screen\n\n";

        cout << " [file_options]\n";
        cout << "  -o\t--output\tfilename for output BLS values [default = stdout]\n\n";

        cout << " File \"tree.newick\" should contain the phylogenetic tree in Newick format, e.g.:\n";
        cout << "  \"(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);\"\n\n";

        cout << " File \"leafs.input\" should contain a list of space or tab separated\n";
        cout << "  leaf names for which the BLS should be computed (one line per entry)\n\n";

        cout << " Example:\n";
        cout << "  blstools bls -o output.bls tree.newick leafs.input\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

string gene2species(const string& gene)
{
        char c = gene[0];
        if (c == 'B')
                return "bdi";
        if (c == 'H')
                return "bdi";
        if (c == 'M')
                return "mac";
        if (c == 'Z')
                return "zma";
        if (c == 'O')
                return gene[2] == 'I' ? "osindica" : "osa";
        if (c == 'S')
                return gene[1] == 'B' ? "sbi" : "sit";
        cerr << "Error converting gene to species" << endl;
        return "ERROR";
}

Ortho::Ortho(int argc, char ** argv)
{
        // check for sufficient arguments
        if (argc < 7) {
                printUsage();
                exit(EXIT_FAILURE);
        }

        // process optional arguments
       /* for (int i = 2; i < argc-2; i++) {
                string arg(argv[i]);

                if ((arg == "-h") || (arg == "--help")) {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else if ((arg == "-p") || (arg == "--print")) {
                        printTree = true;
                } else if (((arg == "-o") || (arg == "--output")) && (i+1 < argc-2)) {
                        outputFilename = string(argv[i+1]);
                        i++;
                } else {
                        printUsage();
                        exit(EXIT_FAILURE);
                }
        }*/

        string occFilename = argv[argc-5];
        string seqIdxFilename = argv[argc-4];
        string motifIdxFilename = argv[argc-3];
        string orthoFilename = argv[argc-2];
        string phyloFilename = argv[argc-1];

        // A) Read the ortho groups
        ifstream ifs(orthoFilename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + orthoFilename);

        string temp, orthoname, seqName, motifName;
        while (ifs) {
                ifs >> temp >> orthoname >> temp >> seqName;
                if (!ifs)
                        break;
                orthoGroups[orthoname].insert(seqName);
                seq2ortho.insert(pair<string, string>(seqName, orthoname));
        }

        // also insert the zma gene into the ortho group
        for (auto it : orthoGroups) {
                it.second.insert(it.first);
                seq2ortho.insert(pair<string, string>(orthoname, orthoname));
        }

        cout << "Read " << orthoGroups.size() << " ortho groups." << endl;

        ifs.close();

        // B) Read the index file
        ifs.open(seqIdxFilename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + seqIdxFilename);

        while (ifs) {
                ifs >> seqName >> temp >> temp >> temp;
                if (!ifs)
                        break;
                seqIndex.push_back(seqName);
        }

        ifs.close();

        // C) Read the motif index file
        ifs.open(motifIdxFilename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + motifIdxFilename);

        while (ifs) {
                ifs >> motifName;
                if (!ifs)
                        break;
                motifIndex.push_back(motifName);
        }

        ifs.close();

        // process the input tree
        ifs.open(phyloFilename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + string(argv[argc-2]));
        string newickStr;
        getline(ifs, newickStr);
        ifs.close();

        PhylogeneticTree pt(newickStr);
        pt.normalizeBranchLength();

        // C) Read the occurrence file
        ofstream ofs("BLS.txt");

        ifs.open(occFilename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + occFilename);

        int motifID, prevMotifID = -1, seqID;
        map<string, set<string> > orthoSpecComb;
        map<string, set<string> > orthoGeneComb;

        while (ifs) {
                ifs >> motifID >> seqID >> temp >> temp;
                if (!ifs)
                        break;

                if (prevMotifID == -1)
                        prevMotifID = motifID;

                if (motifID != prevMotifID) {
                        array<size_t, 10> counts = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
                        for (auto it : orthoSpecComb) {
                                float BLS = pt.getBLS(it.second);

                                ofs << prevMotifID << "\t" << motifIndex[prevMotifID] << "\t" << it.first << "\t" << BLS;

                                for (auto gene : orthoGeneComb[it.first])
                                        ofs << "\t" << gene;
                                ofs << endl;

                                int end = 0;
                                if (BLS >= 0.00)
                                        end = 1;
                                if (BLS >= 0.10)
                                        end = 2;
                                if (BLS >= 0.20)
                                        end = 3;
                                if (BLS >= 0.30)
                                        end = 4;
                                if (BLS >= 0.40)
                                        end = 5;
                                if (BLS >= 0.50)
                                        end = 6;
                                if (BLS >= 0.60)
                                        end = 7;
                                if (BLS >= 0.70)
                                        end = 8;
                                if (BLS >= 0.80)
                                        end = 9;
                                if (BLS >= 0.90)
                                        end = 10;

                                for (size_t i = 0; i < end; i++)
                                        counts[i]++;
                        }

                        cout << prevMotifID << "\t" << motifIndex[prevMotifID];
                        for (int i = 0; i < 10; i++)
                                cout << "\t" << counts[i];
                        cout << "\n";

                        orthoSpecComb.clear();
                        orthoGeneComb.clear();
                        prevMotifID = motifID;
                }

                string seq = seqIndex[seqID];
                auto range = seq2ortho.equal_range(seq);

                for (auto it = range.first; it != range.second; it++) {
                        orthoSpecComb[it->second].insert(gene2species(seq));
                        orthoGeneComb[it->second].insert(seq);
                }
        }


        ifs.close();
        ofs.close();
}
