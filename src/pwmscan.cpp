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
#include <iomanip>
#include <thread>
#include <algorithm>

#include "pwmscan.h"
#include "sequence.h"

using namespace std;

void PWMScan::printUsage() const
{
        cout << "Usage: blstools scan [options] sequences.input motifs.input\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help message\n";
        cout << "  -rc\t--revcompl\talso search the reverse strand for occurrences [default = no]\n";
        cout << "  -at\t--absthreshold\tset the minimal absolute score for a motif occurrence\n";
        cout << "  -rt\t--relthreshold\tset the minimal relative score [0..1] for a motif occurrence\n";
        cout << "  -t\t--numthreads\tset the number of parallel threads [default = 1]\n";
        cout << "  -l\t--list\t\tprovide a list of input fasta filenames\n\n";

        cout << " [file_options]\n";
        cout << "  -o\t--output\tfilename for the motif occurrences [default = stdout]\n\n";

        cout << " File \"motifs.input\" should contain the motifs in Jaspar format\n";
        cout << "  (see documentation for specification)\n\n";

        cout << " Example:\n";
        cout << "  blstools scan -o output.txt sequences.fasta motifs.input\n";
        cout << "  blstools scan -o output.txt -l sequencelist.txt motifs.input\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

void PWMScan::loadFasta(const string& filename, string& sequence)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + filename);

        string temp;
        ifs >> temp;

        while (ifs.good()) {
                ifs >> temp;
                sequence.append(temp);
        }

        ifs.close();
}

void PWMScan::countFrequencies(const string& sequence,
                               array<float, 5>& frequencies)
{
        frequencies.fill(0);

        for (size_t i = 0; i < sequence.size(); i++) {
                if (sequence[i] == 'A')
                        frequencies[0]++;
                else if (sequence[i] == 'C')
                        frequencies[1]++;
                else if (sequence[i] == 'G')
                        frequencies[2]++;
                else if (sequence[i] == 'T')
                        frequencies[3]++;
                else
                        frequencies[4]++;
        }

        for (int i = 0; i < 5; i++)
                frequencies[i] /= sequence.size();
}

void PWMScan::scanThread(size_t myID, const MotifContainer& motifs,
                         const Matrix<float>& P, FastaBatch& seqBatch, ostream& os)
{
        cout << "Thread with threadID " << myID << " launched" << endl;

        SeqMatrix sm(K, overlap, n);
        while (sm.getNextSeqMatrix(seqBatch)) {

                vector<MotifOccurrence> motifOcc;
                sm.findOccurrences(P, motifOcc, motifs);

                lock_guard<mutex> lock(myMutex);
                for (size_t i = 0; i < motifOcc.size(); i++) {
                        size_t motifIdx = motifOcc[i].getMotifID();
                        os << motifOcc[i] << "\n";
                }

                totMatches += motifOcc.size();

                cout << ".";
                cout.flush();
        }

        cout << "Thread with threadID " << myID << " terminates" << endl;
}

PWMScan::PWMScan(int argc, char ** argv) : listProvided(false),
        absThSpecified(false), absThreshold(0.0), relThSpecified(false),
        relThreshold(1.0), numThreads(1), revCompl(false)
{
        // check for sufficient arguments
        if (argc < 4) {
                printUsage();
                exit(EXIT_FAILURE);
        }

        // process optional arguments
        for (int i = 2; i < argc-2; i++) {
                string arg(argv[i]);

                if ((arg == "-h") || (arg == "--help")) {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else if (((arg == "-at") || (arg == "--absthreshold")) && (i+1 < argc-2)) {
                        absThSpecified = true;
                        absThreshold = atof(argv[i+1]);
                        i++;
                } else if (((arg == "-rt") || (arg == "--relthreshold")) && (i+1 < argc-2)) {
                        relThSpecified = true;
                        relThreshold = atof(argv[i+1]);
                        i++;
                } else if (((arg == "-t") || (arg == "--numthreads")) && (i+1 < argc-2)) {
                        numThreads = atoi(argv[i+1]);
                        i++;
                }else if ((arg == "-rc") || (arg == "--revcompl")) {
                        revCompl = true;
                } else if ((arg == "-l") || (arg == "--list")) {
                        listProvided = true;
                } else if (((arg == "-o") || (arg == "--output")) && (i+1 < argc-2)) {
                        outputFilename = string(argv[i+1]);
                        i++;
                } else {
                        printUsage();
                        exit(EXIT_FAILURE);
                }
        }

        cout << "Welcome to fast PWM scan" << endl;

        // We cannot specify both the absolute and relative threshold
        if (absThSpecified && relThSpecified)
                throw runtime_error("Specify either the absolute or relative threshold, not both.");
        if (relThSpecified) {
                if ((relThreshold < 0.0) || (relThreshold > 1.0))
                        throw runtime_error("The relative threshold should be in range [0..1].");
        }

        if (numThreads < 1)
                throw runtime_error("Number of threads must be a non-zero positive number");

        // A) LOAD THE SEQUENCES AND COMPUTE THE BACKGROUND FREQUENCIES
        FastaBatch seqBatch;
        if (listProvided)
                seqBatch.addList(string(argv[argc-2]));
        else
                seqBatch.addFile(string(argv[argc-2]));

        cout << "Computing nucleotide frequencies for "
             << seqBatch.size() << " input files..." << endl;
        seqBatch.calcBGFrequencies();
        seqBatch.writeSeqNames(string(argv[argc-2]) + ".idx");
        cout << "Read " << seqBatch.getNumSequences()
             << " sequences with a total length of "
             << seqBatch.getTotalSeqLength() << endl;

        array<float, 5> bgFreq = seqBatch.getNucleotideFreq();

        cout << "Background frequencies: [" << setprecision(3)
             << "A:" << 100*bgFreq[0] << "%, "
             << "C:" << 100*bgFreq[1] << "%, "
             << "G:" << 100*bgFreq[2] << "%, "
             << "T:" << 100*bgFreq[3] << "%, "
             << "non-ACGT: " << 100*bgFreq[4] << "%]" << endl;

        // B) LOAD THE MOTIFS AND CONVERT TO POSITION WEIGHT MATRICES
        MotifContainer motifs(argv[argc-1], bgFreq);
        cout << "Loaded " << motifs.size() << " motifs from disk..." << endl;
        cout << "Maximum motif size: " << motifs.getMaxMotifLen() << endl;

        if (revCompl) {
                cout << "Adding reverse complementary motifs" << endl;
                motifs.addReverseCompl();
        }

        if (absThSpecified) {
                cout << "Absolute motif score threshold set to: " << absThreshold << endl;
                motifs.setAbsThreshold(absThreshold);
        } else {
                cout << "Relative motif score threshold set to: " << relThreshold << endl;
                motifs.setRelThreshold(relThreshold);
        }

        motifs.writeMotifNames(string(argv[argc-1]) + ".idx");

        // result matrix
        m = motifs.size();       // number of motifs
        n = 1000;                // choose freely
        Matrix<float> R(m, n);

        // pattern matrix
        k = 4 * motifs.getMaxMotifLen();
        Matrix<float> P(m, k);
        motifs.generateMatrix(P);

        // sequence matrix
        overlap = motifs.getMaxMotifLen() - 1;
        K = 250;                // choose freely

        ofstream ofs;
        if (!outputFilename.empty())
                ofs.open(outputFilename.c_str());
        ostream& os = (outputFilename.empty()) ? cout : ofs;

        seqBatch.reset();

        cout << "Scanning sequences (every dot corresponds to " << K*n
             << " nucleotides of sequence processed)" << endl;

        totMatches = 0;

        cout << "Number of threads: " << numThreads << endl;

        // start worker threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&PWMScan::scanThread, this,
                                          i, cref(motifs), cref(P), ref(seqBatch), ref(os));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        if (!outputFilename.empty())
                ofs.close();

        cout << "\nFound " << totMatches << " matches\n";
        cout << "Bye!" << endl;
}
