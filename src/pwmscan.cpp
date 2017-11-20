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
#include "species.h"

using namespace std;

void PWMScan::printUsage() const
{
        cout << "Usage: blstools scan [options] motifs.input sequences.input\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help message\n";
        cout << "  -s\t--simple\tenable simple (slower) scan algorithm\n";
        cout << "  -rc\t--revcompl\talso search the reverse strand for occurrences\n\n";

        cout << " [options arg]\n";
        cout << "  -at\t--absthreshold\tset the minimal absolute score for a motif occurrence\n";
        cout << "  -rt\t--relthreshold\tset the minimal relative score [0..1] for a motif occurrence (default = 0.95)\n";
        cout << "  -pt\t--pthreshold\tcompute the motif score threshold from p-value [0..1]\n";
        cout << "  -t\t--numthreads\tset the number of parallel threads [default =  #cores]\n\n";

        cout << " [file_options]\n";
        cout << "  -H\t--histdir\tdirectory where the histogram file(s) are stored [default = .]\n";
        cout << "  -o\t--output\tfilename for the motif occurrences [default = occurences.txt]\n\n";

        cout << " File \"motifs.input\" should contain the motifs in Jaspar format\n";
        cout << "  (see documentation for specification)\n\n";

        cout << " File \"sequences.input\" should contain a list of input fasta files in the following format:\n";
        cout << "   speciesID_1\tspecies1_sequences.fasta\n";
        cout << "   speciesID_2\tspecies2_sequences.fasta\n";
        cout << "   speciesID_3\tspecies2_sequences.fasta\n";
        cout << "   ...\n";
        cout << "  (refer to the documentation for more info)\n\n";

        cout << " Example:\n";
        cout << "  blstools scan -o occurences.txt motifs.input sequences.input\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

void PWMScan::extractOccurrences(const Matrix<float>& R, size_t offset,
                                 SeqMatrix& sm,
                                 const MotifContainer& motifContainer,
                                 vector<MotifOccurrence>& motifOcc)
{
        for (size_t j = 0; j < sm.getNumOccCol(); j++) {
                for (size_t i = 0; i < R.nRows(); i++) {
                        float thisScore = R(i,j);
                        size_t motifIdx = motifContainer.getMotifIDAtRow(i);
                        const Motif& m = motifContainer[motifIdx];

                        if (thisScore < (m.getThreshold() - 1e-5))
                                continue;

                        // at this point an occurrence is found
                        SeqPos seqPos = sm.getSeqPos(offset, j);
                        size_t remSeqLen = sm.getRemainingSeqLen(offset, j);

                        if (m.size() > remSeqLen)
                                continue;

                        char strand = m.isRevCompl() ? '-' : '+';

                        motifOcc.push_back(MotifOccurrence(m.getID(), seqPos.getSeqIndex(),
                                                           seqPos.getSeqPos(), strand, thisScore));
                }
        }
}

void PWMScan::scanThread(size_t speciesID, const MotifContainer& motifContainer,
                         FastaBatch& seqBatch, ostream& os)
{
        size_t overlap = motifContainer.getMaxMotifLen() - 1;
        size_t K = 250;                 // choose freely
        size_t W = 4*K;                 // choose freely

        // pattern matrix
        const Matrix<float> &P = motifContainer.getMatrix();

        // sequence matrix
        SeqMatrix sm(K, overlap, W);

        // result matrix
        Matrix<float> R(P.nRows(), W);

        vector<MotifOccurrence> occurrences;

        while (sm.getNextSeqMatrix(seqBatch)) {
                for (int offset = 0; offset < K; offset++) {
                        SubMatrix<float> subS = sm.getSubMatrix(offset);
                        R.gemm(P, subS);
                        extractOccurrences(R, offset, sm, motifContainer, occurrences);
                }

                // write the occurrences to disk
                lock_guard<mutex> lock(myMutex);
                for (auto o : occurrences)
                        os << o.getMotifID() << "\t" << speciesID << "\t"
                           << o.getSequenceID() << "\t" << o.getSequencePos()
                           << "\t" << o.getStrand() << "\t"  << o.getScore() << "\n";

                totMatches += occurrences.size();
                occurrences.clear();

                cout << "."; cout.flush();
        }
}

void PWMScan::scanThreadSimple(size_t speciesID, const MotifContainer& motifContainer,
                               FastaBatch& seqBatch, std::ostream& os)
{
        vector<MotifOccurrence> occurrences;

        string line; SeqPos seqPos;
        while (seqBatch.getNextFilteredLine(line, seqPos)) {
                for (size_t i = 0; i < line.size(); i++) {
                        for (size_t j = 0; j < motifContainer.size(); j++) {
                                const Motif& m = motifContainer[j];
                                if (m.size() > (line.size() - i))
                                        continue;
                                float thisScore = m.getScore(line.substr(i, m.size()));
                                if (thisScore < (m.getThreshold() - 1e-5))
                                        continue;

                                char strand = m.isRevCompl() ? '-' : '+';

                                occurrences.push_back(MotifOccurrence(m.getID(), seqPos.getSeqIndex(),
                                                                      seqPos.getSeqPos() + i, strand, thisScore));
                        }
                }

                // write the occurrences to disk
                lock_guard<mutex> lock(myMutex);
                for (auto o : occurrences)
                        os << o.getMotifID() << "\t" << speciesID << "\t"
                           << o.getSequenceID() << "\t" << o.getSequencePos()
                           << "\t" << o.getStrand() << "\t"  << o.getScore() << "\n";

                totMatches += occurrences.size();
                occurrences.clear();
        }
}


void PWMScan::scanPWM(size_t speciesID, const MotifContainer& motifContainer,
                      FastaBatch& seqBatch, std::ostream& os)
{
        cout << "Using " << numThreads << " threads" << endl;

        // start histogram threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&PWMScan::scanThread, this, speciesID,
                                          cref(motifContainer),
                                          ref(seqBatch), ref(os));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        cout << endl;
}

void PWMScan::scanPWMSimple(size_t speciesID, const MotifContainer& motifContainer,
                            FastaBatch& seqBatch, std::ostream& os)
{
        cout << "Using " << numThreads << " threads" << endl;

        // start histogram threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&PWMScan::scanThreadSimple, this, speciesID,
                                          cref(motifContainer),
                                          ref(seqBatch), ref(os));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        cout << endl;
}

PWMScan::PWMScan(int argc, char ** argv) : simpleMode(false),
        outputFilename("occurrences.txt"),
        absThSpecified(false), absThreshold(0.0), relThSpecified(false),
        relThreshold(0.95), pvalueSpecified(false), pvalue(0.001),
        numThreads(thread::hardware_concurrency()), revCompl(false),
        pseudocounts(1), totMatches(0)
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
                } else if ((arg == "-rc") || (arg == "--revcompl")) {
                        revCompl = true;
                } else if ((arg == "-s") || (arg == "--simple")) {
                        simpleMode = true;
                } else if (((arg == "-at") || (arg == "--absthreshold")) && (i+1 < argc-2)) {
                        absThSpecified = true;
                        absThreshold = atof(argv[i+1]);
                        i++;
                } else if (((arg == "-rt") || (arg == "--relthreshold")) && (i+1 < argc-2)) {
                        relThSpecified = true;
                        relThreshold = atof(argv[i+1]);
                        if ((relThreshold < 0.0) || (relThreshold > 1.0))
                                throw runtime_error("The relative threshold should be in range [0..1].");
                        i++;
                } else if (((arg == "-pt") || (arg == "--pthreshold")) && (i+1 < argc-2)) {
                        pvalueSpecified = true;
                        pvalue = atof(argv[i+1]);
                        if ((pvalue < 0.0) || (pvalue > 1.0))
                                throw runtime_error("The p-value should be in range [0..1].");
                        i++;
                } else if (((arg == "-t") || (arg == "--numthreads")) && (i+1 < argc-2)) {
                        numThreads = atoi(argv[i+1]);
                        if (numThreads < 1)
                                throw runtime_error("Number of threads must be a non-zero positive number");
                        i++;
                } else if (((arg == "-H") || (arg == "--histdir")) && (i+1 < argc-2)) {
                        histdir = string(argv[i+1]);
                        if (histdir.back() != '/')
                                histdir.push_back('/');
                        i++;
                } else if (((arg == "-o") || (arg == "--output")) && (i+1 < argc-2)) {
                        outputFilename = string(argv[i+1]);
                        i++;
                } else {
                        printUsage();
                        exit(EXIT_FAILURE);
                }
        }

        cout << "Welcome to fast PWM scan" << endl;

        // If none of the thresholds is specified, default to relative threshold
        if (!(absThSpecified || relThSpecified || pvalueSpecified))
                relThSpecified = true;

        // Don't specify multiple thresholds
        if (absThSpecified && relThSpecified)
                throw runtime_error("Specify either the absolute or relative threshold, not both.");
        if (absThSpecified && pvalueSpecified)
                throw runtime_error("Specify either the absolute or p-value threshold, not both.");
        if (relThSpecified && pvalueSpecified)
                throw runtime_error("Specify either the relative or p-value threshold, not both.");

        // A) Load the manifest file
        string manifestFilename(argv[argc-1]);
        string dictFilename = manifestFilename + ".dict";

        SpeciesContainer speciesContainer;
        speciesContainer.load(dictFilename);

        // B) Load the motifs
        string motifFilename = string(argv[argc-2]);

        MotifContainer motifContainer(motifFilename, true);
        cout << "Loaded " << motifContainer.size() << " motifs from disk\n";
        cout << "Maximum motif size: " << motifContainer.getMaxMotifLen() << endl;

        if (revCompl) {
                motifContainer.addReverseComplements();
                cout << "Added reverse complementary motifs" << endl;
        }

       // motifContainer.writePossumFile("motifs.possum");
        motifContainer.writeMOODSFiles();

        // C) Scan the sequences
        ofstream ofsCutoff("motifCutoff.txt");

        ofstream ofs;
        if (!outputFilename.empty())
                ofs.open(outputFilename.c_str());
        ostream& os = (outputFilename.empty()) ? cout : ofs;

        size_t speciesID = 0;
        for (auto species : speciesContainer) {
                cout << "Scanning species: " << species.getName() << endl;

                // generate the PWM using the background counts for those species
                for (auto& motif : motifContainer)
                        motif.PFM2PWM(species.getNuclCounts(), pseudocounts);

                // set or compute the thresholds
                if (absThSpecified) {
                        cout << "Absolute motif score threshold set to: " << absThreshold << endl;
                        for (auto& motif : motifContainer) {
                                motif.setThreshold(absThreshold);
                                ofsCutoff << species.getName() << "\t" << motif.getName() << "\t" << motif.getThreshold() << endl;
                        }
                } else if (relThSpecified) {
                        cout << "Relative motif score threshold set to: " << relThreshold << endl;
                        for (auto& motif : motifContainer) {
                                float maxScore = motif.getMaxScore();
                                float minScore = motif.getMinScore();
                                float threshold = relThreshold * (maxScore - minScore) + minScore;
                                motif.setThreshold(threshold);
                                ofsCutoff << species.getName() << "\t" << motif.getName() << "\t" << motif.getThreshold() << endl;
                        }
                } else if (pvalueSpecified) {
                        cout << "P-value motif score threshold set to: " << pvalue << endl;

                        for (auto& motif : motifContainer) {
                                ScoreHistogram hist(motif.getMinScore(), motif.getMaxScore(), 250);
                                hist.loadHistogram(histdir, "hist_" + species.getName() + "_" +  motif.getBaseName());
                                motif.setThreshold(hist.getScoreCutoff(pvalue));
                                ofsCutoff << species.getName() << "\t" << motif.getName() << "\t" << motif.getMinScore() << "\t" << motif.getThreshold() << "\t" << motif.getMaxScore() << endl;
                        }
                }

                // from these PWMs, generate the pattern matrix
                motifContainer.generateMatrix();

                vector<string> filenames = species.getSequenceFilenames();
                FastaBatch seqBatch(filenames);

                // scan the sequences for PWM occurrences
                if (simpleMode) {
                        scanPWMSimple(speciesID++, motifContainer, seqBatch, os);
                } else {
                        scanPWM(speciesID++, motifContainer, seqBatch, os);
                }
        }

        if (!outputFilename.empty())
                ofs.close();

        ofsCutoff.close();

        cout << "\nWrote " << totMatches << " matches to " << outputFilename << ".\n";
}
