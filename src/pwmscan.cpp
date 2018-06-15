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
#include <sstream>
#include <future>

#include "pwmscan.h"
#include "sequence.h"
#include "species.h"

#ifdef HAVE_CUDA
        #include <cuda_runtime.h>
        #include <cublas_v2.h>
        #include "helper_cuda.h"
#endif

using namespace std;

extern void kernel_wrapper(float *R, int m, int n, float *threshold, int* occIdx, float* occScore, int *nOcc);

void PWMScan::printUsage() const
{
        cout << "Usage: blstools scan [options] motifs.input sequences.input\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help message\n";
        cout << "  -s\t--simple\tenable simple (slower) scan algorithm\n";
        cout << "  -c\t--cuda\t\tenable the CUDA scan algorithm\n";
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
                for (size_t i = 0; i < (size_t)R.nRows(); i++) {
                        float thisScore = R(i,j);
                        size_t motifIdx = motifContainer.getMotifIDAtRow(i);
                        const Motif& m = motifContainer[motifIdx];

                        if (thisScore < m.getThreshold())
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

void extractOccurrences2(int m, const map<int, int>& offset_v, int *occIdx, float *occScore,
			 SeqMatrix& sm, const MotifContainer& motifContainer,
			 vector<MotifOccurrence>& motifOcc)
{
	int idx = 0;
	for (const auto& it: offset_v) {
		int offset = it.first;
		int thisOcc = it.second;

		for (int c = 0; c < thisOcc; c++, idx++) {
			float thisScore = occScore[idx];
			int i = occIdx[idx] % m;
			int j = occIdx[idx] / m;

			size_t motifIdx = motifContainer.getMotifIDAtRow(i);
	                const Motif& m = motifContainer[motifIdx];

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

void PWMScan::scanThreadBLAS(size_t speciesID, const MotifContainer& motifContainer,
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
        const auto matBlock = motifContainer.getMatrixBlock();

        vector<MotifOccurrence> occurrences;

        while (sm.getNextSeqMatrix(seqBatch)) {
                for (size_t offset = 0; offset < K; offset++) {
                        SubMatrix<float> subS = sm.getSubMatrix(offset);
                        //R.gemm(P, subS);
                        R.gemm(P, subS, matBlock);
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

void PWMScan::scanThreadNaive(size_t speciesID, const MotifContainer& motifContainer,
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
                                if (thisScore < m.getThreshold())
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


void PWMScan::scanPWMBLAS(size_t speciesID, const MotifContainer& motifContainer,
                      FastaBatch& seqBatch, std::ostream& os)
{
        cout << "Using " << numThreads << " threads" << endl;

        // start histogram threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&PWMScan::scanThreadBLAS, this, speciesID,
                                          cref(motifContainer),
                                          ref(seqBatch), ref(os));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        cout << endl;
}

void PWMScan::scanPWMNaive(size_t speciesID, const MotifContainer& motifContainer,
                            FastaBatch& seqBatch, std::ostream& os)
{
        cout << "Using " << numThreads << " threads" << endl;

        // start histogram threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&PWMScan::scanThreadNaive, this, speciesID,
                                          cref(motifContainer),
                                          ref(seqBatch), ref(os));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        cout << endl;
}

void writeToDisk(size_t speciesID, const vector<MotifOccurrence> occurrences, mutex& m, std::ostream& os)
{
	ostringstream oss;
        for (auto o : occurrences)
        	oss << o.getMotifID() << "\t" << speciesID << "\t"
                    << o.getSequenceID() << "\t" << o.getSequencePos()
                    << "\t" << o.getStrand() << "\t"  << o.getScore() << "\n";

        // write the occurrences to disk
        lock_guard<mutex> lock(m);
        os << oss.str();
}

#ifdef HAVE_CUDA
void PWMScan::scanThreadCUBLAS(int devID, size_t speciesID,
                               const MotifContainer& motifContainer,
                               FastaBatch& seqBatch, std::ostream& os)
{
        const float alpha = 1.0f;
        const float beta = 0.0f;

        float *d_P = 0, *d_S = 0, *d_R = 0;
        float *d_threshold = 0, *d_occScore = 0;
	int *d_occIdx = 0, *d_nOcc = 0;

        size_t overlap = motifContainer.getMaxMotifLen() - 1;
        size_t K = 250;                 // choose freely
        size_t W = 4*K;                 // choose freely

        // set CUDA device for this thread
        cudaSetDevice(devID);

        // create a CUBLAS handle
        cublasHandle_t handle;
        cublasCreate(&handle);

        // pattern matrix
        const Matrix<float> &P = motifContainer.getMatrix();
        if (cudaMalloc((void **)&d_P, P.nRows() * P.nCols() * sizeof(float)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device\n");
        cublasSetVector(P.nRows() * P.nCols(), sizeof(float), P.data, 1, d_P, 1);

        // sequence matrix
        SeqMatrix sm(K, overlap, W);
        if (cudaMalloc((void **)&d_S, 4 * (K + overlap) * W * sizeof(float)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device\n");

        // result matrix
        Matrix<float> R(P.nRows(), W);
        if (cudaMalloc((void **)&d_R, R.nRows() * R.nCols() * sizeof(float)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device\n");

        // set the thresholds
        float *threshold = new float[R.nRows()];
        for (size_t i = 0; i < R.nRows(); i++) {
                size_t motifID = motifContainer.getMotifIDAtRow(i);
                const Motif& m = motifContainer[motifID];
                threshold[i] = m.getThreshold();
        }

        if (cudaMalloc((void **)&d_threshold, R.nRows() * sizeof(float)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device\n");
        cublasSetVector(R.nRows(), sizeof(float), threshold, 1, d_threshold, 1);
        delete [] threshold;

        // data structures for the occurrences
        vector<MotifOccurrence> occurrences;
        float *occScore = new float[2*R.nRows() * R.nCols()];
        if (cudaMalloc((void **)&d_occScore, 2 * R.nRows() * R.nCols() * sizeof(float)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device\n");

        int *occIdx = new int[2*R.nRows() * R.nCols()];
        if (cudaMalloc((void **)&d_occIdx, 2 * R.nRows() * R.nCols() * sizeof(int)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device\n");

        int nOcc = 0;
        if (cudaMalloc((void **)&d_nOcc, sizeof(int)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device\n");
        cublasSetVector(1, sizeof(int), &nOcc, 1, d_nOcc, 1);

	const auto matBlock = motifContainer.getMatrixBlock();

        map<int, int> offset_v;

	size_t nOutputTasks = 10;	// FIXME: derive from available number of threads
	vector<future<void> > outputTask(nOutputTasks);
	size_t currOutputTask = 0;

	while (sm.getNextSeqMatrix(seqBatch)) {
		cublasSetVector(4 * (K + overlap) * W, sizeof(float), sm.S.data, 1, d_S, 1);

                for (size_t offset = 0; offset < K; offset++) {

		        for (size_t i = 0; i < matBlock.size(); i++) {
		                size_t start = (i == 0) ? 0 : matBlock[i-1].first;
                		size_t end = matBlock[i].first;

		                int m = end-start;
                		int k = matBlock[i].second;

				cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, R.nCols(), P.nCols(), &alpha, d_P + start, P.nRows(), d_S + 4 * offset, sm.S.nRows(), &beta, d_R + start, R.nRows());
			}
			//cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, R.nRows(), R.nCols(), P.nCols(), &alpha, d_P, P.nRows(), d_S + 4 * offset, sm.S.nRows(), &beta, d_R, R.nRows());
			kernel_wrapper(d_R, R.nRows(), R.nCols(), d_threshold, d_occIdx, d_occScore, d_nOcc);
			int prevOcc = nOcc;
			cublasGetVector(1, sizeof(int), d_nOcc, 1, &nOcc, 1);
			offset_v[offset] = nOcc - prevOcc;

			// don't retreive vector just yet if insufficient results are on the GPU
			if ( (nOcc <= R.nRows() * R.nCols()) && (offset < K-1) )
				continue;

			// get the results
			cublasGetVector(nOcc, sizeof(float), d_occScore, 1, occScore, 1);
			cublasGetVector(nOcc, sizeof(int), d_occIdx, 1, occIdx, 1);
			extractOccurrences2(R.nRows(), offset_v, occIdx, occScore, sm, motifContainer, occurrences);
			offset_v.clear();
			nOcc = 0;
			cublasSetVector(1, sizeof(int), &nOcc, 1, d_nOcc, 1);
                }

		// write the output to disk
		if (outputTask[currOutputTask].valid())
			outputTask[currOutputTask].get();
		outputTask[currOutputTask] = async(launch::async, writeToDisk, speciesID, occurrences, ref(myMutex), ref(os));
		currOutputTask = (currOutputTask + 1) % nOutputTasks;

                totMatches += occurrences.size();
                occurrences.clear();

                cout << "."; cout.flush();
        }

	// wait for all output to be written
	for (size_t i = 0; i < outputTask.size(); i++)
		if (outputTask[i].valid())
			outputTask[i].get();

        delete [] occIdx;
        delete [] occScore;

        cudaFree(d_P);
        cudaFree(d_S);
        cudaFree(d_R);
        cudaFree(d_threshold);
        cudaFree(d_occScore);
        cudaFree(d_occIdx);
        cudaFree(d_nOcc);

        cublasDestroy(handle);
}

void PWMScan::scanPWMCUBLAS(size_t speciesID, const MotifContainer& motifContainer,
                            FastaBatch& seqBatch, std::ostream& os)
{
        int numDevices;
        cudaGetDeviceCount(&numDevices);

        if (numDevices == 0) {
                cerr << "CUDA error: no devices found. Aborting..." << endl;
                return;
        }

        // limit the number of devices to the number of threads specified
        if (numThreads < numDevices)
                numDevices = numThreads;

        cout << "Using " << numDevices << " GPU devices" << endl;

        // start one thread per GPU device
        vector<thread> workerThreads(numDevices);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&PWMScan::scanThreadCUBLAS, this, i, speciesID,
                                          cref(motifContainer),
                                          ref(seqBatch), ref(os));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        cout << endl;
}
#endif

PWMScan::PWMScan(int argc, char ** argv) : simpleMode(false), cudaMode(false),
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
                } else if ((arg == "-c") || (arg == "--cuda")) {
                        cudaMode = true;
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
       // motifContainer.writeMOODSFiles();

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
                        scanPWMNaive(speciesID++, motifContainer, seqBatch, os);
                } else if (cudaMode) {
#ifdef HAVE_CUDA
                        scanPWMCUBLAS(speciesID++, motifContainer, seqBatch, os);
#else
                        cerr << "ERROR: CUDA support not enabled in blstools\n";
                        cerr << "Please recompile with CUDA support enabled" << endl;
#endif
                } else {
                        scanPWMBLAS(speciesID++, motifContainer, seqBatch, os);
                }
        }

        if (!outputFilename.empty())
                ofs.close();

        ofsCutoff.close();

        cout << "\nWrote " << totMatches << " matches to " << outputFilename << ".\n";
}
