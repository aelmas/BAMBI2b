/*=========================================================================
 
 Program:   BAMBI-2b
 Language:  C++
 Date:      $Date: 2013/6/17$
 Version:   $Revision: 2.0 $
 
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.
 
 =========================================================================*/
#include <iostream>
#include <istream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "database.h"
#include <time.h>
#include <algorithm>
#include <string>
#include <cmath>
using namespace std;
/*
 ./BAMBI2b -i sequences.fasta  -o bambi_motifs.txt -P 15000 -LM 19 -UM 21 -LB 1 -UB 5 -LF 0 -UF 3 -r0 1 -r1 300 -d 0.25,0.26,0.25,0.24 -n 5 -s0 1 -s1 1 -s2 1 -sm
 -i  : REQUIRED (default: sequences.fasta)
 -o  : optional (default: bambi_motifs.txt)
 -P  : optional (50 * (max sequence length) )
 -LM : optional (default: 16)
 -UM : optional (default: 26)
 -LB : optional (default: 1)
 -UB : optional (default: UM/2)
 -LF : optional (default: 0)
 -UF : optional (default: UM-2)
 -r0 : optional (default: 1)
 -r1 : optional (default: average length of input sequences)
 -s0 : optional (default: 1 (expected number of instances with palindromic symmetry per input sequence) )
 -s1 : optional (default: 1 (expected number of instances with inverted repeat symmetry per input sequence) )
 -s2 : optional (default: 1 (expected number of instances with direct repeat symmetry per input sequence) )
 -d  : optional (default: 0.25,0.25,0.25,0.25)
 -n  : optional (default: 1)
 -sm : optional (default: false)
 -op : optional (default: 1)
 */

int main(int argc, char * argv[])
{
    // initialize input variables
    std::string inputFileName = "sequences.fasta"; //default input file name
    for(int i = 0; i < argc-1; i++){
        if(string(argv[i])=="-i"){
            inputFileName = string(argv[i+1]).c_str();
            break;
        }
    }
    // ---------------------------------------------------
    // Read fasta file
    // ---------------------------------------------------
    ifstream inputFileStream (inputFileName.c_str());
    int totalsize(0);
    int nSeq(0);
    int maxSeq(0);
    if (inputFileStream.is_open()){
        bool flagSeq(false);
        string dna;
        dna.clear();
        while ( inputFileStream.good() ){
            string line;
            getline (inputFileStream,line);
            if (line.find(">")==0){
                if (dna.empty()) {
                    flagSeq = true;
                    // read next line
                } else{
                    if (flagSeq){
                        // write
                        if(int(dna.size())>maxSeq){
                            maxSeq = int(dna.size());
                        }
                        totalsize += int(dna.size());
                        nSeq ++;
                        dna.clear();
                        flagSeq = true;
                        // read next line
                    } else{
                        flagSeq = true;
                        // read next line
                    }
                }
            } else if (line.find(">")>0 && line.find(">")<line.size()){
                // empty, read next line
            } else if (line.empty()){
                // empty, read next line
            } else {
                if (flagSeq) {
                    // deblank line
                    string dbline;
                    dbline.clear();
                    for (int i=0; i<int(line.size()); i++){
                        switch (line[i]) {
                            case 'a':
                            case 'c':
                            case 'g':
                            case 't':
                            case 'A':
                            case 'C':
                            case 'G':
                            case 'T':
                                dbline += line[i];
                                break;
                            default:
                                break;
                        }
                    }
                    // append current line
                    dna += dbline;
                }
                // read next line
            }
        }
        if (dna.empty()){
            // continue
        } else {
            if (flagSeq){
                // write
                if(int(dna.size())>maxSeq){
                    maxSeq = int(dna.size());
                }
                totalsize += int(dna.size());
                nSeq ++;
                dna.clear();
                // continue
            }
        }
    }
    inputFileStream.close();
    int avgSeq(0);
    if(nSeq>0){
        avgSeq = int(totalsize/float(nSeq));
    }
    // ---------------------------------------------------
    std::string outputFileName = "bambi_motifs.txt";
    int P  = 50*maxSeq;
    int LM = 16;
    int UM = 26;
    int LB = 0;
    int UB = int(floor(UM/2));
    int LF = 0;
    int UF = UM-2;
    int r0 = 1;
    int r1 = avgSeq;
    int s0 = 1;
    int s1 = 1;
    int s2 = 1;
    int op = 1;
    float d[] = {0.25, 0.25, 0.25, 0.25};
    int n = 1;
    string strCommand;
    // parse command line arguments
    for(int i = 0; i < argc-1; i++){
        strCommand += argv[i];
        strCommand += ' ';
        if(string(argv[i])=="-o"){
            outputFileName[string(argv[i+1]).size()]=0;
            outputFileName = string(argv[i+1]).c_str();
        }else if(string(argv[i])=="-P"){
            P = atoi(argv[i+1]);
        }else if(string(argv[i])=="-LM"){
            LM = atoi(argv[i+1]);
        }else if(string(argv[i])=="-UM"){
            UM = atoi(argv[i+1]);
        }else if(string(argv[i])=="-LB"){
            LB = atoi(argv[i+1]);
        }else if(string(argv[i])=="-UB"){
            UB = atoi(argv[i+1]);
        }else if(string(argv[i])=="-LF"){
            LF = atoi(argv[i+1]);
        }else if(string(argv[i])=="-UF"){
            UF = atoi(argv[i+1]);
        }else if(string(argv[i])=="-r0"){
            r0 = atoi(argv[i+1]);
        }else if(string(argv[i])=="-r1"){
            r1 = atoi(argv[i+1]);
        }else if(string(argv[i])=="-s0"){
            s0 = atoi(argv[i+1]);
        }else if(string(argv[i])=="-s1"){
            s1 = atoi(argv[i+1]);
        }else if(string(argv[i])=="-s2"){
            s2 = atoi(argv[i+1]);
        }else if(string(argv[i])=="-op"){
            op = atoi(argv[i+1]);
        }else if(string(argv[i])=="-d"){
            int j=0;
            string val;
            string temp = string(argv[i+1]).c_str();
            for(int k=0; k<int(temp.size()); k++){
                if (temp[k] == ','){
                    d[j] = (float)atof(val.c_str());
                    val = " ";
                    j++;
                }else{
                    val += temp[k];
                }
            }
            d[3] = (float)atof(val.c_str());
        }else if(string(argv[i])=="-n"){
            n = atoi(argv[i+1]);
        }
    }
    strCommand += argv[argc-1];
    if(strCommand.find("-UM")!=std::string::npos){
        if(strCommand.find("-UB")==std::string::npos){
            UB = int(floor(UM/2));
        }
    }
    if(strCommand.find("-UM")!=std::string::npos || strCommand.find("-LB")!=std::string::npos){
        if(strCommand.find("-UF")==std::string::npos){
            UF = UM-2;
        }
    }
    if(strCommand.find("-h")!=std::string::npos){ // || argc<2){
        cout << "\nMissing input arguments. Example Usage: \n \n";
        cout << "./BAMBI2b -i sequences.fasta  -o bambi_motifs.txt -P 15000 -LM 19 -UM 21 -LB 1 -UB 5 -LF 0 -UF 3 -r0 1 -r1 300 -d 0.25,0.26,0.25,0.24 -n 5 -s0 1 -s1 1 -s2 1 -sm -op 1 \n \n";
        cout << "Options: \n";
        cout << "-i  : REQUIRED (default: sequences.fasta) \n";
        cout << "-o  : optional (default: bambi_motifs.txt) \n";
        cout << "-P  : optional (50 * (max sequence length) ) \n";
        cout << "-LM : optional (default: 16) \n";
        cout << "-UM : optional (default: 26) \n";
        cout << "-LB : optional (default: 1) \n";
        cout << "-UB : optional (default: UM/2) \n";
        cout << "-LF : optional (default: 0) \n";
        cout << "-UF : optional (default: UM - 2) \n";
        cout << "-r0 : optional (default: 1) \n";
        cout << "-r1 : optional (default: average length of input sequences) \n";
        cout << "-s0 : optional (default: 1 (expected number of instances with palindromic symmetry per input sequence) ) \n";
        cout << "-s1 : optional (default: 1 (expected number of instances with inverted repeat symmetry per input sequence) ) \n";
        cout << "-s2 : optional (default: 1 (expected number of instances with direct repeat symmetry per input sequence) ) \n";
        cout << "-d  : optional (default: 0.25,0.25,0.25,0.25) \n";
        cout << "-n  : optional (default: 1) \n";
        cout << "-sm : optional (default: false) \n";
        cout << "-op : optional (default: 1) \n \n";
        return EXIT_FAILURE;
    }
    ofstream fileOut;
    fileOut.open(outputFileName.c_str());
    fileOut << "========================================\n";
    fileOut << "BAMBI - Motif Finding Algorithm\n";
    fileOut << "========================================\n";
    fileOut << "BAMBI version 2.0 (BAMBI-2b: Release date: Oct 13 2013)\n";
    fileOut << "\nFor further information please visit: http://genomics.lbl.gov/BAMBI/";
    fileOut << "\n\n";
    fileOut << "========================================\n";
    fileOut << "COMMAND LINE SUMMARY\n";
    fileOut << "========================================\n";
    fileOut << "command: " << strCommand << endl;
    fileOut << "\n";
    fileOut << "----------------------------------------\n";
    fileOut << "Input file name: " << inputFileName << endl;
    fileOut << "Number of sequences: " << nSeq << endl;
    fileOut << "Maximum sequence length: " << maxSeq << endl;
    fileOut << "Average sequence length: " << avgSeq << endl;
    // ---------------------------------------------------
    // Re-initialize parameters for the example data
    // ---------------------------------------------------
    if (nSeq==0){
        nSeq = 10;
        maxSeq = 100;
        avgSeq = 100;
        if (P==0){P = 50*maxSeq;}
        if (r1==0){r1 = avgSeq;}
        fileOut << "----------------------------------------\n";
        fileOut << "Using the example data: " << endl;
        fileOut << "Number of sequences: " << nSeq << endl;
        fileOut << "Maximum sequence length: " << maxSeq << endl;
        fileOut << "Average sequence length: " << avgSeq << endl;
    }
    // ---------------------------------------------------
    ofstream classesOut;
    classesOut.open("classes.txt");
    classesOut << "M\tB\tlF\trF\tTemplate\n";
    int nclasses(0);
    bool noblock(false);
    int lambda;
    for(int m=0; m<int(UM-LM+1); m++){
        int b=0;
        while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
            noblock = false;
            for(int f=0; f<int(UF-LF+1); f++){
//                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                lambda = 0;
                if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                // class(M,B,F,lambda)
                nclasses++;
                vector<string> seq(LM+m);
                for (int i=0; i < LM+m; i++) {
                    if (i < LF+f) { //left flank
                        seq[i] = 'F';
                    } else if (i <= LF+f+LB+b-1){ //1st block
                        seq[i] = 'B';
                    } else if (i < LM+m-1-(LF+f+lambda)-(LB+b)+1){ //gap
                        seq[i] = '/';
                    } else if (i <= LM+m-1-(LF+f+lambda)){ //second block
                        seq[i] = 'B';
                    } else { //right flank
                        seq[i] = 'F';
                    }
                }
                classesOut << LM+m << "\t" << LB+b << "\t" << LF+f << "\t" << LF+f+lambda << "\t";
                for (int i=0; i<int(seq.size()); i++) {
                    classesOut << seq[i];
                }
                classesOut << endl;
                // class(M,B,F,lambda)
                if (LB+b==0){ noblock=true;}
//                } //end lambda
            } //end f
            b++;
        } // end b
    } // end m
    cout << "\nNumber of classes: " << nclasses << endl;
    P = max(P,nclasses);
    // ---------------------------------------------------
    fileOut << "----------------------------------------\n";
    fileOut << "Number of runs (motifs): " << n << endl;
    if(strCommand.find("-sm")!=std::string::npos){
        fileOut << "Search minor sites: true\n";
    }else{
        fileOut << "Search minor sites: false\n";
    }
    fileOut << "Background nucleotide distribution [A,C,G,T]: [" << d[0] << "," << d[1] << "," << d[2] << "," << d[3] << "]\n";
    fileOut << "LM: " << LM << ", UM: " << UM << endl;
    fileOut << "LB: " << LB << ", UB: " << UB << endl;
    fileOut << "LF: " << LF << ", UF: " << UF << endl;
    fileOut << "Number of particles: " << P << endl;
    fileOut << "r0: " << r0 << ", r1: " << r1 << endl;
    fileOut << "s0: " << s0 << ", s1: " << s1 << ", s2: " << s2 << endl;
    fileOut << "op: " << op << endl;
    fileOut << "Number of classes: " << nclasses << endl;
    fileOut << "Number of particles per class: " << int(ceil(float(P)/float(nclasses))) << endl;
    fileOut << "----------------------------------------\n";
    fileOut << "\n";
    for(int a=0; a<n; a++){
        std::cout << "Run #: " << a+1 <<endl;
        fileOut << endl;
        fileOut << "========================================\n";
        fileOut << "MOTIF " << a+1 << ": \n";
        fileOut << "========================================\n";
        ifstream myfile (inputFileName.c_str());
        Database RealData(myfile);
        myfile.close();
        RealData.SetNumberOfParticles(P);
        RealData.SetRho0(r0);
        RealData.SetRho1(r1);
        RealData.SetSigma0(s0);
        RealData.SetSigma1(s1);
        RealData.SetSigma2(s2);
        RealData.SetSymOption(op);
        RealData.SetNullDist(d);
        RealData.SetNclasses(nclasses);
        RealData.SetLM(LM);
        RealData.SetUM(UM);
        RealData.SetLB(LB);
        RealData.SetUB(UB);
        RealData.SetLF(LF);
        RealData.SetUF(UF);
        //
        // main program
        //
        /* initialize random seed: */
        srand ( int(time(NULL)) );
        RealData.uninformativePWM();
        RealData.ComputeParticleFilter(fileOut);
        if(strCommand.find("-sm")!=std::string::npos){
            // search minor sites
            RealData.DetectionZeroAgainstMany();
            RealData.DetectionOneAgainstMore();
            RealData.ComputeFinalSolution(fileOut);
        }
    }
    cout << "End of the program. Please see the output file " << outputFileName << endl;
    return 0;
}
