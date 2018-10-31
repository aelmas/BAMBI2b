#include "database.h"
#include <fstream>
#include <iostream>
#include <istream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <math.h>
#define w1 double(1)
#define w2 double(1)
#define tau double(0.01)
#define cB double(1)
#define cW double(1)
#define verbose int(0)
// Detection and Estimation Test
#define Pfa 0.01
#define NUMITERATIONS 100
#define initcount double(1)

Database::Database()
{
    NumSequences = 0;
}

Database::Database(ifstream & myfile)
{
    // ---------------------------------------------------
    // Read fasta file
    // ---------------------------------------------------
    cout << "Reading input file...\n";
    NumSequences = 0;
    if (myfile.is_open()){
        bool flagSeq(false);
        string dna;
        dna.clear();
        while ( myfile.good() ){
            string line;
            getline (myfile,line);
            if (line.find(">")==0){
                if (dna.empty()) {
                    flagSeq = true;
                    // read next line
                } else{
                    if (flagSeq){
                        // write
                        cout << "Sequence " << NumSequences+1 << endl;
                        cout << dna << endl;
                        Sequences.push_back(dna);
                        NumSequences ++;
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
                cout << "Sequence " << NumSequences+1 << endl;
                cout << dna << endl;
                Sequences.push_back(dna);
                NumSequences ++;
                dna.clear();
                // continue
            }
        }
    }
    // ---------------------------------------------------
    // Read example data
    // ---------------------------------------------------
    if (NumSequences==0){
        cout << "\nFile not found! Please provide an input file in proper format (fasta) ... \n" << endl;
        cout << "\nUsing example data \n" << endl;
        Sequences.push_back("ATGCTAATAGATAGAGCTGCTTGGGCTTACGCCCGTTAAAGCTGTGTATAAATACAGTAAAGAAAGGAAACATGTAAAACAACCTTGCCAAGATAAAGCA");//1
        NumSequences++;
        Sequences.push_back("CTTACTGTTTTGATAAACAGATAGAATAATGTTGTTTTTTTTTAGCGGTACATGTTCACTTCCATAATTGACATTCAGTTATGTTACATTTCTTATGCAA");//2
        NumSequences++;
        Sequences.push_back("CGCTTTACTGTCAATACTGGAGCTGTATTTCTATCCAGTAACCTTCACCTCGCTATATATTTGTGTTTGAACAATACATGATATGGATTTTCTTTTAATC");//3
        NumSequences++;
        Sequences.push_back("TATTTCCCTGTTTTTTCTGATTCTCCAGGTACTGTAAATTTATACAGTAACGGTACTGTTAGTTGTAAACTATGCATACAAACTCAAAGGTGTAAAGACC");//4
        NumSequences++;
        Sequences.push_back("TGTATTTCTAGCCTTATGCTGCTATTTAAATAAATAAACTTAAGCTACTGAATTTTCATACAGTATATTGACAACTAAAGACGATCATAGTAAAAGAGTG");//5
        NumSequences++;
        Sequences.push_back("ATCGTCATCTTTAAATTATAATATTAAAAATAGCGTATGGCAGTCTTGTACTGGATGTTTATCCAGTAGCTCGAATACACATGGTCTGTTTCATATCCTT");//6
        NumSequences++;
        Sequences.push_back("ATTCTTTTTTTTTAGGTGCTGTTTATATATACAGGAAAATGAATCAATTTACTGGTTGATAGTCCTTTTTTTTTCTAAAGTATCGATGTCAATTTTTTAT");//7
        NumSequences++;
        Sequences.push_back("TTTTTAATACTGTCCATTACTGTATAATCATACAGTTGCGCCAAAATTTCTGTAAGTGCTATGGGACTTGTTTAATACTAGAAGCAGTTATTTTTGATAA");//8
        NumSequences++;
        Sequences.push_back("GAAAGAGAGAGAAATGTGTAGCGCATGTATTTTGACCGGATCAGCATTAGTGTATGACTGTATATTTTAACAGTTCACAGACCTTTTCAAAAGATTCCCG");//9
        NumSequences++;
        Sequences.push_back("GCATAAGAGCTTGCGCAAATTTTGCAGCCTGGATGTTTATACAGCTCTGGTTGAAACCACTGTGCGATAAAACTAATGAAGTTAACTGATCGTGAAACTT");//10
        NumSequences++;
    }
    oSequences = Sequences;
}

void Database::insertSequence(string NewSequence)
{
    Sequences.push_back(NewSequence);
    NumSequences ++;
}

int Database::getSize()
{
    return NumSequences;
}

Database::~Database()
{
    //dtor
}

void Database::showPWM()
{
    for (int c=0; c<4; c++)
    {
        for (int j=0; j<int(PWM.size()); j++)
        {
            cout << PWM[j][c] << " ";
        }
        cout <<endl;
    }
}

void Database::uninformativePWM(){
    bool noblock(false);
    int lambda;
    vector<double> column_vec;
    for (int i=0; i<4; i++) {
        column_vec.push_back(double(NullDist[i]));
    }
    for(int m=0; m<int(UM-LM+1); m++){
        int b=0;
        vector<vector<vector<vector< double> > > > tempBG;
        while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
            noblock = false;
            vector<vector<vector< double> > > tempG;
            for(int f=0; f<int(UF-LF+1); f++){
                //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                lambda = 0;
                if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                // class(M,B,F,lambda)
                vector<vector<double> > tempPWM;
                for(int ir=0; ir<4; ir++){
                    vector<double> ColDist((LM+m),NullDist[ir]);
                    tempPWM.push_back(ColDist);
                }
                // class(M,B,F,lambda)
                if (LB+b==0){ noblock=true;}
                //                } //end lambda
                tempG.push_back(tempPWM);
            }
            tempBG.push_back(tempG);
            b++;
        }
        ALPHA.push_back(tempBG);
    }
    vector< vector<double> > tempPWM ((LM+0), column_vec);
    PWM = tempPWM;
    for(int c=0; c<(LM+0); c++)
    {
        float normalizing(0);
        for (int j=0; j<4; j++)
        {
            normalizing+=ALPHA[0][0][0][j][c];
        }
        for (int j=0; j<4; j++)
        {
            PWM[c][j]=double(ALPHA[0][0][0][j][c])/normalizing;
        }
    }
}

void Database::DetectionZeroAgainstMany()
{
    vector<bool> tIndicator(oSequences.size(), false);
    IndicatorAtLeastOneMotif = tIndicator;
    IndicatorAtLeastTwoMotif = tIndicator;
    /////////////////////////////////////////////////////////////////////////////
    ///////////// Perform the detection and estimation task!!!! /////////////////
    /////////////////////////////////////////////////////////////////////////////
    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++\n";
    cout << "Detection And Estimation Test"<<std::endl;
    cout << "++++++++++++++++++++++++++++++++++++++++\n\n";
    cout << "Running DetectionZeroAgainstMany"<<endl;
    for (int row = 0; row<4; row++)
    {
        for (int column=0; column<(LM+estiM); column++)
        {
        }
    }
    for (int row = 0; row<4; row++)
    {
        for (int column=0; column<(LM+estiM); column++)
        {
        }
    }
    // fix Alpha and Beta
    int MotifLength((LM+estiM));
    vector<float> Tempfix(MotifLength,1);
    vector<vector<float> > AlphaDirichlet(4, Tempfix);
    vector<vector<float> > fixAlpha(3, Tempfix);
    vector<vector<float> > fixBeta(3, Tempfix);
    for (int row = 0; row<4; row++)
    {
        for (int column=0; column<MotifLength; column++)
        {
            AlphaDirichlet[row][column] = ALPHA[estiM][estiB][estiF][row][column];
        }
    }
    for (int column=0; column<MotifLength; column++)
    {
        fixBeta[2][column]=AlphaDirichlet[2+1][column];
        fixAlpha[2][column]=AlphaDirichlet[2][column];
    }
    for (int row = 1; row>-1; row--)
    {
        for (int column=0; column<MotifLength; column++)
        {
            fixBeta[row][column]=AlphaDirichlet[row+1][column]+fixBeta[row+1][column];
            fixAlpha[row][column]=AlphaDirichlet[row][column];
        }
    }
    ////////////////////////////////////////////////////////////////////////
    // Find the test to decide whether there is one motif or no motif. ///
    ///////////////////////////////////////////////////////////////////////
    vector<int> tempxy2(2,0);
    vector<vector<int> > BestPositioning (oSequences.size(),tempxy2);
    vector<long double> likelihoodratios(oSequences.size(),0);
    for(int iteratorS=0; iteratorS<int(oSequences.size()); iteratorS++)
    {
        string CurrentSequence = oSequences[iteratorS];
        
        float UnderH1(0);
        for (int pos=0; pos< int(CurrentSequence.size()-MotifLength+1); pos++)
        {
            float Posterior;
            Posterior = ComputePosterior1Motif(CurrentSequence, pos, fixAlpha, fixBeta, MotifLength);
            UnderH1 += Posterior;
        }
        likelihoodratios[iteratorS] = UnderH1;
    }
    vector<long double> Threshold(oSequences.size(),0);
    int numIterations(NUMITERATIONS);
    vector<double>p0(4,0); p0[0]=(float)NullDist[0]; p0[1]=(float)NullDist[1]; p0[2]=(float)NullDist[2]; p0[3]=(float)NullDist[3];
    vector<float> nn(4,0);
    nn[0]=0;
    nn[1]=1;
    nn[2]=2;
    nn[3]=3;
    for (int indSeq=0; indSeq<int(oSequences.size()); indSeq++)
    {
        vector<long double> ratios(numIterations,0);
        for (int iter=0; iter<numIterations; iter++)
        {
            string Sequence;
            Sequence.clear();
            for (int posSeq = 0; posSeq<int(oSequences[indSeq].size()); posSeq++)
            {
                int currentNN(0);
                currentNN = int(pmfsim(nn, p0));
                switch(currentNN)
                {
                    case 0:
                        Sequence+="A";
                        break;
                    case 1:
                        Sequence+="C";
                        break;
                    case 2:
                        Sequence+="G";
                        break;
                    case 3:
                        Sequence+="T";
                        break;
                }
            }
            float UnderH0(0);
            for (int pos=0; pos< int(Sequence.size()-MotifLength+1); pos++)
            {
                float Posterior;
                Posterior = ComputePosterior1Motif(Sequence, pos, fixAlpha, fixBeta, MotifLength);
                UnderH0 += Posterior;
            }
            ratios[iter] = UnderH0;
            
        }
        cout << "Finding the threshold..."<<endl;
        sort(ratios.begin(), ratios.end());
        for(int jo=0; jo<int(ratios.size()); jo++)
        {
        }
        int tempPos(0);
        tempPos = numIterations-int(numIterations*Pfa);
        Threshold[indSeq] = (ratios[tempPos]+ratios[tempPos+1])/2;
        cout << "For sequence " << indSeq+1<<", the threshold is "<<Threshold[indSeq]<<endl;
        ratios.clear();
    }
    
    for(int iteratorS=0; iteratorS<int(oSequences.size()); iteratorS++)
    {
        string CurrentSequence = oSequences[iteratorS];
        
        if (likelihoodratios[iteratorS]>Threshold[iteratorS] && likelihoodratios[iteratorS]!=0)
        {
            cout << iteratorS+1<<"  has an instance of the motif"<<endl;
            IndicatorAtLeastOneMotif[iteratorS]=true;
        }
        else
        {
            cout << iteratorS+1<< " does NOT have an instance of the motif" <<endl;
            IndicatorAtLeastOneMotif[iteratorS]=false;
        }
    }
}

void Database::DetectionOneAgainstMore()
{
    vector<bool> tIndicator(oSequences.size(), false);
    IndicatorAtLeastTwoMotif = tIndicator;
    cout << "Running DetectionOneAgainstMore"<<endl;
    // fix Alpha and Beta
    int MotifLength((LM+estiM));
    vector<float> Tempfix(MotifLength,1);
    vector<vector<float> > AlphaDirichlet(4, Tempfix);
    vector<vector<float> > fixAlpha(3, Tempfix);
    vector<vector<float> > fixBeta(3, Tempfix);
    for (int row = 0; row<4; row++)
    {
        for (int column=0; column<MotifLength; column++)
        {
            AlphaDirichlet[row][column] = ALPHA[estiM][estiB][estiF][row][column];
        }
    }
    for (int column=0; column<MotifLength; column++)
    {
        fixBeta[2][column]=AlphaDirichlet[2+1][column];
        fixAlpha[2][column]=AlphaDirichlet[2][column];
    }
    for (int row = 1; row>-1; row--)
    {
        for (int column=0; column<MotifLength; column++)
        {
            fixBeta[row][column]=AlphaDirichlet[row+1][column]+fixBeta[row+1][column];
            fixAlpha[row][column]=AlphaDirichlet[row][column];
        }
    }
    ////////////////////////////////////////////////////////////////////////
    // Find the test to decide whether there is one motif or more.      ///
    ///////////////////////////////////////////////////////////////////////
    cout << "Finding the likelihoods NOW!..."<<endl;
    vector<long double> likelihoodratios(oSequences.size(),0);
    for(int iteratorS=0; iteratorS<int(oSequences.size()); iteratorS++)
    {
        cout <<"Processing sequence "<<iteratorS+1<<" of length "<<oSequences[iteratorS].size()<<endl;
        if(IndicatorAtLeastOneMotif[iteratorS]!=true)
        {
            IndicatorAtLeastTwoMotif[iteratorS]=false;
            likelihoodratios[iteratorS] = 0;
        }
        else
        {
            string CurrentSequence = oSequences[iteratorS];
            float UnderH1(0);
            for (int pos=0; pos< int(CurrentSequence.size()-MotifLength+1); pos++)
            {
                for (int pos2=pos+MotifLength; pos2< int(CurrentSequence.size()-MotifLength+1); pos2++)
                {
                    float Posterior;
                    Posterior = ComputePosterior2Motifs(CurrentSequence.substr(pos,MotifLength), CurrentSequence.substr(pos2,MotifLength), fixAlpha, fixBeta, MotifLength);
                    UnderH1 += Posterior;
                }
            }
            float UnderH0(0);
            for (int pos=0; pos< int(CurrentSequence.size()-MotifLength+1); pos++)
            {
                float Posterior;
                Posterior = ComputePosterior1Motif(CurrentSequence, pos, fixAlpha, fixBeta, MotifLength);
                UnderH0 += Posterior;
            }
            likelihoodratios[iteratorS] = UnderH1/UnderH0;
        }
    }
    vector<long double> Threshold(oSequences.size(),0);
    int numIterations(NUMITERATIONS);
    vector<double>p0(4,0); p0[0]=(float)NullDist[0]; p0[1]=(float)NullDist[1]; p0[2]=(float)NullDist[2]; p0[3]=(float)NullDist[3];
    vector<float> nn(4,0);
    nn[0]=0;
    nn[1]=1;
    nn[2]=2;
    nn[3]=3;
    cout << "generating sequences..."<<endl;
    for (int indSeq=0; indSeq<int(oSequences.size()); indSeq++)
    {
        cout << "Procesing sequence "<<indSeq+1<<", of length "<<oSequences[indSeq].size()<<endl;
        if(int(oSequences[indSeq].size())>=int((LM+estiM)))
        {
            vector<long double> ratios(numIterations,0);
            for (int iter=0; iter<numIterations; iter++)
            {
                // generate sequence with H0
                string Sequence;
                Sequence.clear();
                for (int posSeq = 0; posSeq<int(oSequences[indSeq].size()); posSeq++)
                {
                    int currentNN(0);
                    currentNN = int(pmfsim(nn, p0));
                    switch(currentNN)
                    {
                        case 0:
                            Sequence+="A";
                            break;
                        case 1:
                            Sequence+="C";
                            break;
                        case 2:
                            Sequence+="G";
                            break;
                        case 3:
                            Sequence+="T";
                            break;
                    }
                }
                // put a motif in the beginning (the position is not important)
                for (int posSeq = 0; posSeq<(LM+estiM); posSeq++)
                {
                    int currentNN(0);
                    vector<double> distr(4,0);
                    for (int c=0; c<4; c++)
                    {
                        distr[c]=PWM[posSeq][c];
                    }
                    currentNN = int(pmfsim(nn, distr));
                    switch(currentNN)
                    {
                        case 0:
                            Sequence[posSeq]='A';
                            break;
                        case 1:
                            Sequence[posSeq]='C';
                            break;
                        case 2:
                            Sequence[posSeq]='G';
                            break;
                        case 3:
                            Sequence[posSeq]='T';
                            break;
                    }
                }
                float UnderH1(0);
                for (int pos=0; pos< int(Sequence.size()-MotifLength+1); pos++)
                {
                    for (int pos2=pos+MotifLength; pos2< int(Sequence.size()-MotifLength+1); pos2++)
                    {
                        float Posterior;
                        Posterior = ComputePosterior2Motifs(Sequence.substr(pos,MotifLength), Sequence.substr(pos2,MotifLength), fixAlpha, fixBeta, MotifLength);
                        UnderH1 += Posterior;
                    }
                }
                float UnderH0(0);
                for (int pos=0; pos< int(Sequence.size()-MotifLength+1); pos++)
                {
                    float Posterior;
                    Posterior = ComputePosterior1Motif(Sequence, pos, fixAlpha, fixBeta, MotifLength);
                    UnderH0 += Posterior;
                }
                ratios[iter] = UnderH1/UnderH0;
            }
            cout << "Finding the threshold..."<<endl;
            sort(ratios.begin(), ratios.end());
            for(int jo=0; jo<int(ratios.size()); jo++)
            {
                //                cout << ratios[jo]<<" ";
            }
            cout <<endl;
            int tempPos(0);
            tempPos = numIterations-int(numIterations*Pfa);
            Threshold[indSeq] = (ratios[tempPos]+ratios[tempPos+1])/2;
            cout << "For sequence "<<indSeq+1<<", the threshold is "<<Threshold[indSeq]<<endl;
        }
        else
        {
            Threshold[indSeq] = 0;
            cout << "For sequence "<<indSeq+1<<", the threshold is "<<Threshold[indSeq]<<endl;
        }
    }//end of indSeq
    for(int iteratorS=0; iteratorS<int(oSequences.size()); iteratorS++)
    {
        if(IndicatorAtLeastOneMotif[iteratorS]!=true)
        {
            cout << iteratorS+1<<"  has 0 instances of the motif"<<endl;
            IndicatorAtLeastOneMotif[iteratorS]=false;
            IndicatorAtLeastTwoMotif[iteratorS]=false;
        }
        else
        {
            string CurrentSequence = oSequences[iteratorS];
            if (likelihoodratios[iteratorS]>Threshold[iteratorS] && likelihoodratios[iteratorS]!=0)
            {
                cout << iteratorS+1<<"  has 2 instances of the motif"<<endl;
                IndicatorAtLeastOneMotif[iteratorS]=true;
                IndicatorAtLeastTwoMotif[iteratorS]=true;
            }
            else
            {
                cout << iteratorS+1<< " has 1 instance of the motif" <<endl;
                IndicatorAtLeastOneMotif[iteratorS]=true;
                IndicatorAtLeastTwoMotif[iteratorS]=false;
            }
        }
    }
}


void Database::ComputeFinalSolution(ofstream & fileOut)
{
    // fix Alpha and Beta
    int MotifLength((LM+estiM));
    vector<float> Tempfix(MotifLength,1);
    vector<vector<float> > AlphaDirichlet(4, Tempfix);
    vector<vector<float> > fixAlpha(3, Tempfix);
    vector<vector<float> > fixBeta(3, Tempfix);
    // Plotting AlphaDirichlet
    for (int row = 0; row<4; row++)
    {
        for (int column=0; column<MotifLength; column++)
        {
            AlphaDirichlet[row][column] = ALPHA[estiM][estiB][estiF][row][column];
        }
    }
    
    for (int column=0; column<MotifLength; column++)
    {
        fixBeta[2][column]=AlphaDirichlet[2+1][column];
        fixAlpha[2][column]=AlphaDirichlet[2][column];
    }
    for (int row = 1; row>-1; row--)
    {
        for (int column=0; column<MotifLength; column++)
        {
            fixBeta[row][column]=AlphaDirichlet[row+1][column]+fixBeta[row+1][column];
            fixAlpha[row][column]=AlphaDirichlet[row][column];
        }
    }
    fileOut << "Sequence No \t Minor Instance Start \t Minor Instance Sequence"<<endl;
    fileOut << "----------------------------------------\n";
    vector<int> tempxy(2,0);
    vector<vector<int> > tInferredPositions(oSequences.size(),tempxy);
    for(int iteratorS=0; iteratorS<int(oSequences.size()); iteratorS++)
    {
        string CurrentSequence;
        CurrentSequence = oSequences[iteratorS];
        if(IndicatorAtLeastOneMotif[iteratorS]!=true)
        {
            tempxy[0]=0;
            tempxy[1]=0;
            tInferredPositions[iteratorS]=tempxy;
        }
        else
        {
            if(IndicatorAtLeastTwoMotif[iteratorS]!=true)
            {
                // USE ML estimator
                float MaxML(0);
                for (int pos=0; pos< int(CurrentSequence.size()-MotifLength+1); pos++)
                {
                    float Posterior;
                    Posterior = ComputePosterior1Motif(CurrentSequence, pos, fixAlpha, fixBeta, MotifLength);
                    if(Posterior>MaxML)
                    {
                        MaxML = Posterior;
                        tempxy[0]=pos+1;
                        tempxy[1]=0;
                    }
                }
                tInferredPositions[iteratorS]=tempxy;
            }
            else
            {
                // USE ML estimator
                float MaxML(0);
                for (int pos=0; pos< int(CurrentSequence.size()-MotifLength+1); pos++)
                {
                    for (int pos2=pos+MotifLength; pos2< int(CurrentSequence.size()-MotifLength+1); pos2++)
                    {
                        float Posterior;
                        Posterior = ComputePosterior2Motifs(CurrentSequence.substr(pos,MotifLength), CurrentSequence.substr(pos2,MotifLength), fixAlpha, fixBeta, MotifLength);
                        if(Posterior>MaxML)
                        {
                            MaxML = Posterior;
                            tempxy[0]=pos+1;
                            tempxy[1]=pos2+1;
                        }
                    }
                }
                tInferredPositions[iteratorS]=tempxy;
            }
        }
        cout << "ML estimated site locations at " << iteratorS << "th sequence\n";
        cout << tInferredPositions[iteratorS][0]<<" "<<tInferredPositions[iteratorS][1]<<endl;
        if(tInferredPositions[iteratorS][0]!=InferredPositions1[iteratorS] && tInferredPositions[iteratorS][0]!=0){
            fileOut << iteratorS+1 << "\t" << tInferredPositions[iteratorS][0] << "\t" <<
            CurrentSequence.substr(tInferredPositions[iteratorS][0]-1,(LM+estiM))<< endl;
        }
        if(tInferredPositions[iteratorS][1]!=InferredPositions1[iteratorS] && tInferredPositions[iteratorS][1]!=0){
            fileOut << iteratorS+1 << "\t" << tInferredPositions[iteratorS][1] << "\t" <<
            CurrentSequence.substr(tInferredPositions[iteratorS][1]-1,(LM+estiM))<< endl;
        }
        
    }//end iteratorS
    fileOut << "----------------------------------------\n";
    InferredPositions=tInferredPositions;
}

void Database::ComputeParticleFilter(ofstream & fileOut)
{
    srand ( int(time(NULL)) );
    // Main algorithm
    int NumParticles;
    NumParticles = NumberOfParticles;
    // Set the prior for the PWM
    SetEstiM(0); SetEstiB(0); SetEstiF(0);
    cout << "Average number of particles per class: " << int(ceil(float(NumParticles)/float(nclasses))) << endl;
    //
    // Initialize variables: Alpha, rho0, rho1, sigma0, sigma1, weights
    //
    int MotifLength;
    vector<vector<vector<vector<vector<vector<double> > > > > >Alpha;
    vector<vector<vector<vector<int> > > > rho0;
    vector<vector<vector<vector<int> > > > rho1;
    vector<vector<vector<vector<int> > > > sigma0;
    vector<vector<vector<vector<int> > > > sigma1;
    vector<vector<vector<vector<int> > > > sigma2;
    vector<vector<vector<vector<long double> > > > weights;
    vector<vector<vector<vector<long double> > > > weightsPicked;
    vector<vector<vector<vector<vector<int> > > > > EstimatedStatesSeq;
    vector<vector<vector<vector<vector<int> > > > > Trace;
    bool noblock(false);
    int lambda;
    for(int m=0; m<int(UM-LM+1); m++){
        int b=0;
        vector<double> column_vec((LM+m),initcount);
        vector< vector<double> > transpwm(4, column_vec);
        vector<vector<vector<vector<vector< double > > > > > tempBG;
        vector<vector<vector<int> > > temprho0BG;
        vector<vector<vector<int> > > temprho1BG;
        vector<vector<vector<int> > > tempsigma0BG;
        vector<vector<vector<int> > > tempsigma1BG;
        vector<vector<vector<int> > > tempsigma2BG;
        vector<vector<vector<long double> > > tempweightsBG;
        while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
            noblock = false;
            vector<vector<vector<vector< double> > > > tempG;
            vector<vector<int> > temprho0G;
            vector<vector<int> > temprho1G;
            vector<vector<int> > tempsigma0G;
            vector<vector<int> > tempsigma1G;
            vector<vector<int> > tempsigma2G;
            vector<vector<long double> > tempweightsG;
            for(int f=0; f<int(UF-LF+1); f++){
                //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                lambda = 0;
                if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                // class(M,B,F,lambda)
                tempG.push_back(vector<vector<vector<double> > > (int(ceil(float(NumParticles)/float(nclasses))),transpwm));
                temprho0G.push_back(vector<int> (int(ceil(float(NumParticles)/float(nclasses))),InitialRho0));
                temprho1G.push_back(vector<int> (int(ceil(float(NumParticles)/float(nclasses))),InitialRho1));
                tempsigma0G.push_back(vector<int>(int(ceil(float(NumParticles)/float(nclasses))),InitialSigma0));
                tempsigma1G.push_back(vector<int>(int(ceil(float(NumParticles)/float(nclasses))),InitialSigma1));
                tempsigma2G.push_back(vector<int>(int(ceil(float(NumParticles)/float(nclasses))),InitialSigma2));
                tempweightsG.push_back(vector<long double> (int(ceil(float(NumParticles)/float(nclasses))),1/double(ceil(float(NumParticles)/float(nclasses)))));
                // class(M,B,F,lambda)
                if (LB+b==0){ noblock=true;}
                //                } //end lambda
            }//end f
            tempBG.push_back(tempG);
            temprho0BG.push_back(temprho0G);
            temprho1BG.push_back(temprho1G);
            tempsigma0BG.push_back(tempsigma0G);
            tempsigma1BG.push_back(tempsigma1G);
            tempsigma2BG.push_back(tempsigma2G);
            tempweightsBG.push_back(tempweightsG);
            b++;
        }//end b
        Alpha.push_back(tempBG);
        rho0.push_back(temprho0BG);
        rho1.push_back(temprho1BG);
        sigma0.push_back(tempsigma0BG);
        sigma1.push_back(tempsigma1BG);
        sigma2.push_back(tempsigma2BG);
        weights.push_back(tempweightsBG);
    }//end m
    //
    // initialize variables: Instance locations, Particle numbers per class
    //
    vector<vector<vector<vector<int> > > > EstimatedStatesMBG;
    vector<vector<vector<int> > > Pdist;
    for(int m=0; m<int(UM-LM+1); m++){
        int b=0;
        vector<vector<vector<int> > > tempEsBG;
        vector<vector<int> > tempPdistBG;
        while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
            noblock = false;
            vector<vector<int> > tempEsG;
            vector<int> tempPdistG;
            for(int f=0; f<int(UF-LF+1); f++){
                //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                lambda = 0;
                if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                // class(M,B,F,lambda)
                tempEsG.push_back(vector<int>(int(ceil(float(NumParticles)/float(nclasses)))));
                tempPdistG.push_back(int(ceil(float(NumParticles)/float(nclasses))));
                // class(M,B,F,lambda)
                if (LB+b==0){ noblock=true;}
                //                } //end lambda
            }//end f
            tempEsBG.push_back(tempEsG);
            tempPdistBG.push_back(tempPdistG);
            b++;
        }//end b
        EstimatedStatesMBG.push_back(tempEsBG);
        Pdist.push_back(tempPdistBG);
    }//end m
    for(int indexSequence=0; indexSequence<int(NumSequences); indexSequence++)
    {
        vector<vector<vector<vector<int> > > > TraceMBG;
        for(int m=0; m<int(UM-LM+1); m++){
            int b=0;
            vector<vector<vector<int> > > tempTraceBG;
            while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
                noblock = false;
                vector<vector<int> > tempTraceG;
                for(int f=0; f<int(UF-LF+1); f++){
                    //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                    lambda = 0;
                    if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                    // class(M,B,F,lambda)
                    vector<int> tempvec;
                    for( int p=0; p<Pdist[m][b][f]; p++){
                        tempvec.push_back(p);
                    }
                    tempTraceG.push_back(tempvec);
                    // class(M,B,F,lambda)
                    if (LB+b==0){ noblock=true;}
                    //                } //end lambda
                }//end f
                tempTraceBG.push_back(tempTraceG);
                b++;
            }//end b
            TraceMBG.push_back(tempTraceBG);
        }//end m
        int newIndesSequence = (indexSequence % NumSequences);
        cout <<"Processing sequence "<<indexSequence+1<<" of length "<<Sequences[newIndesSequence].size()<<endl;
        string currentSequence;
        currentSequence = Sequences[newIndesSequence];
        bool noblock(false);
        int lambda;
        for(int m=0; m<int(UM-LM+1); m++){
            MotifLength = LM+m;
            int Lm(int(currentSequence.size()-MotifLength+1));
            int b=0;
            while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
                noblock = false;
                for(int f=0; f<int(UF-LF+1); f++){
                    //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                    lambda = 0;
                    if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                    // class(M,B,F,lambda)
                    if(Lm<0){
                        cout <<"MOTIF<LENGTH"<<endl;
                        for (int indexParticle=0; indexParticle<Pdist[m][b][f]; indexParticle++)
                        {
                            EstimatedStatesMBG[m][b][f][indexParticle]=0;
                            weights[m][b][f][indexParticle] = weights[m][b][f][indexParticle]*1;
                            rho0[m][b][f][indexParticle]++;
                        }
                    }else{
                        vector<long double> q(Lm+1,0); // forward strands
                        for (int indexParticle=0; indexParticle<Pdist[m][b][f]; indexParticle++){
                            // Find importance distribution
                            // i=0
                            q[0] = 1;
                            for (int indexCurrentSeq = 0; indexCurrentSeq<int(currentSequence.length()); indexCurrentSeq++){
                                switch (currentSequence[indexCurrentSeq]) {
                                    case 'a':
                                    case 'A':
                                        q[0] *= NullDist[0];
                                        break;
                                    case 'c':
                                    case 'C':
                                        q[0] *= NullDist[1];
                                        break;
                                    case 'g':
                                    case 'G':
                                        q[0] *= NullDist[2];
                                        break;
                                    case 't':
                                    case 'T':
                                        q[0] *= NullDist[3];
                                        break;
                                }
                            }
                            q[0] = q[0]*double(rho0[m][b][f][indexParticle])/double(rho0[m][b][f][indexParticle]+rho1[m][b][f][indexParticle]);
                            // i>0
                            for (int indexLm=1; indexLm<Lm+1; indexLm++){
                                q[indexLm]=1;
                                for (int indexNoMotif=0; indexNoMotif<int(currentSequence.length()); indexNoMotif++){
                                    if( (indexNoMotif<indexLm-1) || (indexNoMotif>indexLm-1+MotifLength-1) ){
                                        switch(currentSequence[indexNoMotif]){
                                            case 'a':
                                            case 'A':
                                                q[indexLm] *= NullDist[0];
                                                break;
                                            case 'c':
                                            case 'C':
                                                q[indexLm] *= NullDist[1];
                                                break;
                                            case 'g':
                                            case 'G':
                                                q[indexLm] *= NullDist[2];
                                                break;
                                            case 't':
                                            case 'T':
                                                q[indexLm] *= NullDist[3];
                                                break;
                                        }
                                    }
                                }
                                //
                                // Measure block similarity
                                //
                                double palScore(0);
                                double invScore(0);
                                double dirScore(0);
                                int blocksize(0);
                                int ib(0); //block index
                                int paltri(0);
                                int invtri(0);
                                int dirtri(0);
                                vector< int> vecPalInvDir((LB+b)*3,0);
                                for (int indexMotif=0; indexMotif<=(LF+f+LB+b-1) ; indexMotif++){ //indexMotif < GapStart
                                    if (indexMotif>=LF+f){ // 1st-block
                                        // Compare 1st-block elements with those in the second block of the same length.
                                        blocksize++;
                                        switch(currentSequence[indexLm-1+indexMotif]){
                                            case 'a':
                                            case 'A':
                                                switch(currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]){
                                                    case 't':
                                                    case 'T':
                                                        palScore = palScore + 1;
                                                        vecPalInvDir[ib] = 1;
                                                        break;
                                                    case 'a':
                                                    case 'A':
                                                        invScore = invScore + 1;
                                                        vecPalInvDir[ib+LB+b] = 1;
                                                        break;
                                                }
                                                switch(currentSequence[indexLm-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]){
                                                    case 'a':
                                                    case 'A':
                                                        dirScore = dirScore + 1;
                                                        vecPalInvDir[ib+LB+b+LB+b] = 1;
                                                        break;
                                                }
                                                break;
                                            case 'c':
                                            case 'C':
                                                switch(currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]){
                                                    case 'g':
                                                    case 'G':
                                                        palScore = palScore + 1;
                                                        vecPalInvDir[ib] = 1;
                                                        break;
                                                    case 'c':
                                                    case 'C':
                                                        invScore = invScore + 1;
                                                        vecPalInvDir[ib+LB+b] = 1;
                                                        break;
                                                }
                                                switch(currentSequence[indexLm-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]){
                                                    case 'c':
                                                    case 'C':
                                                        dirScore = dirScore + 1;
                                                        vecPalInvDir[ib+LB+b+LB+b] = 1;
                                                        break;
                                                }
                                                break;
                                            case 'g':
                                            case 'G':
                                                switch(currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]){
                                                    case 'c':
                                                    case 'C':
                                                        palScore = palScore + 1;
                                                        vecPalInvDir[ib] = 1;
                                                        break;
                                                    case 'g':
                                                    case 'G':
                                                        invScore = invScore + 1;
                                                        vecPalInvDir[ib+LB+b] = 1;
                                                        break;
                                                }
                                                switch(currentSequence[indexLm-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]){
                                                    case 'g':
                                                    case 'G':
                                                        dirScore = dirScore + 1;
                                                        vecPalInvDir[ib+LB+b+LB+b] = 1;
                                                        break;
                                                }
                                                break;
                                            case 't':
                                            case 'T':
                                                switch(currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]){
                                                    case 'a':
                                                    case 'A':
                                                        palScore = palScore + 1;
                                                        vecPalInvDir[ib] = 1;
                                                        break;
                                                    case 't':
                                                    case 'T':
                                                        invScore = invScore + 1;
                                                        vecPalInvDir[ib+LB+b] = 1;
                                                        break;
                                                }
                                                switch(currentSequence[indexLm-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]){
                                                    case 't':
                                                    case 'T':
                                                        dirScore = dirScore + 1;
                                                        vecPalInvDir[ib+LB+b+LB+b] = 1;
                                                        break;
                                                }
                                                break;
                                        }// end switch current block index
                                        if (ib>=2) {
                                            if (vecPalInvDir[ib]==1) {
                                                if (vecPalInvDir[ib-1]==1) {
                                                    if (vecPalInvDir[ib-2]==1) {
                                                        paltri++;
                                                    }
                                                }
                                            }
                                            if (vecPalInvDir[ib+LB+b]==1) {
                                                if (vecPalInvDir[ib+LB+b-1]==1) {
                                                    if (vecPalInvDir[ib+LB+b-2==1]) {
                                                        invtri++;
                                                    }
                                                }
                                            }
                                            if (vecPalInvDir[ib+LB+b+LB+b]==1) {
                                                if (vecPalInvDir[ib+LB+b+LB+b-1]==1) {
                                                    if (vecPalInvDir[ib+LB+b+LB+b-2==1]) {
                                                        dirtri++;
                                                    }
                                                }
                                            }
                                        }
                                        ib++;
                                    } // 1st block
                                }// end index motif
                                float symScore(0);
                                int symbases0(0);
                                int symbases1(0);
                                int symbases2(0);
                                double hXY0(0);
                                double hXY1(0);
                                double hXY2(0);
                                double hX0(0);
                                double hX1(0);
                                double hX2(0);
                                double KLXY0(0);
                                double KLXY1(0);
                                double KLXY2(0);
                                if (SymOption==1){ //default
                                    if (b==0) {
                                        //                                    hXY0 = 1;
                                        //                                    hXY1 = 1;
                                        //                                    hXY2 = 1;
                                    }else{
                                        for (int indexMotif=0; indexMotif<=(LF+f+LB+b-1); indexMotif++){ //indexMotif < GapStart
                                            if (indexMotif>=LF+f){ // 1st-block
                                                // Cross-entropy b/w 1st-block elements with those in the 2nd-block
                                                vector<double> p1(4,0);
                                                vector<double> p2(4,0);
                                                vector<double> p3(4,0);
                                                for (int row=0; row<4; row++) {
                                                    p1[row] = (Alpha[m][b][f][indexParticle][row][indexMotif]/double(Alpha[m][b][f][indexParticle][0][indexMotif]+Alpha[m][b][f][indexParticle][1][indexMotif]+Alpha[m][b][f][indexParticle][2][indexMotif]+Alpha[m][b][f][indexParticle][3][indexMotif]));
                                                    p2[row] = (Alpha[m][b][f][indexParticle][row][(LM+m-1-lambda-indexMotif)]/double(Alpha[m][b][f][indexParticle][0][(LM+m-1-lambda-indexMotif)]+Alpha[m][b][f][indexParticle][1][(LM+m-1-lambda-indexMotif)]+Alpha[m][b][f][indexParticle][2][(LM+m-1-lambda-indexMotif)]+Alpha[m][b][f][indexParticle][3][(LM+m-1-lambda-indexMotif)]));
                                                    p3[row] = (Alpha[m][b][f][indexParticle][row][(LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]/double(Alpha[m][b][f][indexParticle][0][(LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]+Alpha[m][b][f][indexParticle][1][(LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]+Alpha[m][b][f][indexParticle][2][(LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]+Alpha[m][b][f][indexParticle][3][(LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]));
                                                }
                                                if (currentSequence[indexLm-1+indexMotif]=='a' ||
                                                    currentSequence[indexLm-1+indexMotif]=='A'){
                                                    if (currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='t' ||
                                                        currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='T') {
                                                        symbases0++;
                                                        for (int row=0; row<4; row++) {
                                                            if (p1[row]==0 || p2[row]==0 || p1[3-row]==0 || p2[3-row]==0) {
                                                                cout << "warning zero probability!\n";
                                                                continue;
                                                            }
                                                            hX0 += -(p1[row])*log(p1[row]);
                                                            hX0 += -(p2[row])*log(p2[row]);
                                                            KLXY0 += (p1[row])*log((p1[row])/(p2[3-row]));
                                                            KLXY0 += (p2[row])*log((p2[row])/(p1[3-row]));
                                                        }
                                                    }else if (currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='a' ||
                                                              currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='A') {
                                                        symbases1++;
                                                        for (int row=0; row<4; row++) {
                                                            if (p1[row]==0 || p2[row]==0) {
                                                                cout << "warning zero probability!\n";
                                                                continue;
                                                            }
                                                            hX1 += -(p1[row])*log(p1[row]);
                                                            hX1 += -(p2[row])*log(p2[row]);
                                                            KLXY1 += (p1[row])*log((p1[row])/(p2[row]));
                                                            KLXY1 += (p2[row])*log((p2[row])/(p1[row]));
                                                        }
                                                    }
                                                    if (currentSequence[indexLm-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]=='a' ||
                                                        currentSequence[indexLm-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]=='A') {
                                                        symbases2++;
                                                        for (int row=0; row<4; row++) {
                                                            if ((Alpha[m][b][f][indexParticle][row][indexMotif])==0
                                                                || (Alpha[m][b][f][indexParticle][row][(LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)])==0) {
                                                                cout << "warning zero probability!\n";
                                                                continue;
                                                            }
                                                            hX2 += -(p1[row])*log(p1[row]);
                                                            hX2 += -(p3[row])*log(p3[row]);
                                                            KLXY2 += (p1[row])*log((p1[row])/(p3[row]));
                                                            KLXY2 += (p3[row])*log((p3[row])/(p1[row]));
                                                        }
                                                    }
                                                }else if (currentSequence[indexLm-1+indexMotif]=='c' ||
                                                          currentSequence[indexLm-1+indexMotif]=='C'){
                                                    if (currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='g' ||
                                                        currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='G') {
                                                        symbases0++;
                                                        for (int row=0; row<4; row++) {
                                                            if (p1[row]==0 || p2[row]==0 || p1[3-row]==0 || p2[3-row]==0) {
                                                                cout << "warning zero probability!\n";
                                                                continue;
                                                            }
                                                            hX0 += -(p1[row])*log(p1[row]);
                                                            hX0 += -(p2[row])*log(p2[row]);
                                                            KLXY0 += (p1[row])*log((p1[row])/(p2[3-row]));
                                                            KLXY0 += (p2[row])*log((p2[row])/(p1[3-row]));
                                                        }
                                                    }else if (currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='c' ||
                                                              currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='C') {
                                                        symbases1++;
                                                        for (int row=0; row<4; row++) {
                                                            if (p1[row]==0 || p2[row]==0) {
                                                                cout << "warning zero probability!\n";
                                                                continue;
                                                            }
                                                            hX1 += -(p1[row])*log(p1[row]);
                                                            hX1 += -(p2[row])*log(p2[row]);
                                                            KLXY1 += (p1[row])*log((p1[row])/(p2[row]));
                                                            KLXY1 += (p2[row])*log((p2[row])/(p1[row]));
                                                        }
                                                    }
                                                    if (currentSequence[indexLm-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]=='c' ||
                                                        currentSequence[indexLm-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]=='C') {
                                                        symbases2++;
                                                        for (int row=0; row<4; row++) {
                                                            if ((Alpha[m][b][f][indexParticle][row][indexMotif])==0
                                                                || (Alpha[m][b][f][indexParticle][row][(LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)])==0) {
                                                                cout << "warning zero probability!\n";
                                                                continue;
                                                            }
                                                            hX2 += -(p1[row])*log(p1[row]);
                                                            hX2 += -(p3[row])*log(p3[row]);
                                                            KLXY2 += (p1[row])*log((p1[row])/(p3[row]));
                                                            KLXY2 += (p3[row])*log((p3[row])/(p1[row]));
                                                        }
                                                    }
                                                }else if (currentSequence[indexLm-1+indexMotif]=='g' ||
                                                          currentSequence[indexLm-1+indexMotif]=='G'){
                                                    if (currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='c' ||
                                                        currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='C') {
                                                        symbases0++;
                                                        for (int row=0; row<4; row++) {
                                                            if (p1[row]==0 || p2[row]==0 || p1[3-row]==0 || p2[3-row]==0) {
                                                                cout << "warning zero probability!\n";
                                                                continue;
                                                            }
                                                            hX0 += -(p1[row])*log(p1[row]);
                                                            hX0 += -(p2[row])*log(p2[row]);
                                                            KLXY0 += (p1[row])*log((p1[row])/(p2[3-row]));
                                                            KLXY0 += (p2[row])*log((p2[row])/(p1[3-row]));
                                                        }
                                                    }else if (currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='g' ||
                                                              currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='G') {
                                                        symbases1++;
                                                        for (int row=0; row<4; row++) {
                                                            if (p1[row]==0 || p2[row]==0) {
                                                                cout << "warning zero probability!\n";
                                                                continue;
                                                            }
                                                            hX1 += -(p1[row])*log(p1[row]);
                                                            hX1 += -(p2[row])*log(p2[row]);
                                                            KLXY1 += (p1[row])*log((p1[row])/(p2[row]));
                                                            KLXY1 += (p2[row])*log((p2[row])/(p1[row]));
                                                        }
                                                    }
                                                    if (currentSequence[indexLm-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]=='g' ||
                                                        currentSequence[indexLm-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]=='G') {
                                                        symbases2++;
                                                        for (int row=0; row<4; row++) {
                                                            if (p1[row]==0 || p3[row]==0) {
                                                                cout << "warning zero probability!\n";
                                                                continue;
                                                            }
                                                            hX2 += -(p1[row])*log(p1[row]);
                                                            hX2 += -(p3[row])*log(p3[row]);
                                                            KLXY2 += (p1[row])*log((p1[row])/(p3[row]));
                                                            KLXY2 += (p3[row])*log((p3[row])/(p1[row]));
                                                        }
                                                    }
                                                }else if (currentSequence[indexLm-1+indexMotif]=='t' ||
                                                          currentSequence[indexLm-1+indexMotif]=='T'){
                                                    if (currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='a' ||
                                                        currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='A') {
                                                        symbases0++;
                                                        for (int row=0; row<4; row++) {
                                                            if (p1[row]==0 || p2[row]==0 || p1[3-row]==0 || p2[3-row]==0) {
                                                                cout << "warning zero probability!\n";
                                                                continue;
                                                            }
                                                            hX0 += -(p1[row])*log(p1[row]);
                                                            hX0 += -(p2[row])*log(p2[row]);
                                                            KLXY0 += (p1[row])*log((p1[row])/(p2[3-row]));
                                                            KLXY0 += (p2[row])*log((p2[row])/(p1[3-row]));
                                                        }
                                                    }else if (currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='t' ||
                                                              currentSequence[indexLm-1+ (LM+m-1-lambda-indexMotif)]=='T') {
                                                        symbases1++;
                                                        for (int row=0; row<4; row++) {
                                                            if (p1[row]==0 || p2[row]==0) {
                                                                cout << "warning zero probability!\n";
                                                                continue;
                                                            }
                                                            hX1 += -(p1[row])*log(p1[row]);
                                                            hX1 += -(p2[row])*log(p2[row]);
                                                            KLXY1 += (p1[row])*log((p1[row])/(p2[row]));
                                                            KLXY1 += (p2[row])*log((p2[row])/(p1[row]));
                                                        }
                                                    }
                                                    if (currentSequence[indexLm-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]=='t' ||
                                                        currentSequence[indexLm-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]=='T') {
                                                        symbases2++;
                                                        for (int row=0; row<4; row++) {
                                                            if (p1[row]==0 || p3[row]==0) {
                                                                cout << "warning zero probability!\n";
                                                                continue;
                                                            }
                                                            hX2 += -(p1[row])*log(p1[row]);
                                                            hX2 += -(p3[row])*log(p3[row]);
                                                            KLXY2 += (p1[row])*log((p1[row])/(p3[row]));
                                                            KLXY2 += (p3[row])*log((p3[row])/(p1[row]));
                                                        }
                                                    }
                                                }
                                            } // 1st block
                                        }// end index motif
                                        hXY0 = w1*hX0 + w2*KLXY0;
                                        hXY1 = w1*hX1 + w2*KLXY1;
                                        hXY2 = w1*hX2 + w2*KLXY2;
                                    }
                                }
                                vector<double> symScores(3,0);
                                vector<double> symScoresnorm(3,0);
                                if(SymOption==1) { //default
                                    // Option-1. symmetrized cross-entropy
                                    if (symbases0==0) {
                                        symScores[0] = (sigma0[m][b][f][indexParticle]/double(sigma0[m][b][f][indexParticle]+sigma1[m][b][f][indexParticle]+sigma2[m][b][f][indexParticle])) * log(1+ tau);
                                    }else {
                                        symScores[0] = (sigma0[m][b][f][indexParticle]/double(sigma0[m][b][f][indexParticle]+sigma1[m][b][f][indexParticle]+sigma2[m][b][f][indexParticle]))
                                        *log(1+ tau + double(symbases0/double(1+b))*double(symbases0/double(1+hXY0)) );
                                    }
                                    if (symbases1==0){
                                        symScores[1] = (sigma1[m][b][f][indexParticle]/double(sigma0[m][b][f][indexParticle]+sigma1[m][b][f][indexParticle]+sigma2[m][b][f][indexParticle])) * log(1+ tau);
                                    }else {
                                        symScores[1] = (sigma1[m][b][f][indexParticle]/double(sigma0[m][b][f][indexParticle]+sigma1[m][b][f][indexParticle]+sigma2[m][b][f][indexParticle]))
                                        *log(1+ tau + double(symbases1/double(1+b))*double(symbases1/double(1+hXY1)) );
                                    }
                                    if (symbases2==0){
                                        symScores[2] = (sigma2[m][b][f][indexParticle]/double(sigma0[m][b][f][indexParticle]+sigma1[m][b][f][indexParticle]+sigma2[m][b][f][indexParticle])) * log(1+ tau);
                                    }else {
                                        symScores[2] = (sigma2[m][b][f][indexParticle]/double(sigma0[m][b][f][indexParticle]+sigma1[m][b][f][indexParticle]+sigma2[m][b][f][indexParticle]))
                                        *log(1+ tau + double(symbases2/double(1+b))*double(symbases2/double(1+hXY2)) );
                                    }
                                }else if(SymOption==2) {
                                    // Option-2. symmetric base(word) ratio
                                    symScores[0] = (sigma0[m][b][f][indexParticle]/double(sigma0[m][b][f][indexParticle]+sigma1[m][b][f][indexParticle]+sigma2[m][b][f][indexParticle])) *log(1+ tau + (cB*palScore+cW*double(paltri))/double(1+cB*(LB+b)+cW*max(0,LB+b-2)) ); // palScore<=LB+b, paltri<=LB+b-2
                                    symScores[1] = (sigma1[m][b][f][indexParticle]/double(sigma0[m][b][f][indexParticle]+sigma1[m][b][f][indexParticle]+sigma2[m][b][f][indexParticle])) *log(1+ tau + (cB*invScore+cW*double(invtri))/double(1+cB*(LB+b)+cW*max(0,LB+b-2)) ); // invScore<=LB+b, invtri<=LB+b-2
                                    symScores[2] = (sigma2[m][b][f][indexParticle]/double(sigma0[m][b][f][indexParticle]+sigma1[m][b][f][indexParticle]+sigma2[m][b][f][indexParticle])) *log(1+ tau + (cB*dirScore+cW*double(dirtri))/double(1+cB*(LB+b)+cW*max(0,LB+b-2)) ); // invScore<=LB+b, invtri<=LB+b-2
                                }
                                //
                                //sample from scores -- favor sites that will exclusively fit to a certain symmetry type
                                //
                                symScoresnorm[0] = symScores[0]/double(symScores[0]+symScores[1]+symScores[2]);
                                symScoresnorm[1] = symScores[1]/double(symScores[0]+symScores[1]+symScores[2]);
                                symScoresnorm[2] = symScores[2]/double(symScores[0]+symScores[1]+symScores[2]);
                                vector<float> srange(3,0); srange[0]=0; srange[1]=1; srange[2]=2;
                                int estSym = int(pmfsim(srange,symScoresnorm));
                                symScore = symScores[estSym];
                                q[indexLm]=q[indexLm]*symScore;
                                //
                                if (double(symScore)==double(0)){continue;}
                                //
                                for (int indexMotif=0; indexMotif<MotifLength; indexMotif++){
                                    if ((indexMotif>=LF+f && indexMotif<=LF+f+LB+b-1)||(indexMotif>=LM+m-1-(LF+f+lambda)-(LB+b)+1 && indexMotif<=LM+m-1-(LF+f+lambda))){
                                        // block index
                                        double normalizingAlpha(1);
                                        normalizingAlpha = Alpha[m][b][f][indexParticle][0][indexMotif]+Alpha[m][b][f][indexParticle][1][indexMotif]+Alpha[m][b][f][indexParticle][2][indexMotif]+Alpha[m][b][f][indexParticle][3][indexMotif];
                                        switch(currentSequence[indexLm-1+indexMotif]){
                                            case 'a':
                                            case 'A':
                                                q[indexLm]=q[indexLm]*double(Alpha[m][b][f][indexParticle][0][indexMotif])/normalizingAlpha;
                                                break;
                                            case 'c':
                                            case 'C':
                                                q[indexLm]=q[indexLm]*double(Alpha[m][b][f][indexParticle][1][indexMotif])/normalizingAlpha;
                                                break;
                                            case 'g':
                                            case 'G':
                                                q[indexLm]=q[indexLm]*double(Alpha[m][b][f][indexParticle][2][indexMotif])/normalizingAlpha;
                                                break;
                                            case 't':
                                            case 'T':
                                                q[indexLm]=q[indexLm]*double(Alpha[m][b][f][indexParticle][3][indexMotif])/normalizingAlpha;
                                                break;
                                        }
                                    }else {
                                        // overhang index
                                        switch(currentSequence[indexLm-1+indexMotif])
                                        {
                                            case 'a':
                                            case 'A':
                                                q[indexLm] *= NullDist[0];
                                                break;
                                            case 'c':
                                            case 'C':
                                                q[indexLm] *= NullDist[1];
                                                break;
                                            case 'g':
                                            case 'G':
                                                q[indexLm] *= NullDist[2];
                                                break;
                                            case 't':
                                            case 'T':
                                                q[indexLm] *= NullDist[3];
                                                break;
                                        }
                                    }// end if-else
                                }// end for-indexmotif
                                q[indexLm]=q[indexLm]*(double(rho1[m][b][f][indexParticle])/double(rho0[m][b][f][indexParticle]+double(rho1[m][b][f][indexParticle])))/double(Lm);
                            }// end for-indexLm
                            //
                            // Sample from the importance distribution
                            //
                            vector<float> range(q.size(),0);
                            long double normalizeQ(0);
                            for (int i=0; i<int(q.size()); i++)
                            {
                                range[i]=i;
                                normalizeQ += q[i];
                            }
                            for (int i=0; i<int(q.size()); i++)
                            {
                                q[i] = q[i]/normalizeQ;
                            }
                            int tmpES = int(pmfsimL(range, q));
                            EstimatedStatesMBG[m][b][f][indexParticle]= tmpES;
                            weights[m][b][f][indexParticle] = weights[m][b][f][indexParticle]*normalizeQ;
                            if (tmpES==int(0)){
                                rho0[m][b][f][indexParticle]++;
                            }else{
                                rho1[m][b][f][indexParticle]++;
                                if (tmpES < Lm+1){
                                    for (int indexMotif=0; indexMotif<MotifLength; indexMotif++)
                                    {
                                        switch(currentSequence[tmpES-1+indexMotif])
                                        {
                                            case 'a':
                                            case 'A':
                                                Alpha[m][b][f][indexParticle][0][indexMotif]++;
                                                // update sigma0 sigma1 sigma2
                                                if(indexMotif>=LF+f && indexMotif<=LF+f+LB+b-1){ //1st-block index
                                                    switch(currentSequence[tmpES-1+ (LM+m-1-lambda-indexMotif)]){
                                                        case 't':
                                                        case 'T':
                                                            sigma0[m][b][f][indexParticle]++;
                                                            break;
                                                        case 'a':
                                                        case 'A':
                                                            sigma1[m][b][f][indexParticle]++;
                                                            break;
                                                    }
                                                    switch(currentSequence[tmpES-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]){
                                                        case 'a':
                                                        case 'A':
                                                            sigma2[m][b][f][indexParticle]++;
                                                            break;
                                                    }
                                                }// end update sigma0 sigma1 sigma2
                                                break;
                                            case 'c':
                                            case 'C':
                                                Alpha[m][b][f][indexParticle][1][indexMotif]++;
                                                // update sigma0 sigma1 sigma2
                                                if(indexMotif>=LF+f && indexMotif<=LF+f+LB+b-1){ //1st-block index
                                                    switch(currentSequence[tmpES-1+ (LM+m-1-lambda-indexMotif)]){
                                                        case 'g':
                                                        case 'G':
                                                            sigma0[m][b][f][indexParticle]++;
                                                            break;
                                                        case 'c':
                                                        case 'C':
                                                            sigma1[m][b][f][indexParticle]++;
                                                            break;
                                                    }
                                                    switch(currentSequence[tmpES-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]){
                                                        case 'c':
                                                        case 'C':
                                                            sigma2[m][b][f][indexParticle]++;
                                                            break;
                                                    }
                                                }// end update sigma0 sigma1 sigma2
                                                break;
                                            case 'g':
                                            case 'G':
                                                Alpha[m][b][f][indexParticle][2][indexMotif]++;
                                                // update sigma0 sigma1 sigma2
                                                if(indexMotif>=LF+f && indexMotif<=LF+f+LB+b-1){ //1st-block index
                                                    switch(currentSequence[tmpES-1+ (LM+m-1-lambda-indexMotif)]){
                                                        case 'c':
                                                        case 'C':
                                                            sigma0[m][b][f][indexParticle]++;
                                                            break;
                                                        case 'g':
                                                        case 'G':
                                                            sigma1[m][b][f][indexParticle]++;
                                                            break;
                                                    }
                                                    switch(currentSequence[tmpES-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]){
                                                        case 'g':
                                                        case 'G':
                                                            sigma2[m][b][f][indexParticle]++;
                                                            break;
                                                    }
                                                }// end update sigma0 sigma1 sigma2
                                                break;
                                            case 't':
                                            case 'T':
                                                Alpha[m][b][f][indexParticle][3][indexMotif]++;
                                                // update sigma0 sigma1 sigma2
                                                if(indexMotif>=LF+f && indexMotif<=LF+f+LB+b-1){ //1st-block index
                                                    switch(currentSequence[tmpES-1+ (LM+m-1-lambda-indexMotif)]){
                                                        case 'a':
                                                        case 'A':
                                                            sigma0[m][b][f][indexParticle]++;
                                                            break;
                                                        case 't':
                                                        case 'T':
                                                            sigma1[m][b][f][indexParticle]++;
                                                            break;
                                                    }
                                                    switch(currentSequence[tmpES-1+ (LM+m-LF-f-lambda-LB-b+indexMotif-LF-f)]){
                                                        case 't':
                                                        case 'T':
                                                            sigma2[m][b][f][indexParticle]++;
                                                            break;
                                                    }
                                                }// end update sigma0 sigma1 sigma2
                                                break;
                                        }
                                    }//end indexMotif
                                }else{
                                    string currentSequence2;
                                    currentSequence2 = ReverseComplement(currentSequence.substr(tmpES-Lm-1,MotifLength));
                                    for (int indexMotif=0; indexMotif<MotifLength; indexMotif++)
                                    {
                                        switch(currentSequence2[indexMotif])
                                        {
                                            case 'a':
                                            case 'A':
                                                Alpha[m][b][f][indexParticle][0][indexMotif]++;
                                                break;
                                            case 'c':
                                            case 'C':
                                                Alpha[m][b][f][indexParticle][1][indexMotif]++;
                                                break;
                                            case 'g':
                                            case 'G':
                                                Alpha[m][b][f][indexParticle][2][indexMotif]++;
                                                break;
                                            case 't':
                                            case 'T':
                                                Alpha[m][b][f][indexParticle][3][indexMotif]++;
                                                break;
                                        }
                                        
                                    }//end indexMotif
                                } // end if-else tmpES<LM+1
                            }//end if-else tmpES==0
                        }// end for indexParticle
                    } // end if (LM < 0)
                    // class(M,B,F,lambda)
                    if (LB+b==0){ noblock=true;}
                    //                } //end lambda
                }//end of f
                b++;
            }//end of b
        }//end of m
        
        long double sumMBGweights(0);
        long double sumweight(0);
        vector<vector<vector<long double> > > MBGweights;
        for(int m=0; m<int(UM-LM+1); m++){
            vector<vector<long double> > BGweights;
            int b=0;
            while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
                noblock = false;
                vector<long double> Gweights;
                for(int f=0; f<int(UF-LF+1); f++){
                    //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                    lambda = 0;
                    if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                    // class(M,B,F,lambda)
                    long double sumweightMBG(0);
                    for(int i=0; i<Pdist[m][b][f]; i++){
                        sumweightMBG += weights[m][b][f][i];
                        sumweight += weights[m][b][f][i];
                    }
                    Gweights.push_back(sumweightMBG/Pdist[m][b][f]);
                    sumMBGweights += sumweightMBG/Pdist[m][b][f];
                    // class(M,B,F,lambda)
                    if (LB+b==0){ noblock=true;}
                    //                } //end lambda
                } // end f
                BGweights.push_back(Gweights);
                b++;
            } // end b
            MBGweights.push_back(BGweights);
        } // end m
        if(sumMBGweights==0){cerr<< "warning: zero total weight" << endl;}
        long double maxweight(weights[0][0][0][0]);
        for(int m=0; m<int(UM-LM+1); m++){
            int b=0;
            while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
                noblock = false;
                for(int f=0; f<int(UF-LF+1); f++){
                    //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                    lambda = 0;
                    if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                    // class(M,B,F,lambda)
                    MBGweights[m][b][f] = MBGweights[m][b][f]/sumMBGweights;
                    for(int i=0; i<Pdist[m][b][f]; i++){
                        weights[m][b][f][i] = weights[m][b][f][i] / sumweight;
                        // check the best class
                        if (weights[m][b][f][i]>maxweight){
                            maxweight = weights[m][b][f][i];
                        }
                    }
                    // class(M,B,F,lambda)
                    if (LB+b==0){ noblock=true;}
                    //                } //end lambda
                }
                b++;
            }
        }
        weightsPicked = weights;
        // find  Keff
        long double Keff(0);
        for(int m=0; m<int(UM-LM+1); m++){
            int b=0;
            while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
                noblock = false;
                for(int f=0; f<int(UF-LF+1); f++){
                    //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                    lambda = 0;
                    if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                    // class(M,B,F,lambda)
                    for(int i=0; i<Pdist[m][b][f]; i++){
                        Keff += (weights[m][b][f][i]*weights[m][b][f][i]);
                    }
                    // class(M,B,F,lambda)
                    if (LB+b==0){ noblock=true;}
                    //                } //end lambda
                }
                b++;
            }
        }
        Keff = 1/Keff;
        int Pthr(ceil( (float(NumParticles)/float(nclasses)) ));
        if(float(Keff)<float(NumParticles)/float(nclasses)){
            cout << "Resampling...";
            // Choose the number of particles for each class, ￼newPdist[m][b][f], according to a multinomial distribution with parameters MBGweights
            vector<double> vecP;
            vector<double> vecPdist;
            vector<float> Prange;
            int pr(0);
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
                        vecP.push_back(Pdist[m][b][f]);
                        vecPdist.push_back(MBGweights[m][b][f]); // posterior class prob
                        Prange.push_back(pr);
                        pr++;
                        // class(M,B,F,lambda)
                        if (LB+b==0){ noblock=true;}
                        //                } //end lambda
                    }
                    b++;
                }
            }
            int sumPdist(0);
            vector<vector<vector<int> > > newPdist;
            for(int m=0; m<int(UM-LM+1); m++){
                vector<vector<int> > newPdistBG;
                int b=0;
                while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
                    noblock = false;
                    vector<int> newPdistG;
                    for(int f=0; f<int(UF-LF+1); f++){
                        //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                        lambda = 0;
                        if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                        // class(M,B,F,lambda)
                        int newsample = vecP[int(pmfsim(Prange, vecPdist))];
                        if (newsample<Pthr){
                            newPdistG.push_back(Pthr);
                            sumPdist += Pthr;
                        }else{
                            newPdistG.push_back(newsample);
                            sumPdist += newsample;
                        }
                        // class(M,B,F,lambda)
                        if (LB+b==0){ noblock=true;}
                        //                } //end lambda
                    }
                    newPdistBG.push_back(newPdistG);
                    b++;
                }
                newPdist.push_back(newPdistBG);
            }
            // Reduce the number of particles until sumPdist = NumParticles, from the class with the most particles (the random class among the equally-most-particled ones)
            while(sumPdist>NumParticles){
                int maxPdist(newPdist[0][0][0]); int pm(0); int pb(0); int pf(0);
                vector<vector<int> > equalInds;
                vector<int> rangeEq; int eq(0);
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
                            if(newPdist[m][b][f]>maxPdist){
                                maxPdist = newPdist[m][b][f];
                                pm = m; pb = b; pf = f;
                            } else if(newPdist[m][b][f]>=maxPdist){
                                vector<int> maxInds(3,0);
                                maxInds[0] = m;
                                maxInds[1] = b;
                                maxInds[2] = f;
                                equalInds.push_back(maxInds);
                                rangeEq.push_back(eq);
                                eq++;
                            }
                            // class(M,B,F,lambda)
                            if (LB+b==0){ noblock=true;}
                            //                } //end lambda
                        }
                        b++;
                    }
                }
                if (eq>0){
                    int maxInd = rand() % eq;
                    pm = equalInds[maxInd][0];
                    pb = equalInds[maxInd][1];
                    pf = equalInds[maxInd][2];
                }
                int tempVal = newPdist[pm][pb][pf]-1;
                newPdist[pm][pb][pf] = max(1, tempVal);
                sumPdist = 0;
                for(int m=0; m<int(UM-LM+1); m++){
                    int b=0;
                    while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
                        noblock = false;
                        for(int f=0; f<int(UF-LF+1); f++){
                            //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                            lambda = 0;
                            if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                            // class(M,B,F,lambda)
                            sumPdist += newPdist[m][b][f];
                            // class(M,B,F,lambda)
                            if (LB+b==0){ noblock=true;}
                            //                } //end lambda
                        }
                        b++;
                    }
                }
            }// end while(sumPdist>NumParticles){
            // Sample "newPdist" new particles from the set of previous particles (and update rho0/1,sigma0/1/2,weight,Alpha,location) of the class, with probabilities proportional to their "weights". Assign equal weight to each new sample within a class.
            vector<vector<vector<vector<int> > > > resRho0;
            vector<vector<vector<vector<int> > > > resRho1;
            vector<vector<vector<vector<int> > > > resSigma0;
            vector<vector<vector<vector<int> > > > resSigma1;
            vector<vector<vector<vector<int> > > > resSigma2;
            vector<vector<vector<vector<vector<vector<double> > > > > > resAlpha;
            vector<vector<vector<vector<int> > > > resEstimatedStatesMBG;
            vector<vector<vector<vector<int> > > > resTraceMBG;
            vector<vector<vector<vector<long double> > > > resWeights;
            vector<vector<vector<vector<long double> > > > picWeights;
            long double sumPicWeights(0);
            for(int m=0; m<int(UM-LM+1); m++){
                vector<vector<vector<int> > > resrho0BG;
                vector<vector<vector<int> > > resrho1BG;
                vector<vector<vector<int> > > ressigma0BG;
                vector<vector<vector<int> > > ressigma1BG;
                vector<vector<vector<int> > > ressigma2BG;
                vector<vector<vector<vector<vector<double> > > > > resAlphaBG;
                vector<vector<vector<int> > > resESBG;
                vector<vector<vector<int> > > resTRBG;
                vector<vector<vector<long double> > > resWeightsBG;
                vector<vector<vector<long double> > > picWeightsBG;
                int b=0;
                while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
                    noblock = false;
                    vector<vector<int> > resrho0G;
                    vector<vector<int> > resrho1G;
                    vector<vector<int> > ressigma0G;
                    vector<vector<int> > ressigma1G;
                    vector<vector<int> > ressigma2G;
                    vector<vector<vector<vector<double> > > > resAlphaG;
                    vector<vector<int> > resESG;
                    vector<vector<int> > resTRG;
                    vector<vector<long double> > resWeightsG;
                    vector<vector<long double> > picWeightsG;
                    for(int f=0; f<int(UF-LF+1); f++){
                        //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                        lambda = 0;
                        if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                        // class(M,B,F,lambda)
                        long double sumClassWeights(0);
                        vector<float> range(Pdist[m][b][f],0);
                        for(int i=0; i<Pdist[m][b][f]; i++){
                            range[i]=i;
                            sumClassWeights += weights[m][b][f][i];
                        }
                        vector <long double> normalizedClassWeights;
                        for(int i=0; i<Pdist[m][b][f]; i++){
                            normalizedClassWeights.push_back(weights[m][b][f][i]/sumClassWeights);
                        }
                        int tmpP = int(newPdist[m][b][f]);
                        vector <int> sample(tmpP,0);
                        vector<int>resvec0;
                        vector<int>resvec1;
                        vector<int>ressvec0;
                        vector<int>ressvec1;
                        vector<int>ressvec2;
                        vector<vector<vector<double> > > resPwmAlphaSample;
                        vector<int> resES;
                        vector<int> resTR;
                        vector<long double> picWeightsvec;
                        vector<long double> resWeightsvec;
                        for (int i=0; i<tmpP; i++){
                            // Sample
                            sample[i] = int(pmfsimL(range, normalizedClassWeights));
                            resvec0.push_back(int(rho0[m][b][f][sample[i]]));
                            resvec1.push_back(int(rho1[m][b][f][sample[i]]));
                            ressvec0.push_back(int(sigma0[m][b][f][sample[i]]));
                            ressvec1.push_back(int(sigma1[m][b][f][sample[i]]));
                            ressvec2.push_back(int(sigma2[m][b][f][sample[i]]));
                            picWeightsvec.push_back(weights[m][b][f][sample[i]]);
                            sumPicWeights += weights[m][b][f][sample[i]];
                            vector<vector<double> > resPwmAlphaM;
                            for (int j=0; j<4; j++)
                            {
                                vector<double> resPwmAlpha;
                                for (int mm=0; mm<(LM+m); mm++){
                                    resPwmAlpha.push_back(Alpha[m][b][f][sample[i]][j][mm]);
                                }
                                resPwmAlphaM.push_back(resPwmAlpha);
                            }
                            resPwmAlphaSample.push_back(resPwmAlphaM);
                            resES.push_back(EstimatedStatesMBG[m][b][f][sample[i]]);
                            resTR.push_back(sample[i]); //save the index of the particle that is sampled
                            // Assign equal weight to each new sample within a class
                            resWeightsvec.push_back(sumClassWeights/double(tmpP));
                        } //end particles
                        resrho0G.push_back(resvec0);
                        resrho1G.push_back(resvec1);
                        ressigma0G.push_back(ressvec0);
                        ressigma1G.push_back(ressvec1);
                        ressigma2G.push_back(ressvec2);
                        resAlphaG.push_back(resPwmAlphaSample);
                        resESG.push_back(resES);
                        resTRG.push_back(resTR);
                        resWeightsG.push_back(resWeightsvec);
                        picWeightsG.push_back(picWeightsvec);
                        // class(M,B,F,lambda)
                        if (LB+b==0){ noblock=true;}
                        //                } //end lambda
                    } //end f
                    resrho0BG.push_back(resrho0G);
                    resrho1BG.push_back(resrho1G);
                    ressigma0BG.push_back(ressigma0G);
                    ressigma1BG.push_back(ressigma1G);
                    ressigma2BG.push_back(ressigma2G);
                    resAlphaBG.push_back(resAlphaG);
                    resESBG.push_back(resESG);
                    resTRBG.push_back(resTRG);
                    resWeightsBG.push_back(resWeightsG);
                    picWeightsBG.push_back(picWeightsG);
                    b++;
                } //end b
                resRho0.push_back(resrho0BG);
                resRho1.push_back(resrho1BG);
                resSigma0.push_back(ressigma0BG);
                resSigma1.push_back(ressigma1BG);
                resSigma2.push_back(ressigma2BG);
                resAlpha.push_back(resAlphaBG);
                resEstimatedStatesMBG.push_back(resESBG);
                resTraceMBG.push_back(resTRBG);
                resWeights.push_back(resWeightsBG);
                picWeights.push_back(picWeightsBG);
            } //end m
            cout << "...";
            rho0 = resRho0;
            rho1 = resRho1;
            sigma0 = resSigma0;
            sigma1 = resSigma1;
            sigma2 = resSigma2;
            Alpha = resAlpha;
            EstimatedStatesMBG = resEstimatedStatesMBG;
            TraceMBG = resTraceMBG;
            weights = resWeights;
            weightsPicked = picWeights;
            Pdist = newPdist;
            for(int m=0; m<int(UM-LM+1); m++){
                int b=0;
                while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
                    noblock = false;
                    for(int f=0; f<int(UF-LF+1); f++){
                        //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                        lambda = 0;
                        if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                        // class(M,B,F,lambda)
                        for (int iP=0; iP<Pdist[m][b][f]; iP++ ) {
                            weightsPicked[m][b][f][iP] = weightsPicked[m][b][f][iP]/sumPicWeights;
                        }
                        // class(M,B,F,lambda)
                        if (LB+b==0){ noblock=true;}
                        //                } //end lambda
                    }
                    b++;
                }
            }
            cout << " done.\n";
        } // end if needed to resample
        //
        // Get the best class: Class-based
        //
        long double maxMBGweight = MBGweights[0][0][0];
        SetEstiM(0); SetEstiB(0); SetEstiF(0);
        for(int m=0; m<int(UM-LM+1); m++){
            int b=0;
            while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
                noblock = false;
                for(int f=0; f<int(UF-LF+1); f++){
                    //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                    lambda = 0;
                    if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                    // class(M,B,F,lambda)
                    if (MBGweights[m][b][f]>maxMBGweight){
                        maxMBGweight = MBGweights[m][b][f];
                        SetEstiM(m); SetEstiB(b); SetEstiF(f);
                    }
                    // class(M,B,F,lambda)
                    if (LB+b==0){ noblock=true;}
                    //                } //end lambda
                }
                b++;
            }
        }
        long double bestweight = weightsPicked[estiM][estiB][estiF][0];
        int bestParticle(0);
        for (int i=0; i<int(weightsPicked[estiM][estiB][estiF].size()); i++){
            if(weightsPicked[estiM][estiB][estiF][i]>bestweight){
                bestweight = weightsPicked[estiM][estiB][estiF][i];
                bestParticle = i;
            }
        }
        SetEstiP(bestParticle);
        // Save sufficient statistics: ALPHA, sigma0est, sigma1est, sigma2est
        for (int ir=0; ir<int(4); ir++) {
            for (int ic=0; ic<int(LM+estiM); ic++) {
                // overhangs and blocks
                ALPHA[estiM][estiB][estiF][ir][ic] = Alpha[estiM][estiB][estiF][estiP][ir][ic]-initcount; //subtract the initial values
            }
        }
        sigma0est = sigma0[estiM][estiB][estiF][estiP];
        sigma1est = sigma1[estiM][estiB][estiF][estiP];
        sigma2est = sigma2[estiM][estiB][estiF][estiP];
        EstimatedStatesSeq.push_back(EstimatedStatesMBG);
        Trace.push_back(TraceMBG); // if resampling occured new set of particle indices are saved, otherwise the same set of particle indices are passed to next sequence.
        cout << "State: " << EstimatedStatesMBG[estiM][estiB][estiF][estiP] << " @ ";
        cout << "MBFP: " << estiM << " " << estiB << " " << estiF << " " << estiP << endl;
    }// end indexSequence
    fileOut << "Sequence No\tMotif Instance Start\tMotif Instance Sequence"<<endl;
    fileOut << "----------------------------------------\n";
    //
    // Trace-back the best particle's states (through resampling events)
    //
    vector<int> ES;
    int iT = estiP;
    int tmpiT;
    for(int indexSequence=NumSequences-1; indexSequence>=int(0); indexSequence--){
        ES.push_back(EstimatedStatesSeq[indexSequence][estiM][estiB][estiF][iT]);
        tmpiT = iT;
        iT = Trace[indexSequence][estiM][estiB][estiF][tmpiT];
    }
    //
    // Write Motif Instances
    //
    cout << "Writing motif instances...\n";
    for(int indexSequence=0; indexSequence<int(NumSequences); indexSequence++){
        string currentSequence2;
        currentSequence2 = Sequences[indexSequence];
        int Lm2 = int(currentSequence2.size()-(LM+estiM)+1);
        fileOut <<indexSequence+1 <<"\t";
        int tempS = ES[NumSequences-1-indexSequence];
        if (tempS < Lm2+1 || Lm2<0){
            fileOut << tempS << "\t";
            InferredPositions1.push_back(tempS);
        }else{
            fileOut << tempS-Lm2-1 << "\t";
            InferredPositions1.push_back(tempS-Lm2-1);
        }
        if(tempS>0){
            fileOut << currentSequence2.substr(tempS-1,(LM+estiM)) << endl;
        }else{
            fileOut << "No instance" << endl;
        }
    }
    fileOut << "----------------------------------------\n";
    fileOut << "Motif Structure\n";
    fileOut << "----------------------------------------\n";
    fileOut << "Motif width: " <<  (LM+estiM) << endl;
    fileOut << "Motif block width: " << (LB+estiB) << endl;
    fileOut << "Motif offset width: " << (LF+estiF) << endl;
    cout << "Writing final weights...";
    long double finalsum(0);
    for(int m=0; m<int(UM-LM+1); m++){
        int b=0;
        while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
            noblock = false;
            for(int f=0; f<int(UF-LF+1); f++){
                //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                lambda = 0;
                if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                // class(M,B,F,lambda)
                long double tmpsum(0);
                for(int iP=0; iP<int(weightsPicked[m][b][f].size()); iP++){
                    tmpsum += weightsPicked[m][b][f][iP];
                }
                finalsum += tmpsum/double(weightsPicked[m][b][f].size());
                // class(M,B,F,lambda)
                if (LB+b==0){ noblock=true;}
                //                } //end lambda
            }
            b++;
        }
    }
    vector<vector<vector<long double> > > finalWeights;
    for(int m=0; m<int(UM-LM+1); m++){
        vector<vector<long double> > tempFWBG;
        int b=0;
        while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
            noblock = false;
            vector<long double> tempFWG;
            for(int f=0; f<int(UF-LF+1); f++){
                //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
                lambda = 0;
                if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
                // class(M,B,F,lambda)
                long double tmpsum(0);
                for(int iP=0; iP<int(weightsPicked[m][b][f].size()); iP++){
                    tmpsum += weightsPicked[m][b][f][iP];
                }
                tempFWG.push_back((tmpsum/double(weightsPicked[m][b][f].size()))/finalsum);
                // class(M,B,F,lambda)
                if (LB+b==0){ noblock=true;}
                //                } //end lambda
            }
            tempFWBG.push_back(tempFWG);
            b++;
        }
        finalWeights.push_back(tempFWBG);
    }
    cout << "done.\n";
    /*
     fileOut << "Distribution of offset width for each motif width (M) and block width (B): \n";
     for(int m=0; m<int(UM-LM+1); m++){
     int b=0;
     while(b<int(UB-LB+1) && (LB+b) <= floor((LM+m)/2)){
     noblock = false;
     fileOut << "M: " << (LM+m) << ", B: " << (LB+b) << endl;
     for(int f=0; f<int(UF-LF+1); f++){
     //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
     lambda = 0;
     if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
     // class(M,B,F,lambda)
     fileOut << (LF+f) << "\t\t";
     // class(M,B,F,lambda)
     if (LB+b==0){ noblock=true;}
     //                } //end lambda
     }
     fileOut << "\n";
     noblock = false;
     for(int f=0; f<int(UF-LF+1); f++){
     //                for (lambda = -(LF+f); lambda <= LM+m-2*(LB+b)-2*(LF+f); lambda++){
     lambda = 0;
     if (LF+f > LM+m-2*(LB+b)-(LF+f+lambda) || noblock){ continue;}
     // class(M,B,F,lambda)
     fileOut << finalWeights[m][b][f] << "\t";
     // class(M,B,F,lambda)
     if (LB+b==0){ noblock=true;}
     //                } //end lambda
     }
     fileOut << "\n";
     b++;
     }
     }
     */
    double information_content(0);
    vector<double> temp(4,0.25);
    vector<vector<double> > tPWM((LM+estiM),temp);
    PWM = tPWM;
    fileOut << "----------------------------------------\n";
    fileOut << "PWM" << endl;
    fileOut << "----------------------------------------\n";
    fileOut << "A\tC\tG\tT\n";
    cout << "----------------------------------------\n";
    cout << "PWM" << endl;
    cout << "----------------------------------------\n";
    cout << "A\tC\tG\tT\n";
    double total_variance(0);
    double total_entropy_variance(0);
    for (int column=0; column<(LM+estiM); column++)
    {
        double Normaling(0);
        // overhangs and blocks
        for (int row = 0; row<4; row++){
            Normaling+=ALPHA[estiM][estiB][estiF][row][column];
        }
        double tmpvarmean(0);
        double tmpentvarmean(0);
        for (int row = 0; row<4; row++)
        {
            // overhangs and blocks
            PWM[column][row]=ALPHA[estiM][estiB][estiF][row][column]/Normaling;
            if (PWM[column][row]==0) {
            }else{
                information_content += PWM[column][row] * log(PWM[column][row]/NullDist[row]);
                tmpentvarmean += PWM[column][row] * log(PWM[column][row]/NullDist[row]);
            }
            tmpvarmean += PWM[column][row];
            fileOut << PWM[column][row] << "\t";
            cout << PWM[column][row];
            if (row+1<4) {
                cout << "\t";
            }
        }
        tmpentvarmean = tmpentvarmean/double(4);
        tmpvarmean = tmpvarmean/double(4);
        double tmpvar(0);
        double tmpentvar(0);
        for (int row = 0; row<4; row++){
            tmpvar += (PWM[column][row]-tmpvarmean)*(PWM[column][row]-tmpvarmean);
            if (PWM[column][row]==0) {
                tmpentvar += (0 - tmpentvarmean) * (0 - tmpentvarmean);
            }else{
                tmpentvar += (PWM[column][row] * log(PWM[column][row]/NullDist[row]) - tmpentvarmean)*(PWM[column][row] * log(PWM[column][row]/NullDist[row]) -tmpentvarmean);
            }
        }
        total_variance += tmpvar/double(3);
        total_entropy_variance += tmpentvar/double(3);
        fileOut << endl;
        cout << endl;
    }// end column
    total_variance = total_variance/double((LM+estiM));
    total_entropy_variance = total_entropy_variance/double((LM+estiM));
    
    fileOut << "----------------------------------------\n";
    cout << "----------------------------------------\n";
    fileOut << "Information content: " << information_content << endl;
    fileOut << "Total variance: " << total_variance << endl;
    fileOut << "Total entropy variance: " << total_entropy_variance << endl;
    cout << "Information content: " << information_content << endl;
    cout << "Total variance: " << total_variance << endl;
    cout << "Total entropy variance: " << total_entropy_variance << endl;
    fileOut << "Symmetry type: " << round(100*sigma0est/float(sigma0est+sigma1est+sigma2est)) << "% palindromic,\t";
    fileOut << round(100*sigma1est/float(sigma0est+sigma1est+sigma2est)) << "% inverted-repeat,\t";
    fileOut << round(100*sigma2est/float(sigma0est+sigma1est+sigma2est)) << "% direct-repeat.\n";
    cout << "Symmetry type: " << round(100*sigma0est/float(sigma0est+sigma1est+sigma2est)) << "% palindromic,\t";
    cout << round(100*sigma1est/float(sigma0est+sigma1est+sigma2est)) << "% inverted-repeat,\t";
    cout << round(100*sigma2est/float(sigma0est+sigma1est+sigma2est)) << "% direct-repeat.\n";
    cout << "----------------------------------------\n";
    fileOut << "----------------------------------------\n";
}// end ComputeParticleFilter()

float Database::pmfsim(const vector<float>& range, const vector<double>& Px)
{
    vector<float> CDF(Px.size(),0);
    CDF[0]=float(Px[0]);
    for (int i=1; i<int(Px.size()); i++)
    {
        CDF[i]=CDF[i-1]+float(Px[i]);
    }
    float randomnumber;
    randomnumber = float(rand() % 10000)/10000;
    bool flag;
    flag = 1;
    int i(0);
    while(flag)
    {
        if(CDF[i]>randomnumber)
        {
            flag = 0;
        }
        else
        {
            i++;
        }
        if(i==int(Px.size()))
        {
            flag = 0;
            i--;
        }
    }
    return range[i];
}

float Database::pmfsimL(const vector<float>& range, const vector<long double>& Px)
{
    vector<float> CDF(Px.size(),0);
    CDF[0]=float(Px[0]);
    for (int i=1; i<int(Px.size()); i++)
    {
        CDF[i]=CDF[i-1]+float(Px[i]);
    }
    float randomnumber;
    randomnumber = float(rand() % 10000)/10000;
    bool flag;
    flag = 1;
    int i(0);
    while(flag)
    {
        if(CDF[i]>randomnumber)
        {
            flag = 0;
        }
        else
        {
            i++;
        }
        if(i==int(Px.size()))
        {
            flag = 0;
            i--;
        }
    }
    return range[i];
}

long double Database::ComputePosterior1Motif (const string &CurrentSequence1, const int pos1, const vector<vector<float> > &Alpha, const vector<vector<float> > &Beta, const int MotifLength)
{
    long double Posterior = pow(0.25,int(MotifLength));
    vector<int> r(4);
    float gamma1(1),gamma2(1),gamma3(1);
    int posMotif;
    for (posMotif=0; posMotif<MotifLength; posMotif++)
    {
        r[0]=0;
        r[1]=0;
        r[2]=0;
        r[3]=0;
        switch (CurrentSequence1[pos1+posMotif])
        {
            case 'A':
                r[0]=1;
                break;
            case 'a':
                r[0]=1;
                break;
            case 'C':
                r[1]=1;
                break;
            case 'c':
                r[1]=1;
                break;
            case 'G':
                r[2]=1;
                break;
            case 'g':
                r[2]=1;
                break;
            case 'T':
                r[3]=1;
                break;
            case 't':
                r[3]=1;
                break;
            default:
                cerr << "Problem with the sequences' format in Sequence 1.1"<<endl;
        }
        switch (r[0])
        {
            case 0:
                gamma1 = 1;
                break;
            case 1:
                gamma1 = Alpha[0][posMotif];
                break;
        }
        gamma2 = Alpha[0][posMotif]+Beta[0][posMotif];
        switch (r[1]+r[2]+r[3])
        {
            case 0:
                gamma3 = 1;
                break;
            case 1:
                gamma3 = Beta[0][posMotif];
                break;
        }
        Posterior=Posterior*gamma1/gamma2*gamma3;
        // For i = 1
        switch (r[1])
        {
            case 0:
                gamma1 = 1;
                break;
            case 1:
                gamma1 = Alpha[1][posMotif];
                break;
        }
        switch (r[1]+r[2]+r[3])
        {
            case 0:
                gamma2 = 1;
                break;
            case 1:
                gamma2 = Alpha[1][posMotif]+Beta[1][posMotif];
                break;
        }
        switch (r[2]+r[3])
        {
            case 0:
                gamma3 = 1;
                break;
            case 1:
                gamma3 = Beta[1][posMotif];
                break;
        }
        Posterior=Posterior*gamma1/gamma2*gamma3;
        // For i = 2
        switch (r[2])
        {
            case 0:
                gamma1 = 1;
                break;
            case 1:
                gamma1 = Alpha[2][posMotif];
                break;
        }
        switch (r[2]+r[3])
        {
            case 0:
                gamma2 = 1;
                break;
            case 1:
                gamma2 = Alpha[2][posMotif]+Beta[2][posMotif];
                break;
        }
        switch (r[3])
        {
            case 0:
                gamma3 = 1;
                break;
            case 1:
                gamma3 = Beta[2][posMotif];
                break;
        }
        Posterior=Posterior*gamma1/gamma2*gamma3;
    }
    return Posterior;
}


long double Database::ComputePosterior2Motifs (const string &CurrentSequence1, const string &CurrentSequence2, const vector<vector<float> > &Alpha, const vector<vector<float> > &Beta, const int MotifLength)
{
    
    int pos1 = 0;
    int pos2 = 0;
    long double Posterior(1);
    vector<int> r(4);
    float gamma1(1),gamma2(1),gamma3(1);
    int posMotif;
    for (posMotif=0; posMotif<MotifLength; posMotif++)
    {
        r[0]=0;
        r[1]=0;
        r[2]=0;
        r[3]=0;
        switch (CurrentSequence1[pos1+posMotif])
        {
            case 'A':
                r[0]=1;
                break;
            case 'a':
                r[0]=1;
                break;
            case 'C':
                r[1]=1;
                break;
            case 'c':
                r[1]=1;
                break;
            case 'G':
                r[2]=1;
                break;
            case 'g':
                r[2]=1;
                break;
            case 'T':
                r[3]=1;
                break;
            case 't':
                r[3]=1;
                break;
            default:
                cerr << "Problem with the sequences' format in Sequence 1.2"<<endl;
        }
        switch (CurrentSequence2[pos2+posMotif])
        {
            case 'A':
            case 'a':
                r[0]=r[0]+1;
                break;
            case 'C':
            case 'c':
                r[1]=r[1]+1;
                break;
            case 'G':
            case 'g':
                r[2]=r[2]+1;
                break;
            case 'T':
            case 't':
                r[3]=r[3]+1;
                break;
            default:
                cerr << "Problem with the sequences' format in Sequence 2.1"<<endl;
        }
        switch (r[0])
        {
            case 0:
                gamma1 = 1;
                break;
            case 1:
                gamma1 = Alpha[0][posMotif];
                break;
            case 2:
                gamma1 = (Alpha[0][posMotif])*(Alpha[0][posMotif]+1);
                break;
        }
        gamma2 = (Alpha[0][posMotif]+Beta[0][posMotif])*(Alpha[0][posMotif]+Beta[0][posMotif]+1);
        switch (r[1]+r[2]+r[3])
        {
            case 0:
                gamma3 = 1;
                break;
            case 1:
                gamma3 = Beta[0][posMotif];
                break;
            case 2:
                gamma3 = (Beta[0][posMotif])*(Beta[0][posMotif]+1);
                break;
        }
        Posterior=Posterior*gamma1/gamma2*gamma3;
        // For i = 1
        switch (r[1])
        {
            case 0:
                gamma1 = 1;
                break;
            case 1:
                gamma1 = Alpha[1][posMotif];
                break;
            case 2:
                gamma1 = (Alpha[1][posMotif])*(Alpha[1][posMotif]+1);
                break;
        }
        switch (r[1]+r[2]+r[3])
        {
            case 0:
                gamma2 = 1;
                break;
            case 1:
                gamma2 = Alpha[1][posMotif]+Beta[1][posMotif];
                break;
            case 2:
                gamma2 =(Alpha[1][posMotif]+Beta[1][posMotif])*(Alpha[1][posMotif]+Beta[1][posMotif]+1);
                break;
        }
        switch (r[2]+r[3])
        {
            case 0:
                gamma3 = 1;
                break;
            case 1:
                gamma3 = Beta[1][posMotif];
                break;
            case 2:
                gamma3 = (Beta[1][posMotif])*(Beta[1][posMotif]+1);
                break;
        }
        Posterior=Posterior*gamma1/gamma2*gamma3;
        // For i = 2
        switch (r[2])
        {
            case 0:
                gamma1 = 1;
                break;
            case 1:
                gamma1 = Alpha[2][posMotif];
                break;
            case 2:
                gamma1 = (Alpha[2][posMotif])*(Alpha[2][posMotif]+1);
                break;
        }
        switch (r[2]+r[3])
        {
            case 0:
                gamma2 = 1;
                break;
            case 1:
                gamma2 = Alpha[2][posMotif]+Beta[2][posMotif];
                break;
            case 2:
                gamma2 =(Alpha[2][posMotif]+Beta[2][posMotif])*(Alpha[2][posMotif]+Beta[2][posMotif]+1);
                break;
        }
        switch (r[3])
        {
            case 0:
                gamma3 = 1;
                break;
            case 1:
                gamma3 = Beta[2][posMotif];
                break;
            case 2:
                gamma3 = (Beta[2][posMotif])*(Beta[2][posMotif]+1);
                break;
        }
        Posterior=Posterior*gamma1/gamma2*gamma3;
    }
    
    return Posterior;
}


long double Database::FindML (const string &CurrentSequence, const vector<vector<float> > &Alpha, const vector<vector<float> > &Beta, const int MotifLength, vector<int> &BestPositioning)
{
    vector<int> tempxy2(2,0);
    long double Posterior;
    long double maxML(0);
    string RCsequence;
    RCsequence = ReverseComplement(CurrentSequence);
    // 1 motif
    for (int pos=0; pos< int(CurrentSequence.size()-MotifLength+1); pos++)
    {
        Posterior = ComputePosterior1Motif(CurrentSequence, pos, Alpha, Beta, MotifLength);
        if (Posterior>maxML)
        {
            tempxy2[0] = pos+1;
            tempxy2[1] = 0;
            BestPositioning=tempxy2;
            maxML = Posterior;
        }
    }
    for (int pos=0; pos< int(CurrentSequence.size()-MotifLength+1); pos++)
    {
        string tempSequence;
        tempSequence = ReverseComplement(CurrentSequence.substr(pos, MotifLength));
        Posterior = 0*ComputePosterior1Motif(tempSequence, 0, Alpha, Beta, MotifLength);
        if (Posterior>maxML)
        {
            tempxy2[0] = pos+1;
            tempxy2[1] = 0;
            BestPositioning=tempxy2;
            maxML = Posterior;
        }
    }
    // 2 motifs
    for (int pos1=0; pos1< int(CurrentSequence.size()-MotifLength+1); pos1++)
    {
        for (int pos2=pos1+MotifLength; pos2< int(CurrentSequence.size()-MotifLength+1); pos2++)
        {
            Posterior = ComputePosterior2Motifs(CurrentSequence.substr(pos1,MotifLength), CurrentSequence.substr(pos2,MotifLength), Alpha, Beta, MotifLength);
            if (Posterior>maxML)
            {
                tempxy2[0] = pos1+1;
                tempxy2[1] = pos2+1;
                BestPositioning=tempxy2;
                maxML = Posterior;
            }
        }
    }
    for (int pos1=0; pos1< int(CurrentSequence.size()-MotifLength+1); pos1++)
    {
        for (int pos2=0; pos2< int(CurrentSequence.size()-MotifLength+1); pos2++)
        {
            if(pos2>pos1+MotifLength || pos2+MotifLength < pos1)
            {
                Posterior = 0*ComputePosterior2Motifs(CurrentSequence.substr(pos1,MotifLength), ReverseComplement(CurrentSequence.substr(pos2,MotifLength)), Alpha, Beta, MotifLength);
                if (Posterior>maxML)
                {
                    tempxy2[0] = pos1+1;
                    tempxy2[1] = pos2+1;
                    BestPositioning=tempxy2;
                    maxML = Posterior;
                }
            }
        }
    }
    for (int pos1=0; pos1< int(CurrentSequence.size()-MotifLength+1); pos1++)
    {
        for (int pos2=pos1+MotifLength; pos2< int(CurrentSequence.size()-MotifLength+1); pos2++)
        {
            Posterior = 0*ComputePosterior2Motifs(ReverseComplement(CurrentSequence.substr(pos1,MotifLength)), ReverseComplement(CurrentSequence.substr(pos2,MotifLength)), Alpha, Beta, MotifLength);
            if (Posterior>maxML)
            {
                tempxy2[0] = pos1+1;
                tempxy2[1] = pos2+1;
                BestPositioning=tempxy2;
                maxML = Posterior;
            }
        }
    }
    return maxML;
}

string Database::ReverseComplement(const string Sequence1)
{
    string Sequence2;
    for (int i= int(Sequence1.size()-1); i>-1; i--)
    {
        switch (Sequence1[i])
        {
            case 'A':
            case 'a':
                Sequence2.push_back('T');
                break;
            case 'C':
            case 'c':
                Sequence2.push_back('G');
                break;
            case 'G':
            case 'g':
                Sequence2.push_back('C');
                break;
            case 'T':
            case 't':
                Sequence2.push_back('A');
                break;
            default:
                cerr << "Problem with the sequences' format in ReverseComplement"<<endl;
        }
        
    }
    return Sequence2;
}

void Database::EliminateInstances()
{
    vector< int > NumberInstances(NumSequences,0);
    for (int i=0; i<int(NumSequences); i++)
    {
        if(IndicatorAtLeastTwoMotif[i])
        {
            NumberInstances[i] = 2;
        }
        else
        {
            if(IndicatorAtLeastOneMotif[i])
            {
                NumberInstances[i] = 1;
            }
            else
            {
                NumberInstances[i] = 0;
            }
        }
        if(IndicatorAtLeastOneMotif[i]==0)
        {
            Sequences2.push_back(Sequences[i]);
        }
        else
        {
            if(IndicatorAtLeastOneMotif[i]==1)
            {
                Sequences2.push_back(Sequences[i].substr(0,InferredPositions[i][0]));
                Sequences2.push_back(Sequences[i].substr(InferredPositions[i][0]+(LM+estiM), Sequences[i].size()));
            }
            else
            {
                if(InferredPositions[i][0]<InferredPositions[i][1])
                {
                    Sequences2.push_back(Sequences[i].substr(0,InferredPositions[i][0]));
                    Sequences2.push_back(Sequences[i].substr(InferredPositions[i][0]+(LM+estiM), InferredPositions[i][1]-InferredPositions[i][1]-(LM+estiM)));
                    Sequences2.push_back(Sequences[i].substr(InferredPositions[i][1]+(LM+estiM), Sequences[i].size()));
                }
                else
                {
                    Sequences2.push_back(Sequences[i].substr(0,InferredPositions[i][1]));
                    Sequences2.push_back(Sequences[i].substr(InferredPositions[i][1]+(LM+estiM), InferredPositions[i][0]-InferredPositions[i][1]-(LM+estiM)));
                    Sequences2.push_back(Sequences[i].substr(InferredPositions[i][0]+(LM+estiM), Sequences[i].size()));
                }
            }
        }
        
    }
}


void Database::FlipSequences()
{
    vector<string> temp;
    temp = Sequences;
    Sequences.clear();
    Sequences = Sequences2;
    Sequences2 = temp;
    vector <bool> tIndicatorAtLeastTwoMotif;
    IndicatorAtLeastTwoMotif2 = IndicatorAtLeastTwoMotif;
    vector <bool> tIndicatorAtLeastOneMotif;
    IndicatorAtLeastOneMotif2 = IndicatorAtLeastOneMotif;
    NumSequences = int(Sequences.size());
}
