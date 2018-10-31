#ifndef DATABASE_H
#define DATABASE_H
#include <vector>
#include <string>
using namespace std;

class Database 
{
public:
    Database();
    Database(string);
    Database(ifstream&);
    virtual ~Database();
    void insertSequence(string);
    int getSize();
    void showPWM();
    void uninformativePWM();
    void ComputeParticleFilter(ofstream&);
    void DetectionZeroAgainstMany();
    void DetectionOneAgainstMore();
    void ComputeFinalSolution(ofstream&);
    void EliminateInstances();
    void FlipSequences();
    void SetNumberOfParticles(int n) {NumberOfParticles=n;}
    int GetNumberOfParticles() {return NumberOfParticles;}
    void SetRho0(int n) {InitialRho0=n;}
    void SetRho1(int n) {InitialRho1=n;}
    void SetSigma0(int n) {InitialSigma0=n;}
    void SetSigma1(int n) {InitialSigma1=n;}
    void SetSigma2(int n) {InitialSigma2=n;}
    void SetSymOption(int n) {SymOption=n;}
    void SetNullDist(float nulldist[]){
        for(int i=0; i<4; i++){
            NullDist.push_back(nulldist[i]);
        }
    }
    void SetEstiM(int n) {estiM = n;}
    void SetEstiB(int n) {estiB = n;}
    void SetEstiF(int n) {estiF = n;}
    void SetEstiP(int n) {estiP = n;}
    void SetNclasses(int n) {nclasses = n;}
    void SetLM(int n) {LM = n;}
    void SetUM(int n) {UM = n;}
    void SetLB(int n) {LB = n;}
    void SetUB(int n) {UB = n;}
    void SetLF(int n) {LF = n;}
    void SetUF(int n) {UF = n;}
protected:
private:
    vector<float> NullDist;
    int NumberOfParticles;
    int InitialRho0;
    int InitialRho1;
    int InitialSigma0;
    int InitialSigma1;
    int InitialSigma2;
    int SymOption;
    int NumSequences;
    vector<string> Sequences;
    vector<string> oSequences;
    vector<string> Sequences2;
    int estiM;
    int estiB;
    int estiF;
    int estiP;
    vector<vector <double> > PWM;
    vector<vector<vector<vector<vector<double> > > > > ALPHA;
    int sigma0est;
    int sigma1est;
    int sigma2est;
    int nclasses;
    int LM;
    int UM;
    int LB;
    int UB;
    int LF;
    int UF;
    vector<int> InferredPositions1;
    vector<vector <int> > InferredPositions;
    vector <bool> IndicatorAtLeastOneMotif;
    vector <bool> IndicatorAtLeastOneMotif2;
    vector <bool> IndicatorAtLeastTwoMotif;
    vector <bool> IndicatorAtLeastTwoMotif2;
    float pmfsim(const vector<float>&, const vector<double>&);
    float pmfsimL(const vector<float>&, const vector<long double>&);
    long double ComputePosterior1Motif (const string &, const int, const vector<vector<float> > &, const vector<vector<float> > &, const int);
    long double ComputePosterior2Motifs (const string &, const string &, const vector<vector<float> > &, const vector<vector<float> > &, const int);
    long double FindML (const string &, const vector<vector<float> > &, const vector<vector<float> > &, const int, vector<int> &);
    string ReverseComplement(const string);
};

#endif // DATABASE_H
