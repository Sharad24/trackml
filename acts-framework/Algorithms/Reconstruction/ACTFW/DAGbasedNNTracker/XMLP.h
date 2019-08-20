#ifndef XMLP_H
#define XMLP_H
// XMLP
//
// Implementation of the Multi-Layer-Perceptron
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

#include <string>

// Define the precision
typedef float NNO_INTYPE;
typedef float NNO_OUTTYPE;

class VNeuralNetPlotter;
class TDataServe;
class TTree;
class TGraph;
class TCanvas;

// Base struct of all networks
class TNeuralNetParameters {
public:
    // ID of transferfunction
    enum TRANSFER {TR_USER,TR_FERMI,TR_LINEAR,TR_LINEAR_BEND,TR_SIGMOID};
    
    char fNetId[9];
    int    fLayers;        // number of perceptron layers
    double fInScale;        // scale input vector
    int    fInNodes;        // number of input nodes
    int    fOutNodes;        // number of output nodes
    double fLearnStep;    // learning step
    double fMu;        // momentum term
    double fFse;        // flat spot elimination
    TRANSFER fTransferId;   // transfer function
    int fPerceptronId;    // ID of perceptron
    double fThreshold;    // Threshold for output
public:
    TNeuralNetParameters();
    virtual ~TNeuralNetParameters() {}
};

// Base class of all networks
class VNeuralNet {
public:
    // Abstract interface for all networks
    
    virtual void AllocNet() = 0;
    virtual void InitNet() = 0;
    virtual void WriteText() = 0;
    virtual void WriteBinary() = 0;
    virtual void ReadText() = 0;
    virtual void ReadBinary() = 0;
    virtual double* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out=0) = 0;
    virtual double Train(NNO_INTYPE* in,NNO_OUTTYPE* out=0) = 0;    // returns squared error
    
    // For backwards compatibility
    double*    Recallstep(NNO_INTYPE* in,NNO_OUTTYPE* out=0) { return Recall(in,out); }
    double    Learnstep(NNO_INTYPE* in,NNO_OUTTYPE* out=0) { return Train(in,out); }
    
    // Training and testing
    double    TrainEpoch(TDataServe *server, int nEpoch=1);
    double    TestEpoch(TDataServe *server);
    double    TrainEpoch(std::string file, int nEpoch=1);
    double    TestEpoch(std::string file);
    double    Test(NNO_INTYPE* in,NNO_OUTTYPE* trn);        // returns squared error
    void    BalanceSamples(bool yesNo = true) { fBalance = yesNo; }
    virtual void SetMomentumTerm(double f);
    virtual void SetFlatSpotElimination(double f);
    
    enum FILE_TYPE {FILE_BINARY,FILE_TEXT,FILE_ROOT};
    
private:
    bool  fBalance;            //!Balance positive and negative samples
    bool  fOwnPlotter;        //!Ownership of plotter
    
    void    WriteNetText();
    void    WriteNetBinary();
    void    ReadNetText();
    void    ReadNetBinary();
    void    WriteNet();
    unsigned int  BalancedTrnIndex(TDataServe *server);
    unsigned int  BalancedTstIndex(TDataServe *server);
    
protected:
    std::string        fFilename;        // Name of network file
    FILE_TYPE        fFiletype;        // Type of file
    TNeuralNetParameters  fParm;    // Topology of network
    VNeuralNetPlotter* fPlotter;    //!Show plots
    double*        fOut;        //!Output nodes
    bool        fShouldSave;    //!Save the network
    FILE*        fFile;        //!File pointer
    
    // Help functions
    void freadvar(NNO_INTYPE Var) { fread(&Var,sizeof(Var),1,fFile); }
    void fwritevar(NNO_OUTTYPE Var) { fwrite(&Var,sizeof(Var),1,fFile); }
    void ReadNet(const char* netID);
    void TestPointer(void* ptr);
    void Errorf(char* format,...);
    void Warningf(FILE* f,char* format,...);
    void Messagef(FILE* f,char* format,...);
    
public:
    VNeuralNet();
    VNeuralNet(std::string netID,int innodes,int outnodes,std::string netFile);
    VNeuralNet(std::string netFile);
    virtual ~VNeuralNet();
    void Save(std::string file);
    void Save();
    
    // Getter functions
    TNeuralNetParameters& GetParameters() { return fParm; }
    VNeuralNetPlotter& GetPlotter() const { return *fPlotter; }
    std::string     GetFilename() const { return fFilename; }
    double*        GetOutput() const { return fOut; }
    double        GetThreshold() const { return fParm.fThreshold; }
    std::string     GetNetID() const { return fParm.fNetId; }
    static double Random(void);
    
    // Setter functions
    void SetThreshold(double t) { fParm.fThreshold = t; }
    void SetPlotter(VNeuralNetPlotter *plotter) { fPlotter = plotter; }
    
    // Plotter functions
    void SetupPlots(VNeuralNetPlotter *plotter=0);
    void FillPlots(double trn=0.0, double tst=0.0);
    void ShowPlots();
};

class VSupervisedNet : public VNeuralNet {
protected:
    float *fTuple;  // Training data
public:
    VSupervisedNet() : VNeuralNet() {};
    VSupervisedNet(std::string netID,int innodes,int outnodes,std::string netFile) :
    VNeuralNet(netID,innodes,outnodes,netFile) {}
    VSupervisedNet(std::string netFile) :
    VNeuralNet(netFile) {}
    virtual ~VSupervisedNet() {};
    virtual long TrainEpoch (float *tuple, bool random=true);    // learn the hits from the ntuple
};

class PerceptronUnit {
public:
    double*  fVector;    // synaptic vector
    double*  fDelta;       // modification vector
    double   fThreshold; // activity threshold
    int      fID;        // ID of this element
};

class PerceptronBase {
public:
};

class TPerceptron : public VSupervisedNet {
public:
    TPerceptron(int fInNodes,       // constructor for first perceptron
                int fOutNodes,
                double fLearnStep,
                TNeuralNetParameters::TRANSFER fTransferId,
                int fPerceptronId);
    
    TPerceptron(TPerceptron* Prev,    // constructor for linked perceptron
                int fOutNodes,
                double learnstep,
                TNeuralNetParameters::TRANSFER fTransferId,
                int fPerceptronId);
    
    TPerceptron(void);                // constructor for first perceptron (data from file)
    
    TPerceptron(TPerceptron* Prev);   // constructor for linked perceptron (data from file)
    
    virtual ~TPerceptron();           // destructor of network
    
    PerceptronUnit* fU;    //! Temp. unit
    void (*Transfer)(double in,double* out,double* deriv);    //! Transfer function
    
    double* fIn;   //!inputvector:  in[fInNodes]
    double* fOut;  //!outputvector: out[fOutNodes];
    double* fDiffSrc; //!error derivation from following net: fDiffSrc[fOutNodes]
    double* fDiffDst; //!error derivation to previous net:    fDiffDst[fInNodes]
    
private:
    TPerceptron* fPrev;    //! previous perceptron
    PerceptronUnit* fUbound;    //! Temp. unit
    void AllocNet(void);
    void InitNet(void);
    
public:
    void SetFile(FILE *f) { fFile = f; }
    virtual void WriteText();
    virtual void WriteBinary();
    virtual void ReadText();
    virtual void ReadBinary();
    
    double Train(NNO_INTYPE* in=0,NNO_OUTTYPE* out=0);  // calling Learnstep, Recallstep must have been performed already
    double* Recall(NNO_INTYPE* in=0,NNO_OUTTYPE* out=0);
    //void CopyData(const TPerceptron& PERC); // copies data from another perceptron
};

class XMLP : public VSupervisedNet {
private:
    void AllocNet(void);
    void InitNet(void) {};
    void WriteText(void);
    void WriteBinary(void);
    void ReadText(void);
    void ReadBinary(void);
public:
    XMLP() {};
    XMLP(int layers,double inputRange,std::string netFile,int innodes,...);
    XMLP(int layers,double inputRange,std::string netFile,int innodes,int n1,int n2,int n3,double s1,double s2,double s3,
          TNeuralNetParameters::TRANSFER f1,TNeuralNetParameters::TRANSFER f2,TNeuralNetParameters::TRANSFER f3);
    XMLP(std::string netFile) : VSupervisedNet(netFile) {ReadNet("XMLP");};
    virtual ~XMLP();
    double Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    double* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out);
    TPerceptron* GetPerceptron(int i) { return fPerc[i]; }
    virtual void SetMomentumTerm(double f);
    
    TPerceptron** fPerc;    //! Temp.unit
};

#endif

