// XMLP
//
// Implementation of the Multi-Layer-Perceptron
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

#include "XMLP.h"

#include <iostream>
#include <cmath>
#include <cstdarg>
#include <cstring>

using namespace std;

XMLP:: XMLP(int layers,double inputRange,string netFile,int innodes,...)
: VSupervisedNet("XMLP",innodes,0,netFile) 
{    
    int I;
    if (layers<1) Errorf((char *)"(XMLP) at least one layer neccessary");
    fParm.fLayers  = layers;
    fParm.fInScale = 1.0/inputRange;
    
    int*    nodes = new int   [fParm.fLayers];      TestPointer(nodes);
    double* step  = new double[fParm.fLayers];      TestPointer(step);
    TNeuralNetParameters::TRANSFER* func  = new TNeuralNetParameters::TRANSFER[fParm.fLayers];      TestPointer(func);
    fPerc           = new TPerceptron*[fParm.fLayers]; TestPointer(fPerc);
    
    va_list ap;
    va_start(ap,innodes);
    
    for (I=0;I<fParm.fLayers;++I) nodes[I] = va_arg(ap,int);
    for (I=0;I<fParm.fLayers;++I) step[I]  = va_arg(ap,double);
    for (I=0;I<fParm.fLayers;++I) func[I]  = (TNeuralNetParameters::TRANSFER) va_arg(ap,int);
    
    va_end(ap);
    
    fParm.fOutNodes = nodes[layers-1];
    
    for (I=0;I<fParm.fLayers;++I) {
        if (I==0)
            fPerc[I] = new TPerceptron(fParm.fInNodes,nodes[I],step[I],func[I],I);
        else
            fPerc[I] = new TPerceptron(fPerc[I-1],nodes[I],step[I],func[I],I);
        
        TestPointer(fPerc[I]);
    }
    fOut = new double[fParm.fOutNodes]; // As we did not know before, allocate here...
    TestPointer(fOut);
    
    delete[] nodes;
    delete[] step;
    delete[] func;
}

XMLP:: XMLP(int layers,double inputRange,string netFile,int innodes,int n0,int n1,int n2,double s0,double s1,double s2,
              TNeuralNetParameters::TRANSFER f0,TNeuralNetParameters::TRANSFER f1,TNeuralNetParameters::TRANSFER f2)
: VSupervisedNet("XMLP",innodes,0,netFile) 
{    
    int I;
    if (layers!=3) Errorf((char *)"(XMLP) Constructor needs 3 Layers");
    fParm.fLayers  = layers;
    fParm.fInScale = 1.0/inputRange;
    
    int*    nodes = new int   [fParm.fLayers];      TestPointer(nodes);
    double* step  = new double[fParm.fLayers];      TestPointer(step);
    TNeuralNetParameters::TRANSFER* func  = new TNeuralNetParameters::TRANSFER[fParm.fLayers];      TestPointer(func);
    fPerc           = new TPerceptron*[fParm.fLayers];  TestPointer(fPerc);
    
    nodes[0] = n0;nodes[1] = n1;nodes[2] = n2;
    step[0]  = s0;step[1]  = s1;step[2]  = s2;
    func[0]  = f0;func[1]  = f1;func[2]  = f2;
    
    fParm.fOutNodes = nodes[layers-1];
    
    for (I=0;I<fParm.fLayers;++I) {
        if (I==0)
            fPerc[I] = new TPerceptron(fParm.fInNodes,nodes[I],step[I],func[I],I);
        else
            fPerc[I] = new TPerceptron(fPerc[I-1],nodes[I],step[I],func[I],I);
        
        TestPointer(fPerc[I]);
    }
    fOut = new double[fParm.fOutNodes]; // As we did not know before, allocate here...
    TestPointer(fOut);
    
    delete[] nodes;
    delete[] step;
    delete[] func;
}

void XMLP:: AllocNet(void) 
{
    fPerc = new TPerceptron*[fParm.fLayers]; TestPointer(fPerc); // MK: Allocation
    int I;
    for (I=0;I<fParm.fLayers;++I) {
        if (I==0)
            fPerc[I] = new TPerceptron();
        else
            fPerc[I] = new TPerceptron(fPerc[I-1]);
        
        TestPointer(fPerc[I]);
    }
}

XMLP::~XMLP() 
{
    if (fFilename!="") if (fShouldSave) Save();
    int I;
    for (I=0;I<fParm.fLayers;++I) delete fPerc[I];
    delete[] fPerc;
}

void XMLP::ReadBinary(void) 
{
    int I;
    fread(&fParm,sizeof(TNeuralNetParameters),1,fFile);
    AllocNet();
    for (I=0;I<fParm.fLayers;++I) {
        fPerc[I]->SetFile(fFile);
        fPerc[I]->ReadBinary();
    }
}

void  XMLP::ReadText(void) 
{
    int I;
    int layers;
    double scale;
    fscanf(fFile,"layers    %i\n",&layers);
    fParm.fLayers = layers;
    fscanf(fFile,"in_scale  %le\n",&scale);
    fParm.fInScale = scale;
    AllocNet();
    for (I=0;I<fParm.fLayers;++I) {
        fPerc[I]->SetFile(fFile);
        fPerc[I]->ReadText();
    }
}

void XMLP::WriteBinary(void) 
{
    int I;
    fwrite(&fParm,sizeof(TNeuralNetParameters),1,fFile);
    for (I=0;I<fParm.fLayers;++I) {
        fPerc[I]->SetFile(fFile);
        fPerc[I]->WriteBinary();
    }
}

void  XMLP::WriteText(void) 
{
    int I;
    fprintf(fFile,"layers    %i\n",fParm.fLayers);
    fprintf(fFile,"in_scale  %le\n",fParm.fInScale);
    for (I=0;I<fParm.fLayers;++I) {
        fPerc[I]->SetFile(fFile);
        fPerc[I]->WriteText();
    }
}


double* XMLP::Recall(NNO_INTYPE* in,NNO_OUTTYPE* out) 
{
    int I;
    
    // convert input
    NNO_INTYPE* i=in;
    double* pi = fPerc[0]->fIn;
    for (I=0;I<fParm.fInNodes;++I) *pi++ = *i++ * fParm.fInScale;
    
    // recallstep of each perceptron
    for (I=0;I<fParm.fLayers;++I) fPerc[I]->Recall();
    for (I=0;I<fParm.fOutNodes;++I) fOut[I] = fPerc[fParm.fLayers-1]->fOut[I];
    
    return fOut;
}

double XMLP::Train(NNO_INTYPE* in,NNO_OUTTYPE* trout) 
{
    int I,J;
    fShouldSave = true;
    
    // convert input
    NNO_INTYPE* i=in;
    double* pi = fPerc[0]->fIn;
    for (I=0;I<fParm.fInNodes;++I) *pi++ = *i++ * fParm.fInScale;
    
    // recallstep of each perceptron
    for (I=0;I<fParm.fLayers;++I) fPerc[I]->Recall();
    for (I=0;I<fParm.fOutNodes;++I) fOut[I] = fPerc[fParm.fLayers-1]->fOut[I];
    
    double S_Err = 0;
    double*   d = fPerc[fParm.fLayers-1]->fDiffSrc;
    double* out = fPerc[fParm.fLayers-1]->fOut;
    NNO_OUTTYPE* tr_out=trout;
    for (J=0;J<fParm.fOutNodes;++J) {
        *d = *tr_out++ - *out++;
        S_Err += *d * *d;
        d++;
    }
    
    for (I=fParm.fLayers-1;I>=0;--I) fPerc[I]->Train();
    
    //if (fPlotter) fPlotter->AddTrainSample(trout[0],trout[0]>fParm.fThreshold);
    
    return S_Err;
}

void XMLP::SetMomentumTerm(double f)
{
    for (int I=0;I<fParm.fLayers;++I) {
        fParm.fMu = f;
        TNeuralNetParameters &parm = fPerc[I]->GetParameters();
        parm.fMu = f;
    }
}

// VNeuralNet
//
// Base classes for unsupervised and supervised networks
// Partof the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

static const char* NNO_VERSION="2.0ROOT";

#include <iostream>
#include <cstdlib>
using namespace std;

TNeuralNetParameters::TNeuralNetParameters()
{
    for (int i=0;i<9;i++) fNetId[i] = 0;
    fLayers = 0;
    fInScale = 1.0;
    fInNodes = 0;
    fOutNodes = 0;
    fLearnStep = 0.01;
    fTransferId = TNeuralNetParameters::TR_FERMI;
    fPerceptronId = 0;
    fThreshold = 0.0;
    fMu = 0.0;
    fFse = 0.0;
}

VNeuralNet::VNeuralNet()
: fBalance(false), fOwnPlotter(false), fParm(), fPlotter(0)
{
    fShouldSave = false;
    fFile  = 0;
    fOut = 0;
}

VNeuralNet::VNeuralNet(string netID,int innodes,int outnodes,string netFile)
: fBalance(false), fOwnPlotter(false), fParm(), fPlotter(0)
{
    fFilename = netFile;
    strncpy(fParm.fNetId,netID.data(),9);
    fShouldSave  = true;
    fFiletype    = FILE_TEXT;
    fParm.fInNodes  = innodes;
    fParm.fOutNodes = outnodes;
    fFile   = 0;
    fOut    = 0;
    if (outnodes>0) {
        fOut = new double[fParm.fOutNodes];
        TestPointer(fOut);
    }
}

VNeuralNet::VNeuralNet(string netFile)
: fBalance(false), fOwnPlotter(false), fParm(), fPlotter(0)
{
    fFilename   = netFile;
    fShouldSave = false;
    fFile    = 0;
    fOut        = 0;
}

VNeuralNet::~VNeuralNet()
{
    if (fOut!=0) { delete[] fOut; fOut = 0; }
    //if (fOwnPlotter) { delete fPlotter; fPlotter = 0; }
}

void VNeuralNet::Save()
{
    fFile = fopen(fFilename.data(),"wb");
    WriteNet();
}

void VNeuralNet::Save(string file)
{
    fFile = fopen(file.data(),"wb");
    WriteNet();
}

void VNeuralNet::WriteNet()
{
    char ftype[16];
    if (fFiletype==FILE_BINARY) strcpy(ftype,"binary"); else strcpy(ftype,"text");
    if (fFile==0) { cerr << "VNeuralNet::WriteNet:: Could not open for writing " << fFilename << endl; return; }
    fprintf(fFile,"C++  NEURAL NETWORK OBJECTS   VERSION %s\nFiletype %s\n",NNO_VERSION,ftype);
    if (fFiletype==FILE_BINARY) WriteNetBinary(); else WriteNetText();
    if (fFiletype==FILE_BINARY) WriteBinary();     else WriteText();
    fclose(fFile);
}

void VNeuralNet::WriteNetText()
{
    fprintf(fFile,"\nnetwork id  %s\n",fParm.fNetId);
    fprintf(fFile,"innodes     %i\n",fParm.fInNodes);
    fprintf(fFile,"outnodes    %i\n",fParm.fOutNodes);
}

void VNeuralNet::WriteNetBinary()
{
    fwrite(&fParm,sizeof(TNeuralNetParameters),1,fFile);
}

void VNeuralNet::ReadNet(const char* netID)
{
    fFile = fopen(fFilename.data(),"rb");
    if (fFile==0) Errorf((char *)"file %s not found",(char *)fFilename.data());
    char ftype[16];
    char Version[16];
    fscanf(fFile,"C++  NEURAL NETWORK OBJECTS   VERSION %s\nFiletype %s\n",Version,ftype);
    if      (!strcmp(ftype,"binary")) fFiletype = FILE_BINARY;
    else if (!strcmp(ftype,"text"))   fFiletype = FILE_TEXT;
    else Errorf((char *)"illegal fileformat: %s",(char *)fFilename.data());
    
    if (fFiletype==FILE_BINARY)
        ReadNetBinary();
    else
        ReadNetText();
    
    fParm.fNetId[4]=0;
    if (strcmp(netID,fParm.fNetId)) {
        fclose(fFile);
        Errorf((char *)"file %s  (incompatible network)\nnetwork ID is %s and should be %s",fFilename.data(),fParm.fNetId,netID);
    }
    
    if (strcmp(Version,NNO_VERSION)) {
        fclose(fFile);
        Errorf((char *)"illegal NNO version number of file %s\nversion number is %s and should be %s",fFilename.data(),Version,(char *)NNO_VERSION);
    }
    
    if (fFiletype==FILE_BINARY)
        ReadBinary();
    else
        ReadText();
    
    fOut = new double[fParm.fOutNodes];
    
    TestPointer(fOut);
    
    fclose(fFile);
}

void VNeuralNet::ReadNetText()
{
    fscanf(fFile,"\nnetwork id  %s\n",fParm.fNetId);
    fscanf(fFile,"innodes     %i\n",&fParm.fInNodes);
    fscanf(fFile,"outnodes    %i\n",&fParm.fOutNodes);
}

void VNeuralNet::ReadNetBinary()
{
    fread(&fParm,sizeof(TNeuralNetParameters),1,fFile);
}

void VNeuralNet::Errorf(char* format,...)
{
    va_list ap;
    va_start(ap, format);
    char MainFormat[256];
    sprintf(MainFormat,"NNO ERROR: %s\n",format);
    vfprintf(stderr,MainFormat,ap);
    exit(1);
}

void VNeuralNet::Warningf(FILE* f,char* format,...)
{
    va_list ap;
    va_start(ap, format);
    char MainFormat[256];
    sprintf(MainFormat,"NNO WARNING: %s\n",format);
    vfprintf(f,MainFormat,ap);
}

void VNeuralNet::Messagef(FILE* f,char* format,...)
{
    va_list ap;
    va_start(ap, format);
    char MainFormat[256];
    sprintf(MainFormat,"NNO INFO: %s\n",format);
    vfprintf(f,MainFormat,ap);
}

double VNeuralNet::Random(void)
{
#ifdef NNORAND
    //  Machine independent random number generator.
    //  Produces uniformly-distributed floating points between 0 and 1.
    //  Identical sequence on all machines of >= 32 bits.
    //  Universal version (Fred james 1985).
    //  Return numbers in the range -0.5..0.5 (MK)
    
    const float kCONS = 4.6566128730774E-10;
    const int kMASK31 = 2147483647;
    static unsigned int fSeed = 65539;
    
    fSeed *= 69069;
    // keep only lower 31 bits
    fSeed &= kMASK31;
    // Set lower 8 bits to zero to assure exact float
    int jy = (fSeed/256)*256;
    double random = kCONS*jy;
    return 1.0 - 2.0*random;
#else
    return rand() / (RAND_MAX + 1.);
#endif
}

void VNeuralNet::TestPointer(void* Ptr)
{
    if (Ptr==0) Errorf((char *)"not enough memory");
}

void VNeuralNet::SetupPlots(VNeuralNetPlotter *plotter)
{
/*
    if (plotter==0) {
        cout << "Instantiating plotter for " << GetName() << endl;
        if (fOwnPlotter) delete fPlotter;
        fPlotter = new TSimpleNeuralNetPlotter(GetName());
        fOwnPlotter = true;
    }
    
    fPlotter->Initialize();
*/
}

void VNeuralNet::FillPlots(double trn, double tst)
{
    if (fPlotter==0) return;
    //fPlotter->AddTrainSample(trn,true);
    //fPlotter->AddTestSample(tst,true);
}

void VNeuralNet::ShowPlots()
{
    if (fPlotter==0) return;
    //fPlotter->ShowPlots();
}

double VNeuralNet::TrainEpoch(TDataServe *server, int nEpoch)
{
/*
    double       error = 0.0;            // squared error collector
    unsigned int classError;        // classification Error
    unsigned int n;            // number of samples
    
    const int samples = server->GetNumTrnvecs();
    const int tests   = server->GetNumTstvecs();
    
    for (int epo=0; epo<nEpoch; epo++){
        error = 0.0;
        classError = 0;
        n = 0;
        
        server->MixTrn(); // Shuffle the dataset
        
        for (int i=0; i<samples; i++){
            
            int trnind = i;
            if (fBalance) trnind = BalancedTrnIndex(server);
            
            float *inv  = (float *) server->GetInvecTrn(trnind);
            float *outv = (float *) server->GetOutvecTrn(trnind);
            
            error += Train(inv,outv);
            n++;
        }
        
        classError = (unsigned int) TestEpoch(server);
        double percentage = 100. * classError;
        if (tests>0) percentage /= tests; else percentage = 0;
        
        // print training info
        cout << GetNetID() << ": Epoch " << epo <<
        ", samples " << n <<
        ", Error " << error <<
        ", classError " << classError <<
        " (" << (int)percentage << "%)" << endl;
        
        // Fill the plots (for a random recall)
        if (fPlotter!=0) {
            fPlotter->AddTrainGraph(error);
            fPlotter->AddTestGraph(classError);
            fPlotter->ShowPlots();
            fPlotter->Reset();
        }
    }
    
    if (fShouldSave) Save(); // Store the net
    
    return error;
*/
    return 0.0;
}

// Check the network performance

double VNeuralNet::TestEpoch(TDataServe *server)
{
/*
    unsigned int classError = 0;    // classification Error
    const int samples = server->GetNumTstvecs();
    TNeuralNetParameters &parm = GetParameters();
    
    int i;
    
    for (i=0; i<samples; i++){
        
        int tstind = i;
        if (fBalance) tstind = BalancedTstIndex(server);
        
        float *inv  = server->GetInvecTst(tstind);
        float *outv = server->GetOutvecTst(tstind);
        
        // compare network recall with server
        Recall(inv,outv);
        
        for (int ii=0;ii<parm.fOutNodes;++ii) {
            double answer = GetOutput()[ii];
            if ((answer>parm.fThreshold && outv[ii]<=parm.fThreshold) ||
                (answer<=parm.fThreshold && outv[ii]>parm.fThreshold) )
                ++classError; // classification ok ?
        }
        
    }
    
    return classError;
*/
    return 0.0;
}

double  VNeuralNet::Test(NNO_INTYPE* in,NNO_OUTTYPE* out)
{
    Recall(in,out);
    int I;
    double* o = fOut;
    double diff,totalError = 0.0;
    for (I=0;I<fParm.fOutNodes;++I) {
        diff = *out++ - *o++;
        totalError += diff * diff;
    }
    return totalError;
}

double VNeuralNet::TrainEpoch(string file, int nEpoch)
{
    FILE* ftrn=fopen(file.data(),"rb");
    if (ftrn==0) {
        cerr << "Training file does not exist:" << file << endl;
        return 0.0;
    }
    
    int epoch=0;        // epoch counter
    double error;        // squared error collector
    unsigned int classError;        // classification Error
    
    NNO_INTYPE   *in  = new NNO_INTYPE[fParm.fInNodes];        // inputvector
    NNO_OUTTYPE  *out = new NNO_OUTTYPE[fParm.fOutNodes];   // outputvector
    
    // begin of training
    do {
        error      = 0.0;
        classError = 0;
        
        while ( fread(in,sizeof(NNO_INTYPE),fParm.fInNodes,ftrn) ) {  // read inputvector
            fread(out,sizeof(NNO_OUTTYPE),fParm.fOutNodes,ftrn);      // read outputvector
            //double output = out[0];
            
            error += Train(in,out);                // perform learnstep
            
            // compare network output with 'Out'
            double *net = GetOutput();
            //double answer = net[0];
            for (int I=0;I<fParm.fOutNodes;++I) {
                if ((net[I]>fParm.fThreshold && out[I]<=fParm.fThreshold) ||
                    (net[I]<=fParm.fThreshold && out[I]>fParm.fThreshold) )
                    ++classError; // classification ok ?
            }
            
        }
        rewind(ftrn);  // epoch completed, rewind filepointer
        ++epoch;
        
        // print training info
        cout << GetNetID() << ": Epoch " << epoch << ", Error " << error << ", classError " << classError << endl;
        
        // Fill the plots (for a random recall)
        /*
        if (fPlotter!=0) {
            fPlotter->AddTrainGraph(error);
            fPlotter->AddTestGraph(classError);
            fPlotter->ShowPlots();
            fPlotter->Reset();
        }
         */
        
    } while (classError>0 && epoch<nEpoch);
    
    fclose(ftrn);
    
    delete[] in; delete[] out;
    
    if (fShouldSave) Save(); // Store the net
    
    return error;
}

double VNeuralNet::TestEpoch(string file)
{
    FILE* ftst=fopen(file.data(),"rb");
    if (ftst==0) {
        cerr << "Test file does not exist:" << file << endl;
        return 0.0;
    }
    
    unsigned int classError = 0;    // classification Error
    NNO_INTYPE   *in  = new NNO_INTYPE[fParm.fInNodes];        // inputvector
    NNO_OUTTYPE  *out = new NNO_OUTTYPE[fParm.fOutNodes];   // outputvector
    
    while ( fread(in,sizeof(NNO_INTYPE),fParm.fInNodes,ftst) ) {   // read inputvector
        fread(out,sizeof(NNO_OUTTYPE),fParm.fOutNodes,ftst);            // read outputvector
        
        // compare network recall with file
        double *net = Recall(in,out);
        for (int I=0;I<fParm.fOutNodes;++I) {
            if ((net[I]>fParm.fThreshold && out[I]<=fParm.fThreshold) ||
                (net[I]<=fParm.fThreshold && out[I]>fParm.fThreshold) )
                ++classError; // classification ok ?
        }
        
    }
    
    fclose(ftst);
    
    delete[] in; delete[] out;
    
    return classError;
}

unsigned int VNeuralNet::BalancedTrnIndex(TDataServe *server)
{
    /*
    static ULong_t ngood=0, nbad=0;
    unsigned int samples = server->GetNumTrnvecs();
    unsigned int index = (unsigned int) (rand()%samples);
    float *outv = server->GetOutvecTrn(index);
    
    if (ngood<nbad)
        while (outv[0]<=0.0) {
            index++;
            outv = server->GetOutvecTrn(index%samples);
        }
    
    if (outv[0]>0.0) ngood++; else nbad++;
    
    return index%samples;
     */
    return 0;
}

unsigned int VNeuralNet::BalancedTstIndex(TDataServe *server)
{
    /*
    static ULong_t ngood=0, nbad=0;
    unsigned int samples = server->GetNumTstvecs();
    unsigned int index = (unsigned int) (rand()%samples);
    float *outv = server->GetOutvecTst(index);
    
    if (ngood<nbad)
        while (outv[0]<=0.0) {
            index++;
            outv = server->GetOutvecTst(index%samples);
        }
    
    if (outv[0]>0.0) ngood++; else nbad++;
    
    return index%samples;
     */
    return 0;
}

void VNeuralNet::SetMomentumTerm(double f)
{
    fParm.fMu = f;
}

void VNeuralNet::SetFlatSpotElimination(double f)
{
    fParm.fFse = f;
}

// VSupervisedNet
//
// Base classes for supervised learning
// Abstract base class of all unsupervised networks
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

long VSupervisedNet::TrainEpoch(float *tuple, bool randomize) {
    /*
    fTuple = tuple;
    if (fTuple == 0) return 0;
    long nhits = fTuple->GetEntries();
    for (int i=0;i<nhits;i++) {
        Long_t index = i;
        if (randomize) index = rand()%nhits;
        fTuple->GetEvent(index,1);
        float *x=fTuple->GetArgs();
        Learnstep(x, &x[fParm.fInNodes]); // the first fInNodes columns hold input data, the following fOutNodes columns hold the output data
    }
    return nhits;
     */
    return 0;
}

// TPerceptron
//
// Implementation of the perceptron (Supervised Learning)
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

// Transferfunctions
void TransferFermi(double in,double* out,double* deriv)
{
    if (in < -10.0) in = -10.0;
    if (in >  10.0) in =  10.0;
    double O = 1. / (1. + exp(-in));
    *out   = O;
    *deriv = O * (1. - O);
}

void TransferSigmoid(double in,double* out,double* deriv)
{
    if (in < -10.0) in = -10.0;
    if (in >  10.0) in =  10.0;
    double O = 1. - 2. / (1. + exp(-in));
    *out   = O;
    *deriv = 2. * O * (1. - O);
}

void TransferLinear(double in,double* out,double* deriv)
{
    *out = in;
    *deriv = 1.0;
}

void TransferLinearBend(double in,double* out,double* deriv)
{
    if (in < -1.0) {
        *out = -0.9 + in * 0.1;
        *deriv = 0.1;
    } else
        if (in > 1.0) {
            *out = 0.9 + in * 0.1;
            *deriv = 0.1;
        } else {
            *out = in;
            *deriv = 1.0;
        }
}

TPerceptron::TPerceptron(int inNodes,
                         int outNodes,
                         double learnStep,
                         TNeuralNetParameters::TRANSFER transferId,
                         int perceptronId)
{
    
    fParm.fInNodes      = inNodes;
    fParm.fOutNodes     = outNodes;
    fParm.fLearnStep    = learnStep;
    fParm.fTransferId   = transferId;
    fParm.fPerceptronId = perceptronId;
    fPrev     = 0;
    AllocNet();
    InitNet();
}

TPerceptron::TPerceptron(TPerceptron* prev,
                         int outNodes,
                         double learnStep,
                         TNeuralNetParameters::TRANSFER transferId,
                         int perceptronId)
{
    
    fParm.fInNodes      = prev->fParm.fOutNodes;
    fParm.fOutNodes     = outNodes;
    fParm.fLearnStep    = learnStep;
    fParm.fTransferId   = transferId;
    fParm.fPerceptronId = perceptronId;
    fPrev             = prev;
    AllocNet();
    InitNet();
}

TPerceptron::TPerceptron(void)
{
    fPrev = 0;
}

TPerceptron::TPerceptron(TPerceptron* prev)
{
    fPrev = prev;
}

void TPerceptron::AllocNet(void)
{
    int I;
    fU       = new PerceptronUnit[fParm.fOutNodes]; TestPointer(fU);
    fOut     = new double[fParm.fOutNodes]; TestPointer(fOut);
    fDiffSrc = new double[fParm.fOutNodes]; TestPointer(fDiffSrc);
    if (fPrev == 0) {
        fIn      = new double[fParm.fInNodes]; TestPointer(fIn);
        fDiffDst = 0;
    } else {
        fIn      = fPrev->fOut;
        fDiffDst = fPrev->fDiffSrc;
    }
    fUbound = &fU[fParm.fOutNodes];
    PerceptronUnit* up = fU;
    for (I=0;I<fParm.fOutNodes;++I) {
        up->fVector = new double[fParm.fInNodes];
        TestPointer(up->fVector);
        up->fDelta = new double[fParm.fInNodes];
        TestPointer(up->fDelta);
        up->fThreshold = 0;
        up->fID = I;
        ++up;
    }
    switch (fParm.fTransferId) {
        case TNeuralNetParameters::TR_FERMI  :      Transfer=TransferFermi;      break;
        case TNeuralNetParameters::TR_SIGMOID:      Transfer=TransferSigmoid;      break;
        case TNeuralNetParameters::TR_LINEAR :      Transfer=TransferLinear;     break;
        case TNeuralNetParameters::TR_LINEAR_BEND:  Transfer=TransferLinearBend; break;
        default:             Transfer=0;
    }
}

void TPerceptron::InitNet(void)
{
    PerceptronUnit* up;
    int J;
    for(up=fU;up<fUbound;++up) {
        for (J=0;J<fParm.fInNodes;++J) {
            up->fVector[J] = Random();
            up->fDelta[J]  = 0.0;
        }
        up->fThreshold = Random();
    }
}


TPerceptron ::~TPerceptron()
{
    PerceptronUnit* up = fU;
    if (fU!=0)  {
        for(up=fU;up<fUbound;++up) {
            delete[] up->fVector;
            delete[] up->fDelta;
        }
        delete[] fU;
    }
    delete[] fOut;
    delete[] fDiffSrc;
    if (fPrev==0) delete[] fIn;
}


void TPerceptron::WriteBinary()
{
    PerceptronUnit* up;
    fwrite(&fParm,sizeof(PerceptronBase),1,fFile);
    for(up=fU;up<fUbound;++up) {
        fwrite(up->fVector,sizeof(double),fParm.fInNodes,fFile);
        fwritevar(up->fThreshold);
        fwritevar(up->fID);
    }
}

void TPerceptron::ReadBinary()
{
    PerceptronUnit* up;
    fread(&fParm,sizeof(PerceptronBase),1,fFile);
    AllocNet();
    for(up=fU;up<fUbound;++up) {
        fread(up->fVector,sizeof(double),fParm.fInNodes,fFile);
        freadvar(up->fThreshold);
        freadvar(up->fID);
    }
}

void TPerceptron::WriteText()
{
    fprintf(fFile,"\nPerceptron ID %i\n",fParm.fPerceptronId);
    fprintf(fFile,"innodes     %i\n",fParm.fInNodes);
    fprintf(fFile,"outnodes    %i\n",fParm.fOutNodes);
    fprintf(fFile,"learn_step  %le\n",fParm.fLearnStep);
    fprintf(fFile,"transfer_id %i\n",fParm.fTransferId);
    PerceptronUnit* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        fprintf(fFile,"\n");
        fprintf(fFile,"unit number      %i\n",up->fID);
        fprintf(fFile,"threshold        %le\n",up->fThreshold);
        fprintf(fFile,"weights\n");
        for (I=0;I<fParm.fInNodes;++I) fprintf(fFile,"%le\n",up->fVector[I]);
        fprintf(fFile,"\n");
    }
}

void TPerceptron::ReadText()
{
    fscanf(fFile,"\nPerceptron ID %i\n",&fParm.fPerceptronId);
    fscanf(fFile,"innodes     %i\n",&fParm.fInNodes);
    fscanf(fFile,"outnodes    %i\n",&fParm.fOutNodes);
    fscanf(fFile,"learn_step  %le\n",&fParm.fLearnStep);
    fscanf(fFile,"transfer_id %i\n",(int *)&fParm.fTransferId);
    AllocNet();
    PerceptronUnit* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        fscanf(fFile,"\n");
        fscanf(fFile,"unit number      %i\n",&up->fID);
        fscanf(fFile,"threshold        %le\n",&up->fThreshold);
        fscanf(fFile,"weights\n");
        for (I=0;I<fParm.fInNodes;++I) fscanf(fFile,"%le\n",&up->fVector[I]);
        fscanf(fFile,"\n");
    }
}

double* TPerceptron::Recall(NNO_INTYPE*,NNO_OUTTYPE*)
{
    int I;
    PerceptronUnit* up;
    double* o = fOut;
    double* ds = fDiffSrc;
    if (Transfer==0) Errorf((char *)"(TPerceptron) undefined transferfunction");
    for(up=fU;up<fUbound;++up) {
        double* v = up->fVector;
        double* i = fIn;
        double  sum = 0.0;
        for (I=0;I<fParm.fInNodes;++I) sum += *i++ * *v++;
        sum -= up->fThreshold;
        Transfer(sum,o++,ds++);
    }
    
    return o;
}

double TPerceptron::Train(NNO_INTYPE*,NNO_OUTTYPE*)
{
    int I;
    PerceptronUnit* up;
    double* ds;
    
    // modify weights
    ds = fDiffSrc;
    for(up=fU;up<fUbound;++up) {
        double* i = fIn;
        double* v = up->fVector;
        double* m = up->fDelta;
        for (I=0;I<fParm.fInNodes;++I) {
            double delta = *i++ * *ds;
            *v++ += (delta + (*m * fParm.fMu)) * fParm.fLearnStep;
            *m++ = delta;
        }
        up->fThreshold -= *ds * fParm.fLearnStep;
        ++ds;
    }
    
    // propagate derivation backward if previous perceptron exists
    if (fDiffDst!=0) {
        double diff;
        double* dd = fDiffDst;
        for (I=0;I<fParm.fInNodes;++I) {
            diff = 0.0;
            ds = fDiffSrc;
            for(up=fU;up<fUbound;++up) diff += up->fVector[I] * *ds++;
            *dd++ *= diff;
        }
    }
    
    return 0.0;
}


