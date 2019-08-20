// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"
#include "Reconstruction.h"
#include "PolarModule.h"
#include "Graph.h"
#include "parallel.h"

#ifdef USETMVA
#ifdef TMVAREADER
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif
#include "TMVAClassification_MLP1.h"
#include "TMVAClassification_MLP2.h"
#include "TMVAClassification_MLP3.h"
#else
#include "XMLP.h"
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <map>
#include <vector>
#include <cmath>
#include <stack>
#include <queue>
#include <mutex>
#ifdef _WIN32
#include "direntwin32.h"
#else
#include <dirent.h>
#endif

using namespace std;

PolarModule *mod[LAYERS];
map<int, Detector> Tracker::detectors;
bool Tracker::eval;
bool Tracker::debug;

Graph<long long> Tracker::pgraph[NSLICES];   // graph to represent path information trained on slices
map<std::pair<int,int>, Graph<long long> > Tracker::tgraph; // graph to represent path information trained on individual tracks

#define PARALLEL_FOR_BEGIN(nb_elements) parallel_for(nb_elements, [&](int start, int end){ for(int i = start; i < end; ++i)
#define PARALLEL_FOR_END()},MULTITHREAD)

void process_mem_usage(double& vm_usage, double& resident_set);


Tracker::Tracker(int n,const char *datapath,const char *workpath) : filenum(n), dataPath(datapath), workPath(workpath), threshold1(THRESHOLD1), threshold2(THRESHOLD2), threshold3(THRESHOLD3), cscore(0.0)
{
    static bool initialized(false);
#ifdef EVAL
    eval = true;
#else
    eval = false;
#endif
    
    cout << endl << "Instantiate Tracker " << filenum <<  "  Evaluation: " << eval << " Threaded: " << MULTITHREAD << endl;
    cout << "datapath: " << dataPath << endl;
    cout << "workpath: " << workPath << endl;
    
    if (!initialized) {
        cout << "PHIRESOLUTION   " << Reconstruction::phiResolution() << endl;
        cout << "THETARESOLUTION " << Reconstruction::thetaResolution() << endl;
        cout << "SWEIGHT         " << SWEIGHT << endl;
        cout << "VERTEXCUTXY     " << VERTEXCUTXY << endl;
        cout << "VERTEXCUTZ      " << VERTEXCUTZ << endl;
        cout << "VERTEXCUT       " << VERTEXCUT << endl;
        
        string filename = workPath + "/detectors.csv";
        readDetectors(filename.c_str());
        
        cout << "Reading slice graphs from " << SLICEDIR << endl;
        readSliceGraphs();
        
        // Read tracking graphs
        cout << "Reading tile graphs from " << TILESDIR << endl;
        string dirname = workPath + "/" + TILESDIR;
        readTileGraphs(dirname.c_str());
        
        // Set up the neural networks
        initNeuralNetworks();
        
        initialized = true;
    }
    
    for (int i=0;i<MAXDIM;i++) {
        assignment[i] = 0;
    }
    
    memoryFootprint();
}


Tracker::~Tracker()
{
    cout << "Delete Tracker " << filenum << endl;
    
    for (int i = 0; i < LAYERS; i++) {
        if (mod[i] != NULL) {
            delete mod[i];
            mod[i] = NULL;
        }
    }
    
    for (int i = 0; i < LAYERS; i++) {
    }
}


// function to run the complete tracking algorithm
int* Tracker::findTracks(int step, const char *workpath)
{
    cout << "findTracks working directory: " << workpath << endl;
    workPath = workpath;
    double accuracy = 0.0;
    
    if (!eval) timer(-1);
    
    long nhits(0), ncells(0);
    nhits  = initGraphData();
    ncells = initHitDir();
    cout << "initGraphData: " << nhits << "  initHitDir: " << ncells << endl;
    
    initPolarModule();
    
    memoryFootprint();
    
    // Instantiate reconstruction objects for parallel processing
    static const int n = NTHREADS;
    initRecoObjects(n);
    
    // Generate a list of tasks
    initTasks();
    
    memoryFootprint();
    
    vector<triple> pairs[n];
    vector<triple> triples[n];
    vector<vector<int> > path[n];
    vector<vector<int> > paths;
    
    static mutex mylock;
    
    //for (int i=0;i<n;i++)
    PARALLEL_FOR_BEGIN(n)
    {
        while (!tasks.empty()) {
            
            // get a task from the stack (beware of race hazards)
            mylock.lock();
            if (tasks.empty()) break;
            int task = tasks.top().first;
            auto *g = tasks.top().second;
            tasks.pop();
            cout << "thread " << i << " task " << task << endl;
            mylock.unlock();
                        
            // Find promising pairs
            auto threadPairs = reco[i]->findCandidatesGraph(*g); // Pair candidates
            auto threadTriples = reco[i]->findTriplesGraph(*g,threadPairs); // Extend pairs to triples
            
            // Extend promising triples to paths
            auto threadPath = reco[i]->findPaths(threadTriples);
            
            // Add duplicates to paths
            threadPath = reco[i]->addDuplicates(threadPath);
            
            path[i].insert(path[i].end(),threadPath.begin(),threadPath.end());
            
            if (!eval) {
                pairs[i].insert(pairs[i].end(),threadPairs.begin(),threadPairs.end());
                triples[i].insert(triples[i].end(),threadTriples.begin(),threadTriples.end());
            }
            
        }
        
        // Assign tracks
        path[i] = reco[i]->findAssignment(path[i],mod,1); // Trash=1 prunes junk tracks in the output (high numbers)
        
    }
    PARALLEL_FOR_END();
    
    memoryFootprint();
    
    // Accumulate the partial results
    
    cout << "Accumulate the partial results" << endl;
    
    paths.clear();
    for (int i=0;i<n;i++) paths.insert(paths.end(),path[i].begin(),path[i].end());
    
    long npaths = paths.size();
    cout <<  npaths << " paths" << endl;
    
    if (!eval) {
        vector<triple> pairstotal;
        for (int i=0;i<n;i++) pairstotal.insert(pairstotal.end(),pairs[i].begin(),pairs[i].end());
        scorePairs(pairstotal);
        
        vector<triple> triplestotal;
        for (int i=0;i<n;i++) triplestotal.insert(triplestotal.end(),triples[i].begin(),triples[i].end());
        scorePairs(triplestotal);
        scoreTriples(triplestotal);
        
        accuracy = scorePaths4(paths);
        cscore = timer(accuracy);
    }
    
    // Find hit assignments to paths
    
    for (int i=0;i<(int)hits.size();i++) assignment[i] = 0;
    if (paths.size()==0) return assignment;
    
    cout << "Find track assignment" << endl;
    paths = reco[0]->findAssignment(paths,mod,1); // Trash=1 prunes junk tracks in the output (high numbers)
    
    if (!eval) {
        accuracy = scorePaths4(paths);
        cscore = timer(accuracy);
    }
    
    map<int, int> assignment_part;
    for (int k = 1; k < (int)paths.size(); k++)
        for (int j : paths[k])
            if (j > 0 && j < (int)hits.size())
                assignment_part[j] = k;
    
    
    for (pair<int, int> p : assignment_part)
        assignment[p.first] = p.second;
    
    map<int, int> map_assignment;
    for (int i = 1; i < (int)hits.size(); i++)
        if (assignment[i]) map_assignment[i] = assignment[i];
    
    if (!eval) {
        accuracy = scoreAssignment(map_assignment);
        scoreAssignment2(map_assignment);
        cscore = timer(accuracy);
    }
    
    cout << "Write out solution to submission" << filenum << ".csv" << endl;
    writeSubmission(map_assignment,workpath,filenum);
    
    
    // Clean up the reconstruction objects
    
    for (int i=0;i<n;i++) delete reco[i];
    
    memoryFootprint();
    
    return assignment;
}


//Convert "dir" to polar (really cylindrical coordinates) coordinates, suffix p  (as in "refp") usually means polar coodinates throughout the code
point topolar(const point&dir, const point&ref, const point&refp) {
    return point(ref.x*dir.x+ref.y*dir.y, ref.x*dir.y-ref.y*dir.x, dir.z*refp.x);
}


// function to print memory consumption
void Tracker::memoryFootprint()
{
#ifdef MEMORY
    double vm, rss;
    process_mem_usage(vm, rss);
    cout << "VM: " << vm << "; RSS: " << rss << " kb" << endl;
#endif
}


// Function to determine slice number for a hit
int Tracker::slice(int ai)
{
#ifdef OCTANTSLICES
    point &p = hits[ai];
    double x = p.x;
    double y = p.y;
    double z = p.z;
    if (x >= 0 && y >= 0 && z >= 0)
        return 0;
    else if (x < 0 && y >= 0 && z >= 0)
        return 1;
    else if (x < 0 && y < 0 && z >= 0)
        return 2;
    else if (x >= 0 && y < 0 && z >= 0)
        return 3;
    else if (x >= 0 && y >= 0 && z < 0)
        return 4;
    else if (x < 0 && y >= 0 && z < 0)
        return 5;
    else if (x < 0 && y < 0 && z < 0)
        return 6;
    else
        return 7;
#else
    // spherical slices
    int index = NSLICES*graphHash(ai).first/Reconstruction::phiResolution(); // phi slicing
    //int index = NSLICES*graphHash(ai).second/THETARESOLUTION; // theta slicing
    if (index>NSLICES-1) index = NSLICES-1;
    if (index<0) index = 0;
    return index;
#endif
}


// hash function to define a volume element for a hit
long long Tracker::voxel(int hitid) {
    long l=metai[hitid];
    long m=meta[hitid].z;
    auto h = graphHash(hitid);
    long long index = (((long)h.first)<<32) | (((long)h.second)<<24) | (l<<16) | m;
    return index;
}


// hash function to define a graph for a hit
std::pair<int,int> Tracker::graphHash(int hitid) {
    point pol = polar[hitid];
    double rt = pol.x;
    double ph = pol.y;
    double z0 = pol.z;
    double th = asinh(z0/rt); // slope (-4...+4)
    int i1 = (int) (Reconstruction::phiResolution()*0.15*(M_PI+ph));
    int i2 = (int) (Reconstruction::thetaResolution()*0.1*(5-th));
    return make_pair(i1,i2);
}


// read tile graphs from a directory
void Tracker::readTileGraphs(const char *directory) {
    
    DIR *dir;
    struct dirent *ent;
    string dirname(directory);
    
    // read all files and fill graph
    if ((dir = opendir(directory)) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            string file(ent->d_name);
            if (file=="." || file=="..") continue;
            int ph,th;
            sscanf(ent->d_name,"%02d.%02d",&ph,&th);
            pair<int,int> index = make_pair(ph,th);
            string filename = dirname + "/" + file;
            readGraph(filename.c_str(), tgraph[index]);
        }
        closedir (dir);
    }
    
    cout << "Reading " << tgraph.size() << " tracking graphs" << endl;
}


// generate PHIRESOLUTION.THETARESOLUTION tile graphs and write them to a directory
void Tracker::generateTileGraphs(const char *directory) {
    
    for (auto &p : particles) {
        
        if (p.hit.size()<3) continue;
        
        int ai = p.hit[0];
        pair<int,int> index = graphHash(ai);
        
        long indexa,indexb;
        int n = (int) p.hits-1;
        
        for (int i=0;i<n;i++) {
            int ai = p.hit[i];
            int bi = p.hit[i+1];
            int la = metai[ai];
            int lb = metai[bi];
            if (la==lb && i<(p.hits-2)) { // avoid double hit (same layer)
                bi = p.hit[i+2];
                lb = metai[bi];
            }
            indexa = voxel(ai);
            indexb = voxel(bi);
            if (la>lb) swap(indexa,indexb);
            point &v = start_pos[truth_part[ai]];
            float vz = v.z;
            tgraph[index].add(indexa,indexb,vz);
        }
    }
    
    cout << "Writing " << tgraph.size() << " tracking graphs" << endl;
    for (auto &t : tgraph) {
        pair<int,int> index = t.first;
        char file[1000];
        sprintf(file, "%s/%02d.%02d", directory, index.first,index.second);
        writeGraph(file, t.second);
    }
    
    
}


// initialize the ANNs
void Tracker::initNeuralNetworks()
{
    string directory(string(workPath)+"/"+XMLPDIR+"/");
    
    cout << "initNeuralNetworks: Reading networks from " << directory << endl;
    cout << NETFILE1 << " " << NETFILE2 << " " << NETFILE3 << endl;
    
#ifdef USETMVA

    // TMVA setup
#ifdef FOLDEDINPUT1
   vector<string> inputVar1 = {"rz1","abs(abs(phi1)-1.57079632679)","abs(z1)","rz2","abs(abs(phi2)-1.57079632679)","abs(z2)","f0","f1"};
#else
    vector<string> inputVar1 = {"rz1","phi1","z1","rz2","phi2","z2","f0","f1"};
#endif
    
#ifdef FOLDEDINPUT2
   vector<string> inputVar2 = {"rz1","abs(abs(phi1)-1.57079632679)","abs(z1)","rz2","abs(abs(phi2)-1.57079632679)","abs(z2)","f0","f1","log(score)"};
#else
    vector<string> inputVar2 = {"rz1","phi1","z1","rz2","phi2","z2","f0","f1","log(score)"};
#endif
    
#ifdef FOLDEDINPUT3
    vector<string> inputVar3 = {"rz1","abs(abs(phi1)-1.57079632679)","abs(z1)","rz2","abs(abs(phi2)-1.57079632679)","abs(z2)","rz3","abs(abs(phi3)-1.57079632679)","abs(z3)","log(score)"};
#else
     vector<string> inputVar3 = {"rz1","phi1","z1","rz2","phi2","z2","rz3","phi3","z3","log(score)"};
#endif

#ifdef TMVAREADER
    string file1 = string(directory)+"/"+TMVAFILE1;
    string file2 = string(directory)+"/"+TMVAFILE2;
    string file3 = string(directory)+"/"+TMVAFILE3;
    
    for (int i=0;i<NSLICES;i++) {
        if (pgraph[i].getReader1()==NULL) pgraph[i].setReader1(new TMVA::Reader( "!Color:!Silent" ));
        if (pgraph[i].getReader2()==NULL) pgraph[i].setReader2(new TMVA::Reader( "!Color:!Silent" ));
        if (pgraph[i].getReader3()==NULL) pgraph[i].setReader3(new TMVA::Reader( "!Color:!Silent" ));
        
        TMVA::Reader *reader1 = pgraph[i].getReader1();
        float *x1 = pgraph[i].getX1();
        for (int j=0;j<8;j++) {
            reader1->AddVariable(inputVar1[j].c_str(),&x1[j]);
        }
        reader1->BookMVA("MLP method",file1.c_str());
        
        TMVA::Reader *reader2 = pgraph[i].getReader2();
        float *x2 = pgraph[i].getX2();
        for (int j=0;j<9;j++) {
            reader2->AddVariable(inputVar2[j].c_str(),&x2[j]);
        }
        reader2->BookMVA("MLP method",file2.c_str());
        
        TMVA::Reader *reader3 = pgraph[i].getReader3();
        float *x3 = pgraph[i].getX3();
        for (int j=0;j<10;j++) {
            reader3->AddVariable(inputVar3[j].c_str(),&x3[j]);
        }
        reader3->BookMVA("MLP method",file3.c_str());
    }

    for (auto &t : tgraph) {
        auto &g = t.second;
        if (g.getReader1()==NULL) g.setReader1(new TMVA::Reader( "!Color:!Silent" ));
        if (g.getReader2()==NULL) g.setReader2(new TMVA::Reader( "!Color:!Silent" ));
        if (g.getReader3()==NULL) g.setReader3(new TMVA::Reader( "!Color:!Silent" ));

        TMVA::Reader *reader1 = g.getReader1();
        float *x1 = g.getX1();
        for (int j=0;j<8;j++) {
            reader1->AddVariable(inputVar1[j].c_str(),&x1[j]);
        }
        reader1->BookMVA("MLP method",file1.c_str());
        
        TMVA::Reader *reader2 = g.getReader2();
        float *x2 = g.getX2();
        for (int j=0;j<9;j++) {
            reader2->AddVariable(inputVar2[j].c_str(),&x2[j]);
        }
        reader2->BookMVA("MLP method",file2.c_str());
        
        TMVA::Reader *reader3 = g.getReader3();
        float *x3 = g.getX3();
        for (int j=0;j<10;j++) {
            reader3->AddVariable(inputVar3[j].c_str(),&x3[j]);
        }
        reader3->BookMVA("MLP method",file3.c_str());
    }

#endif
    
    for (int i=0;i<NSLICES;i++) {
        if (pgraph[i].getNet1()==NULL) pgraph[i].setNet1(new ReadMLP1(inputVar1));
        if (pgraph[i].getNet2()==NULL) pgraph[i].setNet2(new ReadMLP2(inputVar2));
        if (pgraph[i].getNet3()==NULL) pgraph[i].setNet3(new ReadMLP3(inputVar3));
    }

    for (auto &t : tgraph) {
        auto &g = t.second;
        if (g.getNet1()==NULL) g.setNet1(new ReadMLP1(inputVar1));
        if (g.getNet2()==NULL) g.setNet2(new ReadMLP2(inputVar2));
        if (g.getNet3()==NULL) g.setNet3(new ReadMLP3(inputVar3));
    }    

#else

    for (int i=0;i<NSLICES;i++) {
        string name1(NETFILE1);
        string name2(NETFILE2);
        string name3(NETFILE3);
        if (pgraph[i].net1()==NULL) pgraph[i].setNet1(new XMLP((directory+name1).c_str()));
        if (pgraph[i].net2()==NULL) pgraph[i].setNet2(new XMLP((directory+name2).c_str()));
        if (pgraph[i].net3()==NULL) pgraph[i].setNet3(new XMLP((directory+name3).c_str()));
    }
    
    for (auto &t : tgraph) {
        string netfile1(string(workPath)+"/"+XMLPDIR+"/"+NETFILE1);
        string netfile2(string(workPath)+"/"+XMLPDIR+"/"+NETFILE2);
        string netfile3(string(workPath)+"/"+XMLPDIR+"/"+NETFILE3);
        auto &g = t.second;
        if (g.net1()==NULL) g.setNet1(new XMLP(netfile1.c_str()));
        if (g.net2()==NULL) g.setNet2(new XMLP(netfile2.c_str()));
        if (g.net3()==NULL) g.setNet3(new XMLP(netfile3.c_str()));
    }

#endif
}


void Tracker::initPolarModule() {
    for (int i = 0; i < LAYERS; i++)
        if (mod[i] == NULL) mod[i] = new PolarModule(i,this);
}


// generate a task list to be executed by the threads
void Tracker::initTasks()
{
#ifdef SLICES
    // divide the threads into hemispheres
    const int order16[16] = {0,8,1,9,2,10,3,11,4,12,5,13,6,14,7,15};
    const int order8[8]   = {0,4,1,5,2,6,3,7};
    const int order6[6]   = {0,3,1,4,2,5};
    const int order4[4]   = {0,2,1,3};
    cout << "Perform tile tasks" << endl;
    for (int i=0;i<NSLICES;i++) {
        int n = i;
        if (NSLICES==16)
            n = order16[i];
        if (NSLICES==8)
            n = order8[i];
        else if (NSLICES==6)
            n = order6[i];
        else if (NSLICES==4)
            n = order4[i];
        
        tasks.push(make_pair(n,&pgraph[n]));
    }
#else
    cout << "Perform " << tgraph.size() << "graph tasks" << endl;
    int task(0);
    for (auto &t : tgraph) {
        Graph<long long> *g = &tgraph[t.first];
        tasks[task%NTHREADS].push(make_pair(task++,g));
    }
#endif
}


// Instantiate reconstruction objects for parallel processing
Reconstruction *Tracker::initRecoObjects(int n)
{
    static const double t1offsets[15] = {0,0,0,0,0,0,0,0,-0.01,0.01,0.02,0.02,0,0,0};
    static const double t2offsets[15] = {0,0,0,0,0,0,0,0,-0.01,0.01,0.02,0.02,0,0,0};
    static const double t3offsets[15] = {0,0,0,0,0,0,0,0,-0.2,-0.1,-0.1,0,0,0,0};
    threshold1 += t1offsets[(getNumberHits()/10000)%15]; // tune the nerual network cut wrt event size
    threshold2 += t2offsets[(getNumberHits()/10000)%15]; // tune the nerual network cut wrt event size
    threshold3 += t3offsets[(getNumberHits()/10000)%15]; // tune the nerual network cut wrt event size
    cout << "initRecoObjects: Instantiate " << n << " reconstruction objects for parallel processing" << endl;
    cout << "THRESHOLD1: " << threshold1 << "  THRESHOLD2: " << threshold2 << "  THRESHOLD3: " << threshold3 << endl;
    for (int task=0;task<n;task++) {
        reco[task] = new Reconstruction(filenum,workPath.c_str(),tgraph,hits,polar,meta,metai,metaz,layer,disc_z,hit_cells,hit_dir);
        reco[task]->setThreshold1(threshold1);
        reco[task]->setThreshold2(threshold2);
        reco[task]->setThreshold3(threshold3);
    }
    
    return *reco;
}


//does hits a and b correspond to the same particle?
int Tracker::samepart(int a, int b) {
    long long aa = truth_part[a];
    long long bb = truth_part[b];
    return aa==bb;//aa == bb && aa;
}


void Tracker::sortZ(vector<int> &hits)
{
    std::sort(hits.begin(), hits.end(), [this](int a, int b) {return z_cmp(a, b); });
}


void Tracker::sortR(vector<int> &hits)
{
    std::sort(hits.begin(), hits.end(), [this](int a, int b) {return r_cmp(a, b); });
}


void Tracker::sortDist(vector<int> &hits)
{
    std::sort(hits.begin(), hits.end(), [this](int a, int b) {return dist_cmp(a, b); });
}


inline bool Tracker::z_cmp(const int &a, const int &b) {
    return hits[a].z < hits[b].z;
}


inline bool Tracker::r_cmp(const int &a, const int &b) {
    return hits[a].x*hits[a].x+hits[a].y*hits[a].y < hits[b].x*hits[b].x+hits[b].y*hits[b].y;
}


inline bool Tracker::dist_cmp(const int &a, const int &b) {
    return hits[a].x*hits[a].x+hits[a].y*hits[a].y+hits[a].z*hits[a].z < hits[b].x*hits[b].x+hits[b].y*hits[b].y+hits[b].z*hits[b].z;
}


// Prepare the hits in modules
long Tracker::initGraphData()
{
    long nhits = hits.size();
    for (int i=0;i<NSLICES;i++) pgraph[i].clear();
    for (auto &t : tgraph) t.second.clear();
    
    long n(0);
    for (int i=1; i<nhits; i++) {
        if (assignment[i] != 0) continue; // Skip hits
        int l = metai[i];
        if (l<0 || l>=LAYERS) continue;
        long long index = voxel(i);
        int ti = slice(i);
        pgraph[ti].data(index).push_back(i);
        pgraph[ti].hash(l).insert(index);
        pair<int,int> cindex = graphHash(i);
        tgraph[cindex].data(index).push_back(i);
        tgraph[cindex].hash(l).insert(index);
        n++;
        //if (verbose) cout << i << " " << l << " " << m << " " << index << endl;
    }
    
    return n;
}

#ifdef MEMORY
#include <iostream>
#include <fstream>
#include <unistd.h>

// check memory usage
void process_mem_usage(double& vm_usage, double& resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;
    
    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
        >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
        >> ignore >> ignore >> vsize >> rss;
    }
    
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}
#else
void process_mem_usage(double& vm_usage, double& resident_set)
{}
#endif


#include <ctime>
#define MAX_TIME_PER_EVENT 600
#define ACCURACY_INF 0.5

// check execution time
double Tracker::timer(double accuracy_mean)
{
    static clock_t c_start = clock();
    double score = 0.0;
    
    if (accuracy_mean < 0) {
        c_start = clock();
    }
    else {
        clock_t c_end = clock();
        double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
        
        double speed = log(1.0 + ( MAX_TIME_PER_EVENT / (time_elapsed_ms*0.001) ));
        /*if (accuracy_mean>ACCURACY_INF)*/ score = sqrt(speed * pow((accuracy_mean - ACCURACY_INF),2));
        cout << "CodaLab score: " << score << endl;
        cout << "CPU time used: " << time_elapsed_ms << " ms" << endl;
    }
    
    return score;
}



// Print an integer hit id track vector
void Tracker::print(vector<int> const &input)
{
    for (unsigned int i = 0; i < input.size(); i++) {
        cout << input.at(i) << ' ';
    }
    cout << endl;
}


// Print an integer short hit id track vector
void Tracker::printshort(vector<int> const &input)
{
    for (unsigned int i = 0; i < input.size(); i++) {
        cout << shortid(input.at(i)) << ' ';
    }
    cout << endl;
}


// write the result file
void Tracker::writeSubmission(map<int, int>&assignment,string path=WORKPATH, int filenum=21100) {
    char filename[1000];
    sprintf(filename, "%s/submission%d.csv",path.c_str(),filenum);
    FILE*fp = fopen(filename, "w");
    if (!fp) {
        cout << "Could not open " << filename << " for writing submission" << endl;
        return;
    }
    fprintf(fp, "event_id,hit_id,track_id\n");
    
    for (int i = 1; i < (int)hits.size(); i++) {
        int a = 0;
        if (assignment.count(i)) a = assignment[i];
        fprintf(fp, "%d,%d,%d\n", filenum, i, a);
    }
    fclose(fp);
}


// Write graph data to file
void Tracker::writeGraph(const char *file, Graph<pair<int,int> > &g) {
    string filename(file);
    cout << "Writing graph data to " << filename << endl;
    ofstream path(filename);
    if (path.is_open()) {
        path << g;
        path.close();
    }
    else
        cout << "Unable to open file " << filename << endl;
}


// Write graph data to file
void Tracker::writeGraph(const char *file, Graph<long long> &g) {
    string filename(file);
    cout << "Writing graph data to " << filename << endl;
    ofstream path(filename);
    if (path.is_open()) {
        path << g;
        path.close();
    }
    else
        cout << "Unable to open file " << filename << endl;
}



// Read path data from file
void Tracker::readGraph(const char *file, Graph<pair<int,int> > &g) {
    string filename(file);
    //cout << "Reading graph data from " << filename << endl;
    ifstream input;
    input.open(filename);
    if (!input.is_open()) {
        cerr << "couldn't open " << filename << endl;
    }
    else {
        input >> g;
        input.close();
    }
}


// Read path data from file
void Tracker::readGraph(const char *file, Graph<long long> &g) {
    string filename(file);
    //cout << "Reading graph data from " << filename << endl;
    ifstream input;
    input.open(filename);
    if (!input.is_open()) {
        cerr << "couldn't open " << filename << endl;
    }
    else {
        input >> g;
        input.close();
    }
}


// Generate graphs to represent the track hits in a detector slice
void Tracker::generateSliceGraphs() {
    int tracknumber = 1;
    for (auto &track : truth_tracks) {
        vector<int> &t = track.second;
        if (t.size()==0) continue;
        if (debug) cout << "Track " << tracknumber << ", size " << t.size() << " {" << shortid(t[0]) << "," << shortid(t[1]) << "}" << endl;
        long previndex = voxel(t[0]);
        long n = 0;
        
        int ti = slice(t[0]);
        
        for (auto &hit_id : t) {
            long index = voxel(hit_id);
            point &v = start_pos[truth_part[hit_id]];
            float vz = v.z;
            // Add the hit pair to the paths graph
            pgraph[ti].add(previndex,index,vz); // forward direction, note z vertex
            if (debug) cout << shortid(hit_id) << "={" << previndex << ","  << index << "},";
            previndex = index;
            n++;
        }
        
        tracknumber++;
    }
}


// read graphs to represent the track hits in a detector slice
void Tracker::readSliceGraphs()
{
    for (int i=0;i<NSLICES;i++) {
        stringstream ss;
        ss << workPath << "/" << SLICEDIR << "/paths"  << i << ".csv";
        readGraph(ss.str().c_str(),pgraph[i]);
    }
}


// write graphs to represent the track hits in a detector slice
void Tracker::writeSliceGraphs()
{
    for (int i=0;i<NSLICES;i++) {
        stringstream ss;
        ss << workPath << "/" << SLICEDIR << "/paths" << i << ".csv";
        writeGraph(ss.str().c_str(),pgraph[i]);
    }
}


// The following code has been adopted from Johan Sokrates Wind's award winning trackml Kaggle contribution

//Contains functions to estimate score at different steps in the pipeline (assuming the remaining steps are done optimally)

double Tracker::scorePairs(vector<triple>&pairs) {
    set<long long> found;
    
    for (auto &p : pairs) {
        if (samepart(get<0>(p),get<1>(p))) {
            found.insert(truth_part[get<0>(p)]);
        }
    }
    double score = 0;
    for (long long p : found) {
        score += part_weight[p];
    }
    
    cout << "Score: " << score << " from " << pairs.size() << " pairs" << endl;
    
    return score;
}


double Tracker::scoreTriples(vector<triple>&triples) {
    set<long long> found;
    
    for (auto &t : triples) {
        if (samepart(get<0>(t),get<1>(t)) && samepart(get<0>(t),get<2>(t))) {
            found.insert(truth_part[get<0>(t)]);
        }
    }
    double score = 0;
    for (long long p : found) {
        score += part_weight[p];
    }
    cout << "Score: " << score << " from " << triples.size() << " triples" << endl;
    
    return score;
}


double Tracker::scorePaths4(vector<vector<int> > &paths) {
    int total_length = 0;
    map<long long, double> score_part;
    
    for (auto &path : paths) {
        total_length += path.size();
        set<long long> done;
        //int minmeta = 1e9;
        //for (int i : path) minmeta = min(minmeta, metai[i]);
        for (int i : path) {
            if (i<0 || i>(int)hits.size()) continue;
            long long part = truth_part[i];
            if (done.count(part)) continue;
            done.insert(part);
            double w = 0;
            for (int i : truth_tracks[part]) {
                for (int j : path) {
                    if (j<0 || j>(int)hits.size()) continue;
                    if (truth_part[j] == part && metai[j] == metai[i] && point::dist(hits[i]-hits[j]) < 10) {
                        w += truth_weight[i];
                        break;
                    }
                }
            }
            
            score_part[part] = max(score_part[part], w);
        }
    }
    
    double score = 0;
    for (auto &p : score_part) if (p.first) score += p.second;
    cout << "Score: " << score << " from " << paths.size() << " paths with " << total_length << " hits" << endl;
    
    return score;
}


double Tracker::scoreAssignment(map<int, int>&assignment) {
    map<int, int> track_length;
    
    for (auto &p : assignment)
        track_length[p.second]++;
    
    double score = 0, score2 = 0;
    for (auto &it : truth_tracks) {
        map<int, int> c;
        for (int i : it.second)
            if (assignment.count(i))
                c[assignment[i]]++;
        int pick = -1, maxlen = -1;
        for (auto &p : c)
            if (p.second*2 > maxlen && p.second*2 > track_length[p.first]) {
                pick = p.first;
                maxlen = p.second*2;
            }
        if (pick == -1) continue;
        
        for (int i : it.second)
            if (assignment[i] == pick) {
                if (maxlen > (int)it.second.size()) score += truth_weight[i];
                else score2 += truth_weight[i];
            }
    }
    //"Final score: " should be the same as the official score (except for blacklisted electrons)
    cout << "Final score:  " << score << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
    cout << "Short score:  " << score+score2 << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
    
    return score;
}


void Tracker::scoreAssignment2(map<int,int> &assignment) {
    map<int, int> track_length;
    
    for (auto &p : assignment)
        track_length[p.second]++;
    
    double score = 0, score2 = 0, score3 = 0, score4 = 0;
    for (auto &it : truth_tracks) {
        map<int, int> c;
        for (int i : it.second) {
            if (assignment.count(i))
                c[assignment[i]]++;
        }
        
        int pick = -1, maxlen = -1;
        for (auto &p : c)
            if (p.second*2 > maxlen && p.second*2 > track_length[p.first]) {
                pick = p.first;
                maxlen = p.second*2;
            }
        if (pick == -1) continue;
        map<int, set<int> > modules;
        for (int i : it.second)
            if (assignment[i] == pick)
                modules[metai[i]].insert(i);
        
        int len = 0, len2 = 0;
        double s = 0, s2 = 0;
        for (int i : it.second) {
            if (modules.count(metai[i])) {
                int found = 3;
                for (int j : modules[metai[i]]) {
                    if (point::dist(hits[i]-hits[j]) < 10) found = min(found, 1);
                    found = min(found, 2);
                }
                if (found == 1) {
                    s += truth_weight[i];
                    len++;
                }
                else if (found == 2) {
                    s2 += truth_weight[i];
                    len2++;
                }
            }
        }
        len2 += len;
        if (len*2 > (int)it.second.size()) score += s;
        else score3 += s;
        if (len2*2 > (int)it.second.size()) score2 += s2;
        else score4 += s2;
    }
    score2 += score;
    score4 += score3;
    score3 += score;
    score4 += score2;
    cout << "No-dup score: " << score << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
    cout << "Far score:    " << score2 << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
    cout << "Near all:     " << score3 << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
    cout << "Far all:      " << score4 << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
}


//Functions to find good pairs, and evaluate a triple, both using logistic regression

//Decides which pairs to fit logistic model to
int Tracker::good_pair(int a, int b) {
    if (a==b) return 0;
    if (truth_part[a]==0 && truth_part[b]==0) return 0; // NA hits
    if (!samepart(a, b)) return 0;
    int ai = metai[a];
    int bi = metai[b];
    if (ai==bi) return 1; // ai==bi, double hits
    point s = start_pos[truth_part[a]];
    if (s.x*s.x+s.y*s.y < (VERTEXCUTXY*VERTEXCUTXY) && abs(s.z)<VERTEXCUTZ) return 2; // inside cylinder / vertex region (primary tracks)
    return -1;
}



//Contains functions to read the input data, and global variables to hold the data


set<long long> blacklist;
void Tracker::readBlacklist() {
    blacklist.clear();
    char file[1000];
    sprintf(file, "%s/event%09d-blacklist_particles.csv", dataPath.c_str(), filenum);
    cout << "Reading " << file << endl;
    FILE*fp = fopen(file, "r");
    if (!fp) { printf("couldn't open blacklist %s\n", file); return; }
    char tmpstr[1000];
    int tmp = fscanf(fp, "%s", tmpstr);
    cout << tmp << ":" << tmpstr << endl;
    long long particle_id;
    while (fscanf(fp, "%lld", &particle_id) == 1) {
        blacklist.insert(particle_id);
    }
    fclose(fp);
}


set<long long> whitelist;
void Tracker::readWhitelist() {
    whitelist.clear();
    char file[1000];
    sprintf(file, "%s/event%09d-whitelist_particles.csv", dataPath.c_str(), filenum);
    cout << "Reading " << file << endl;
    FILE*fp = fopen(file, "r");
    if (!fp) { printf("couldn't open whitelist %s\n", file); return; }
    char tmpstr[1000];
    int tmp = fscanf(fp, "%s", tmpstr);
    cout << tmp << ":" << tmpstr << endl;
    long long particle_id;
    while (fscanf(fp, "%lld", &particle_id) == 1) {
        whitelist.insert(particle_id);
    }
    fclose(fp);
}


void Tracker::readTruth() {
    const double Bfield = 1673.;
    char file[1000];
    sprintf(file, "%s/event%09d-truth.csv", dataPath.c_str(), filenum);
    cout << "Reading " << file << endl;
    FILE*fp = fopen(file, "r");
    if (!fp) { printf("couldn't open %s\n",file); return; }
    char tmpstr[1000];
    int tmp = fscanf(fp, "%s", tmpstr);
    cout << tmp << ":" << tmpstr << endl;
    
    truth_tracks.clear();
    part_weight.clear();
    
    while (1) {
        int hit_id;
        long long particle_id;
        double tx, ty, tz, tpx, tpy, tpz, weight;
        if (fscanf(fp, "%d,%lld,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &hit_id, &particle_id, &tx, &ty, &tz, &tpx, &tpy, &tpz, &weight) != 9) break;
        if (!particle_id || blacklist.count(particle_id)) continue;
        if (whitelist.size()>0 && !whitelist.count(particle_id)) continue; // Accept only selected particles
        if (whitelist.count(particle_id)) whitelist.insert(hit_id);
        
        truth_tracks[particle_id].push_back(hit_id);
        truth_pos[hit_id] = point(tx, ty, tz);
        truth_mom[hit_id] = point(tpx, tpy, tpz)*Bfield;
        truth_weight[hit_id] = weight;
        part_weight[particle_id] += weight;
        truth_part[hit_id] = particle_id;
        
        map<long long, int>::iterator it = partIDmap.find(particle_id);
        if( it==partIDmap.end() ){
            cout<<"Particle ID not found in map!!!"<<endl;
            cout<<"ID= "<<hit_id<<" hit "<<hit_id<<" iterator at ID "<<it->first<<endl;
            exit(0);
        }
        
        int newID = it->second;
        if( newID < 0 || newID>= (int)particles.size() ){
            cout<<"Mapped particle ID is wrong!!!"<<endl;
            cout<<"ID= "<<hit_id<<" new ID "<<newID<<endl;
            exit(0);
        }
        
        Particle &p = particles[newID];
        p.hit.push_back(hit_id);
    }
    
    fclose(fp);
    
    int n = 0;
    for (auto &t : truth_tracks) {
        for (auto hit_id : t.second) {
            hitIDmap[hit_id] = n;  // Short hit id
            n++; // next hit
        }
    }
    
    cout << truth_tracks.size() << " particles with truth" << endl;
}


void Tracker::readStarts() {
    const double Bfield = 1673.;
    char file[1000];
    sprintf(file, "%s/event%09d-particles.csv", dataPath.c_str(), filenum);
    cout << "Reading " << file << endl;
    FILE*fp = fopen(file, "r");
    if (!fp) { printf("couldn't open %s\n",file); return; }
    char tmpstr[1000];
    int tmp = fscanf(fp, "%s", tmpstr);
    cout << tmp << ":" << tmpstr << endl;
    while (1) {
        long long id;
        int type;
        point p, m;
        int hits;
        int q;
        //if (fscanf(fp, "%lld,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d", &id, &p.x, &p.y, &p.z, &m.x, &m.y, &m.z, &q, &hits) == -1) break;
        if (fscanf(fp, "%lld,%d,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d", &id, &type, &p.x, &p.y, &p.z, &m.x, &m.y, &m.z, &q, &hits) == -1) break;
        if (blacklist.count(id)) continue;
        if (whitelist.size()>0 && !whitelist.count(id)) continue; // Accept only selected particles
        start_pos[id] = p;
        start_mom[id] = m*Bfield;
        part_q[id] = -q;
        part_hits[id] = hits;
        
        Particle part;
        partIDmap[ (long long) id ] = (int)particles.size();
        part.id = id;
        part.type = type;
        part.x = p.x;
        part.y = p.y;
        part.z = p.z;
        part.r = 0;
        part.px = m.x;
        part.py = m.y;
        part.pz = m.z;
        part.q = q;
        part.hits = hits;
        particles.push_back(part);
    }
    fclose(fp);
    cout << start_pos.size() << " particles" << endl;
}


void Tracker::sortT(vector<int> &hits)
{
    std::sort(hits.begin(), hits.end(), [this](int a, int b) {return track_cmp(a, b); });
}


inline
bool Tracker::track_cmp(const int a,const int b) {
    point&ma = truth_mom[a];
    point&mb = truth_mom[b];
    double va = ma*ma, vb = mb*mb;
    if (fabs((va-vb)/va) > 1e-5) return va > vb;
    return (truth_pos[b]-truth_pos[a])*ma > 0;
}


//Sort the hits in each track chronologically
void Tracker::sortTracks() {
    int fails = 0, goods = 0;
    for (auto &p : truth_tracks) {
        vector<int> v;
        for (int hit_id : p.second) {
            v.push_back(hit_id);
        }
        sortT(v);
        //sort(v.begin(), v.end(), track_cmp);
        for (int i = 0; i < (int)v.size(); i++)
            p.second[i] = v[i];
        int bad = 0;
        for (int i = 2; i < (int)v.size(); i++) {
            point&a = truth_pos[p.second[i-2]];
            point&b = truth_pos[p.second[i-1]];
            point&c = truth_pos[p.second[i]];
            if ((c.z-b.z)*(b.z-a.z) < 0) {
                fails++;
                bad++;
            }
            else goods++;
        }
    }
}


void Tracker::initOrder() {
    int handpicked[LAYERS] = {};
    int c = 0;
    for (int i = 0; i < 4; i++) handpicked[c++] = 7+i;
    for (int i = 0; i < 7; i++) handpicked[c++] = 7-1-i;
    for (int i = 0; i < 7; i++) handpicked[c++] = 11+i;
    
    for (int i = 0; i < 4; i++) handpicked[c++] = 24+i;
    for (int i = 0; i < 2; i++) handpicked[c++] = 40+i;
    for (int i = 0; i < 6; i++) {
        handpicked[c++] = 24-1-i;
        handpicked[c++] = 40-1-i;
    }
    
    for (int i = 0; i < 6; i++) {
        handpicked[c++] = 28+i;
        handpicked[c++] = 42+i;
    }
    for (int i = 0; i < LAYERS; i++) {
        topo[i] = handpicked[i];
        itopo[topo[i]] = i;
    }
}


//init layer geometries
void Tracker::initLayers() {
    double avgz1[2][7] = {{-1500,-1300,-1100,-960,-820,-700,-600},
        { 600, 700, 820, 960, 1100, 1300, 1500}};
    for (int k = 0; k < 2; k++)
        for (int i = 0; i < 7; i++) {
            layer[k*11+i].minr = 30;
            layer[k*11+i].maxr = 176.5;
            layer[k*11+i].avgz = avgz1[k][i];
            layer[k*11+i].type = Disc;
        }
    double avgz2[2][6] = {{-2950,-2550,-2150,-1800,-1500,-1220},
        { 1220, 1500, 1800, 2150, 2550, 2950}};
    for (int k = 0; k < 2; k++)
        for (int i = 0; i < 6; i++) {
            layer[k*10+i+18].minr = 240;
            layer[k*10+i+18].maxr = 701;
            layer[k*10+i+18].avgz = avgz2[k][i];
            layer[k*10+i+18].type = Disc;
            
            layer[k*8+i+34].minr = 755;
            layer[k*8+i+34].maxr = 1018;
            layer[k*8+i+34].avgz = avgz2[k][i];
            layer[k*8+i+34].type = Disc;
        }
    
    double avgr1[4] = {32.3, 72.1, 116.1, 172.1};
    double avgr2[4] = {260.3, 360.2, 500.2, 660.2};
    double avgr3[2] = {820.2, 1020.2};
    
    for (int i = 0; i < 4; i++) {
        layer[i+7].minz =-491;
        layer[i+7].maxz = 491;
        layer[i+7].avgr = avgr1[i];
        layer[i+7].type = Tube;
    }
    for (int i = 0; i < 4; i++) {
        layer[i+24].minz =-1084;
        layer[i+24].maxz = 1084;
        layer[i+24].avgr = avgr2[i];
        layer[i+24].type = Tube;
    }
    for (int i = 0; i < 2; i++) {
        layer[i+40].minz =-1084;
        layer[i+40].maxz = 1084;
        layer[i+40].avgr = avgr3[i];
        layer[i+40].type = Tube;
    }
    Layer layer2[LAYERS];
    for (int i = 0; i < LAYERS; i++) layer2[i] = layer[i];
    for (int i = 0; i < LAYERS; i++) layer[i] = layer2[topo[i]];
    
    layer[0].var0 = 1e-3;
    layer[1].var0 = 5e-4;
    for (int i = 2; i < 18; i++) layer[i].var0 = 3e-4;
    for (int i = 18; i < 22; i++) layer[i].var0 = 5e-2;
    for (int i = 22; i < LAYERS; i++) layer[i].var0 = i%2 || i == 22 ? 9 : 0.1;
    
    for (int i = 0; i < 4; i++) layer[i].var1 = 0.5;
    for (int i = 4; i < 18; i++) layer[i].var1 = 5;
    for (int i = 18; i < 24; i++) layer[i].var1 = 7;
    for (int i = 24; i < LAYERS; i++) layer[i].var1 = i%2 ? 19 : 11;
}


void Tracker::readHitsInitData() {
    initOrder();
    
    char file[1000];
    if (filenum >= 1000)
        sprintf(file, "%s/event%09d-hits.csv", dataPath.c_str(), filenum);
    else
        sprintf(file, "%s/test/event%09d-hits.csv", dataPath.c_str(), filenum);
    cout << "Reading " << file << endl;
    FILE*fp = fopen(file, "r");
    if (!fp) {
        printf("couldn't open %s\n",file);
        exit(1);
    }
    char tmpstr[1000];
    int tmp = fscanf(fp, "%s", tmpstr);
    cout << tmp << ":" << tmpstr << endl;
    
    // For one indexing
    hits.clear();
    hits.push_back(point(0,0,0));
    polar.clear();
    polar.push_back(point(0,0,0));
    meta.clear();
    meta.push_back(point(0,0,0));
    metai.clear();
    metai.push_back(0);
    
    int layers[9] = {7,4,7,6,4,6,6,2,6};
    int metai_list[9][7];
    int c = 0;
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < layers[i]; j++)
            metai_list[i][j] = c++;
    }
    cout << "Detectors: " << c << endl;
    
    for (int i = 0; i < LAYERS; i++) {
        layer[i].minr = layer[i].minz = 1e9;
        layer[i].maxr = layer[i].maxz =-1e9;
    }
    
    while (1) {
        long long hit_id;
        double tx, ty, tz;
        int volume_id, layer_id, module_id;
        if (fscanf(fp, "%lld,%lf,%lf,%lf,%d,%d,%d", &hit_id, &tx, &ty, &tz, &volume_id, &layer_id, &module_id) == -1) break;
        if (hit_id != (int)hits.size()) cout << "Hit id's not as expected" << endl;
        meta.push_back(point(volume_id, layer_id, module_id));
        
        if (volume_id <= 9)
            volume_id -= 7;
        else if (volume_id <= 14)
            volume_id -= 9;
        else
            volume_id -= 10;
        
        int mi = itopo[metai_list[volume_id][layer_id/2-1]];
        
        hits.push_back(point(tx, ty, tz));
        polar.push_back(point(sqrt(tx*tx+ty*ty), atan2(ty,tx), tz));
        metai.push_back(mi);
        
        double r = sqrt(tx*tx+ty*ty);
        Layer&l = layer[metai_list[volume_id][layer_id/2-1]];
        l.minr = min(l.minr, r);
        l.avgr += r;
        l.maxr = max(l.maxr, r);
        l.minz = min(l.minz, tz);
        l.avgz += tz;
        l.maxz = max(l.maxz, tz);
        l.count++;
        //cerr << tz << ' ' << r << endl;
        
        assignment[hit_id] = 0.0;
    }
    
    fclose(fp);
    cout << hits.size() << " hits" << endl;
    
    initLayers();
    
    for (int hit_id = 1; hit_id < (int)hits.size(); hit_id++) {
        metai_weight[truth_part[hit_id]][metai[hit_id]] += truth_weight[hit_id];
    }
    
    map<double, double> mir[LAYERS], mar[LAYERS];
    for (int i = 1; i < (int)hits.size(); i++) {
        int mi = metai[i];
        if (layer[mi].type != Disc) continue;
        double&mir_ = mir[mi][polar[i].z];
        double&mar_ = mar[mi][polar[i].z];
        if (!mir_) mir_ = 1e9;
        mir_ = min(mir_, polar[i].x);
        mar_ = max(mar_, polar[i].x);
    }
    
    map<double, int> zi[LAYERS];
    for (int mi = 0; mi < LAYERS; mi++) {
        if (layer[mi].type != Disc) continue;
        int k = 0;
        for (auto &p : mir[mi]) {
            double mir_ = mir[mi][p.first]-1e-5;
            double mar_ = mar[mi][p.first]+1e-5;
            z_minr[mi][k] = mir_;
            z_maxr[mi][k] = mar_;
            disc_z[mi][k] = p.first;
            zi[mi][p.first] = k++;
        }
    }
    
    metaz.resize(hits.size());
    metaz[0] = 0;
    for (int i = 1; i < (int)hits.size(); i++) {
        int mi = metai[i];
        metaz[i] = meta[i].z;
        if (layer[mi].type == Disc)
            metaz[i] = zi[mi][hits[i].z];
    }
    
}


// Generate an index to address the detector layers 0..47
int Tracker::getLayer(int volume_id, int layer_id) {
    
    const int itopo[LAYERS] = {10,9,8,7,6,5,4,0,1,2,3,11,12,13,14,15,16,17,34,32,30,28,26,24,18,19,20,21,36,38,40,42,44,46,35,33,31,29,27,25,22,23,37,39,41,43,45,47};
    const int metai_list[9][7] = {{0,1,2,3,4,5,6},{7,8,9,10,-1,-1,-1},{11,12,13,14,15,16,17},{18,19,20,21,22,23,-1},{24,25,26,27,-1,-1,-1},{28,29,30,31,32,33,-1},{34,35,36,37,38,39,-1},{40,41,-1,-1,-1,-1,-1},{42,43,44,45,46,47}};
    
    if (volume_id<7 || volume_id>18) return -1;
    if (volume_id <= 9)
        volume_id -= 7;
    else if (volume_id <= 14)
        volume_id -= 9;
    else
        volume_id -= 10;
    
    layer_id = layer_id/2-1;
    
    int index = itopo[metai_list[volume_id][layer_id]];
    return index;
}


void Tracker::readDetectors(const char *file) {
    cout << "Reading " << file << endl;
    FILE*fp = fopen(file, "r");
    if (!fp) {
        printf("couldn't open %s\n",file);
        exit(1);
    }
    
    char tmpstr[1000];
    int tmp = fscanf(fp, "%s", tmpstr);
    cout << tmp << ":" << tmpstr << endl;
    
    Detector d;
    while (fscanf(fp, "%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &d.volume_id, &d.layer_id, &d.module_id, &d.c.x, &d.c.y, &d.c.z,
                  // &d.rx.x, &d.rx.y, &d.rx.z,
                  // &d.ry.x, &d.ry.y, &d.ry.z,
                  // &d.rz.x, &d.rz.y, &d.rz.z,
                  &d.rx.x, &d.ry.x, &d.rz.x,
                  &d.rx.y, &d.ry.y, &d.rz.y,
                  &d.rx.z, &d.ry.z, &d.rz.z,
                  &d.d, &d.minw, &d.maxw, &d.h, &d.cell_w, &d.cell_h) == 21) {
        //if (d.module_id >= 10000 || d.layer_id >= 1000 || d.volume_id >= 100) cout << "What!?" << endl;
        detectors[d.volume_id*10000000+d.layer_id*10000+d.module_id] = d;
    }
    fclose(fp);
}


void Tracker::readCells() {
    char file[1000];
    if (filenum >= 1000)
        sprintf(file, "%s/event%09d-cells.csv", dataPath.c_str(), filenum);
    else
        sprintf(file, "%s/test/event%09d-cells.csv", dataPath.c_str(), filenum);
    cout << "Reading " << file << endl;
    FILE*fp = fopen(file, "r");
    if (!fp) {
        printf("couldn't open %s\n",file);
        exit(1);
    }
    char tmpstr[1000];
    int tmp = fscanf(fp, "%s", tmpstr);
    cout << tmp << ":" << tmpstr << endl;
    int hit_id, ch0, ch1;
    double value;
    int c0 = 0;
    int c1 = 0;
    for (auto &it : hit_cells) it.clear();
    while (fscanf(fp, "%d,%d,%d,%lf", &hit_id, &ch0, &ch1, &value) == 4) {
        hit_cells[hit_id].push_back(make_pair(make_pair(ch0, ch1), value));
        //cout << hit_id << ' ' << ch0 << ' ' << ch1 << ' ' << value << endl;
        if (ch0>c0) c0 = ch0;
        if (ch1>c1) c1 = ch1;
    }
    fclose(fp);
    cout << "Max:" << c0 << " " << c1 << endl;
}


void Tracker::initCells(int ncells,int* id,int *ch0,int *ch1,double *value) {
    
    for (auto &it : hit_cells) it.clear();
    
    for (int i=0;i<ncells;i++) {
        int hit_id = id[i];
        hit_cells[hit_id].push_back(make_pair(make_pair(ch0[i], ch1[i]), value[i]));
        //cout << i << " " << id[i] << " " << ch0[i] << " " << ch1[i] << " " << value[i] << endl;
    }
}


void Tracker::initHits(int nhits,double *x,double *y,double *z,int *vol,int *lay,int *mod) {
    
    initOrder();
    
    for (int i = 0; i < LAYERS; i++) {
        layer[i].minr = layer[i].minz = 1e9;
        layer[i].maxr = layer[i].maxz =-1e9;
    }
    
    // Copy data
    hits.clear();
    hits.reserve(nhits);
    hits.push_back(point(0,0,0)); // For one indexing
    polar.clear();
    polar.reserve(nhits);
    polar.push_back(point(0,0,0));
    meta.clear();
    meta.reserve(nhits);
    meta.push_back(point(0,0,0));
    metai.clear();
    metai.reserve(nhits);
    metai.push_back(0);
    
    int layers[9] = {7,4,7,6,4,6,6,2,6};
    int metai_list[9][7];
    int c = 0;
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < layers[i]; j++)
            metai_list[i][j] = c++;
    }
    cout << "Detectors: " << c << endl;
    
    for (int i = 0; i < LAYERS; i++) {
        layer[i].minr = layer[i].minz = 1e9;
        layer[i].maxr = layer[i].maxz =-1e9;
    }
    
    for (int i=0;i<nhits;i++) {
        meta.push_back(point(vol[i], lay[i], mod[i]));
        
        if (vol[i] <= 9)
            vol[i] -= 7;
        else if (vol[i] <= 14)
            vol[i] -= 9;
        else
            vol[i] -= 10;
        
        int mi = itopo[metai_list[vol[i]][lay[i]/2-1]];
        
        hits.push_back(point(x[i],y[i],z[i]));
        polar.push_back(point(sqrt(x[i]*x[i]+y[i]*y[i]), atan2(y[i],x[i]), z[i]));
        metai.push_back(mi);
        
        double r = sqrt(x[i]*x[i]+y[i]*y[i]);
        Layer&l = layer[metai_list[vol[i]][lay[i]/2-1]];
        l.minr = min(l.minr, r);
        l.avgr += r;
        l.maxr = max(l.maxr, r);
        l.minz = min(l.minz, z[i]);
        l.avgz += z[i];
        l.maxz = max(l.maxz, z[i]);
        l.count++;
        
        long n = hits.size();
        assignment[n] = 0.0;
        if (whitelist.size()>0 && !whitelist.count(n)) assignment[n] = -1.0;
    }
    
    initLayers();
    
    for (int hit_id = 1; hit_id < (int)hits.size(); hit_id++) {
        metai_weight[truth_part[hit_id]][metai[hit_id]] += truth_weight[hit_id];
    }
    
    map<double, double> mir[LAYERS], mar[LAYERS];
    for (int i = 1; i < (int)hits.size(); i++) {
        int mi = metai[i];
        if (layer[mi].type != Disc) continue;
        double&mir_ = mir[mi][polar[i].z];
        double&mar_ = mar[mi][polar[i].z];
        if (!mir_) mir_ = 1e9;
        mir_ = min(mir_, polar[i].x);
        mar_ = max(mar_, polar[i].x);
    }
    
    map<double, int> zi[LAYERS];
    for (int mi = 0; mi < LAYERS; mi++) {
        if (layer[mi].type != Disc) continue;
        int k = 0;
        for (auto &p : mir[mi]) {
            double mir_ = mir[mi][p.first]-1e-5;
            double mar_ = mar[mi][p.first]+1e-5;
            z_minr[mi][k] = mir_;
            z_maxr[mi][k] = mar_;
            disc_z[mi][k] = p.first;
            zi[mi][p.first] = k++;
        }
    }
    
    metaz.resize(hits.size());
    metaz[0] = 0;
    for (int i = 1; i < (int)hits.size(); i++) {
        int mi = metai[i];
        metaz[i] = meta[i].z;
        if (layer[mi].type == Disc)
            metaz[i] = zi[mi][hits[i].z];
    }
}


void Tracker::readHits() {
    char file[1000];
    if (filenum >= 1000)
        sprintf(file, "%s/event%09d-hits.csv", dataPath.c_str(), filenum);
    else
        sprintf(file, "%s/test/event%09d-hits.csv", dataPath.c_str(), filenum);
    cout << "Reading " << file << endl;
    FILE*fp = fopen(file, "r");
    if (!fp) {
        printf("couldn't open %s\n",file);
        exit(1);
    }
    char tmpstr[1000];
    int tmp = fscanf(fp, "%s", tmpstr);
    cout << tmp << ":" << tmpstr << endl;
    
    double x[MAXDIM],y[MAXDIM],z[MAXDIM];
    int v[MAXDIM],l[MAXDIM],m[MAXDIM];
    
    int n(0);
    while (n>=0) {
        long long hit_id;
        double tx, ty, tz;
        int volume_id, layer_id, module_id;
        if (fscanf(fp, "%lld,%lf,%lf,%lf,%d,%d,%d", &hit_id, &tx, &ty, &tz, &volume_id, &layer_id, &module_id) == -1) break;
        if (hit_id != n+1) cout << "Hit id's not as expected" << endl;
        x[n] = tx;
        y[n] = ty;
        z[n] = tz;
        v[n] = volume_id;
        l[n] = layer_id;
        m[n] = module_id;
        n++;
    }
    
    fclose(fp);
    cout << n << " hits" << endl;
    
    initHits(n,x,y,z,v,l,m);
}


//Calculate direction of each hit with cell's data
long Tracker::initHitDir() {
    
    long n(0);
    for (int hit_id = 1; hit_id < (int) hits.size(); hit_id++) {
        point m = meta[hit_id];
        Detector&d = detectors[int(m.x)*10000000+int(m.y)*10000+int(m.z)];
        
        //if (!hit_cells[hit_id].size()) cout << "Hit with zero cells" << endl;
        //if (metai[hit_id] < 18) continue;
        
        //Use linear regression for direction
        double mx = 0, my = 0, mw = 0;
        auto&cells = hit_cells[hit_id];
        for (auto&c : cells) {
            double w = c.second;
            double x = c.first.first*d.cell_w;
            double y = c.first.second*d.cell_h;
            mw += w;
            mx += x*w;
            my += y*w;
        }
        mx /= mw;
        my /= mw;
        double mxx = 0, mxy = 0, myy = 0;
        for (auto&c : cells) {
            double w = c.second;
            double x = c.first.first*d.cell_w-mx;
            double y = c.first.second*d.cell_h-my;
            mxx += x*x*w;
            myy += y*y*w;
            mxy += x*y*w;
            n++;
        }
        //Find eigenvector with minimum eigenvalue
        double a = mxx-myy, b = 2*mxy;
        double x = a+b+sqrt(a*a+b*b);
        double y =-a+b+sqrt(a*a+b*b);
        /*        if (0) {
         double lambda = (mxx+myy+sqrt(a*a+b*b))/2;
         cout << lambda << ' ' << (mxx*x+mxy*y)/x << ' ' << (mxy*x+myy*y)/y << endl;
         }
         */
        //Analytical formula for z
        double z = 2*d.d*(fabs(x)/d.cell_w+fabs(y)/d.cell_h+1e-8);
        x *= (cells.size()*1.-1.3);//1.3 != 1 was adjusted empirically
        y *= (cells.size()*1.-1.3);
        point d1(x,y,z), d2(x,y,-z);
        d1 = d.rx*d1.x+d.ry*d1.y+d.rz*d1.z;
        d2 = d.rx*d2.x+d.ry*d2.y+d.rz*d2.z;
        hit_dir[hit_id][0] = normalize(d1);
        hit_dir[hit_id][1] = normalize(d2);
        
        for (int k = 0; k < 2; k++)
            if (hit_dir[hit_id][k]*hits[hit_id] < 0)
                hit_dir[hit_id][k] = hit_dir[hit_id][k]*-1;
        
    }
    
    return n;
}


point Tracker::normalize(point a) {
    point ret = a*(1./point::dist(a));
    if (ret.z < 0) ret = ret*-1;
    return ret;
}
