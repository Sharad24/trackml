#ifndef _TRACKER_H_
#define _TRACKER_H_
// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#include "Graph.h"
#include "Point.h"
#include "Parameters.h"

#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <stack>
#include <string>

class Reconstruction;
class PolarModule;

//Structure for storing promising triples of hits
typedef std::tuple<int,int,int,float,float> triple;

//Geometry of layer
struct Layer {
    double minr, avgr, maxr;
    double minz, avgz, maxz;
    int count;
    int type;
    
    double var0, var1;
};

// Volumes: 7,8,9, 12,13,14, 16,17,18
// Volume 8 is innermost, it has modules 2,4,6,8 concentric cylinders outwards
struct Detector {
    int volume_id, layer_id, module_id;
    point c;
    point rx, ry, rz;
    double d, minw, maxw, h, cell_w, cell_h;
};

struct Particle // structure for truth particle info
{
    long long id;
    int type;
    double x;
    double y;
    double z;
    double r;
    double px;
    double py;
    double pz;
    double q;
    long hits;
    std::vector<int> hit;
};

class TNtuple;

class Tracker {
public:
    Tracker(int n=21001,const char *datapath=DATAPATH,const char *workpath=WORKPATH);
    ~Tracker();
    long getNumberHits() { return hits.size(); }
    double getScore() { return cscore; }
    void setThreshold1(double threshold) { threshold1 = threshold;}
    void setThreshold2(double threshold) { threshold2 = threshold;}
    void setThreshold3(double threshold) { threshold3 = threshold;}
    long long voxel(int hitid);
    void readBlacklist();
    void readWhitelist();
    void readTruth();
    void sortTracks();
    void readStarts();
    void readHits();
    void readHitsInitData();
    void readCells();
    void readDetectors(const char *path=WORKPATH);
    int *findTracks(int step, const char *path=WORKPATH);
    void initHits(int nhits,double *x,double *y,double *z,int *vol,int *lay,int *mod);
    void initCells(int ncells,int* id,int *ch0,int *ch1,double *value);
    long initGraphData();
    inline bool z_cmp(const int &a, const int &b);
    inline bool r_cmp(const int &a, const int &b);
    inline bool dist_cmp(const int &a, const int &b);
    inline bool track_cmp(const int a,const int b);
    void sortZ(std::vector<int> &hits);
    void sortR(std::vector<int> &hits);
    void sortT(std::vector<int> &hits);
    void sortDist(std::vector<int> &hits);
    static const int Tube = 0, Disc = 1;
    void setFilenumber(int n) {filenum = n; }
    int getFilenumber() { return filenum; }
    bool evaluation() { return eval; }
    void setDebug(bool d = true) { debug = d; }
    void setEvaluation(bool d = true) { eval = d; }
    void print(std::vector<int> const &input);
    void printshort(std::vector<int> const &input);
    long numberHits() { return hits.size(); }
    long long truthPart(int hitid) { return truth_part[hitid]; }
    void generateSliceGraphs();
    void readSliceGraphs();
    void writeSliceGraphs();
    void readTileGraphs(const char *path);
    void generateTileGraphs(const char *path);
    void writeGraph(const char *file, Graph<long long> &g);
    void readGraph(const char *file, Graph<long long> &g);
    void writeGraph(const char *file, Graph<std::pair<int,int> > &g);
    void readGraph(const char *file, Graph<std::pair<int,int> > &g);
    std::vector<std::vector<int> > swimmer();
    std::vector<point>& getHits() { return hits; }
    void draw(unsigned long nt,std::vector<point> &hits,std::map<int,std::vector<int> > &tracks);
    int shortid(int hitid) { return hitIDmap.find(hitid)==hitIDmap.end() ? hitid : hitIDmap[hitid]; }
    int slice(int ai);
    std::pair<int,int> graphHash(int hitid);
    void memoryFootprint();
    void initTasks();
    long addHits(std::vector<triple> &triples,int ai,int bi,Graph<long long> &g, double z0=0.0);
    long initHitDir();
    void initPolarModule();
    void initNeuralNetworks();

protected:
    int samepart(int a, int b);
    void initOrder();
    void initLayers();
    int getLayer(int volume_id, int layer_id);
    point normalize(point a);
    void writeSubmission(std::map<int, int>&assignment,std::string path, int filenum);
    double scorePairs(std::vector<triple>&pairs);
    double scoreTriples(std::vector<triple>&triples);
    double scorePaths4(std::vector<std::vector<int> > &paths);
    double scoreAssignment(std::map<int, int>&assignment);
    void scoreAssignment2(std::map<int,int> &assignment);
    int good_pair(int a, int b);
    double timer(double accuracy_mean);
    std::vector<std::pair<int, int> > allPairStarts();
    std::vector<triple> allTripleStarts();
    inline double dist(double x, double y) { return sqrt(x*x+y*y); }
    inline double dist2(double x, double y) { return (x*x+y*y); }
    inline point topolar(const point&dir, const point&ref, const point&refp) {
        return point(ref.x*dir.x+ref.y*dir.y, ref.x*dir.y-ref.y*dir.x, dir.z*refp.x);
    }
    inline point topolar(point p) {
        return point(sqrt(p.x*p.x+p.y*p.y), atan2(p.y,p.x), p.z);
    }
    bool sameTriple(triple a, triple b);
    bool intersection(int A, int B, int C, int D, point& ip);
    point intersection2d(int a, int b, int c, int d);
    Reconstruction *initRecoObjects(int n);
    
    // Data
    
    static bool eval;
    static bool debug;
    int filenum;
    std::string dataPath, workPath;
    double threshold1,threshold2,threshold3; // neural network cuts
    double cscore; // codalab score
    int assignment[MAXDIM]; //Assignment of hits to paths
    std::vector<point> hits; //hit position
    std::vector<point> polar; //hit position in polar / cylindrical coordinates
    std::map<long long, point> start_pos; //start position
    std::map<long long, point> start_mom; //start momentum
    std::map<long long, int> part_q; //start charge
    std::map<long long, int> part_hits; // = truth_tracks[particle_id].size()
    std::vector<Particle> particles; //true tracks
    std::map<long long,int> partIDmap; // create particle ID->index map
    std::vector<std::pair<std::pair<int, int>, double> > hit_cells[MAXDIM];
    point hit_dir[MAXDIM][2]; //The two possible directions of the hit according to the cell's data for each hit
    static std::map<int, Detector> detectors;
    int topo[LAYERS], itopo[LAYERS];
    Layer layer[LAYERS];
    double z_minr[LAYERS][4], z_maxr[LAYERS][4];
    double disc_z[LAYERS][4];
    std::vector<point> meta; //volume_id / layer_id / module_id
    std::vector<int> metai, metaz; //ordered layer id in [0,LAYERS), and classification of z for disc layers in [0,4)
    point truth_pos[MAXDIM], truth_mom[MAXDIM]; //truth position and momentum
    std::map<long long, std::vector<int> > truth_tracks; //truth hit ids in each track
    double truth_weight[MAXDIM]; //weighting of each hit
    long long truth_part[MAXDIM]; //particle this hit belongs to
    std::map<long long, double> part_weight; //weighting of each particle
    std::map<long long, std::map<int, double> > metai_weight; //weighting of each particle hit, also adding duplicates
    std::map<int,int> hitIDmap; // create hit ID->index map
    static Graph<long long> pgraph[NSLICES];   // graph to represent path information trained on slices
    static std::map<std::pair<int,int>, Graph<long long> > tgraph; // graph to represent path information trained on individual tracks
    std::stack<std::pair<int,Graph<long long> *> > tasks;
    
    Reconstruction *reco[NTHREADS];
    
    friend class PolarModule;
    friend class PolarModuleInternal;
};

#endif
