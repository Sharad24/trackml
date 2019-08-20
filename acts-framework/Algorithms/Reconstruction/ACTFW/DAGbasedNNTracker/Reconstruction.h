#ifndef _RECONSTRUCTION_H_
#define _RECONSTRUCTION_H_
// Neural Network based path reconstruction
// M.Kunze, Heidelberg University, 2018

#include "Graph.h"
#include "Point.h"

class PolarModule;
class XMLP;
class TFile;
class TNtuple;

#ifndef _TRACKER_H_
//Geometry of layer
struct Layer {
    double minr, avgr, maxr;
    double minz, avgz, maxz;
    int count;
    int type;
    
    double var0, var1;
};
#endif

//Structure for storing promising triples of hits
typedef std::tuple<int,int,int,float,float> triple;

class Reconstruction {
    
public:
    Reconstruction(int event,const char *workpath,std::map<std::pair<int,int>, Graph<long long> > &tripling,std::vector<point> &hits,std::vector<point> &polar,std::vector<point> &meta,std::vector<int> &metai,std::vector<int> &metaz,Layer (&layer)[LAYERS],double (&disc_z)[LAYERS][4],std::vector<std::pair<std::pair<int, int>, double> > (&hit_cells)[MAXDIM],point (&hit_dir)[MAXDIM][2]);
    ~Reconstruction();
    long long hitHash(int hitid);
    std::pair<int,int> graphHash(int hitid);
    void setThreshold1(double threshold) { threshold1 = threshold;}
    void setThreshold2(double threshold) { threshold2 = threshold;}
    void setThreshold3(double threshold) { threshold3 = threshold;}
    double getThreshold1() { return threshold1; }
    double getThreshold2() { return threshold2; }
    double getThreshold3() { return threshold3; }
    std::vector<triple> findCandidatesGraph(Graph<long long> &g);
    std::vector<triple> findTriplesGraph(Graph<long long> &g,std::vector<triple> &pairs);
    std::vector<std::vector<int> > findPaths(std::vector<triple>&triples);
    int extend3(int ai, int bi, int ci, int li, double target = 0.5);
    int prepareQuadrupleScore(int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&dirp, point target, double sign);
    int prepareQuadrupleScore(int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&bap, double sign);
    void initDensity3();
    double getDensity3(point&dp, point&xp, double tt, int li);
    double evaluateScore(int ci, point&dp, point&xp, point&bap);
    inline double dist(double x, double y) { return sqrt(x*x+y*y); }
    std::vector<std::vector<int> > addDuplicates(std::vector<std::vector<int> >&paths);
    int prepareDuplicateScore(int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&dirp);
    double findDensity(point&dp, point&xp, double target, int li);
    std::vector<std::vector<int> > prunePaths(std::vector<std::vector<int> >&paths);
    std::vector<std::vector<int> > findAssignment(std::vector<std::vector<int> >&paths, PolarModule *mod[LAYERS], int use_trash);
    double scorepathDensity(std::vector<int>&path);
    double scoreTripleDensity(int ai, int bi, int ci);
    double scoreQuadrupleDensity(int ai, int bi, int ci, int di);
    double scoreDuplicateDensity(int ai, int bi, int ci, int di);
    int prepareTripleScore(int ai, int bi, int li, point&d, point&dp, point&xp, point&bap, point target);
    int prepareTripleScore(int ai, int bi, int li, point&d, point&dp, point&xp, point&bap);
    bool getFeatures(int ai, int bi, float* feature, point p=point(0,0,0));
    double dir_miss(int ai, int bi);
    double wdistr(double r1, double dr, double az, double dz, double w);
    double wdist(point&a, point&d, double w);
    double zdist(point&a, point&b);
    double zdist2(point&a, point&b);
    double xyVertex(int ai, int bi);
    double zVertex(int ai, int bi);
    std::pair<double,double> zVertexScore(int ai, int bi);
    point getVertex(int ai, int bi);
    point distBetweenLines(point &p1, point &p2, point &p3, point &p4);
    double recall1(Graph<long long> &g,int a, int b, double f1=0.0, double f2=0.0, double z0=0.0);
    double recall2(Graph<long long> &g,int a, int b, double f1=0.0, double f2=0.0, double f3=0.0, double z0=0.0);
    double recall3(Graph<long long> &g,int a, int b, int c, double f1=0.0, double z0=0.0);
    double recall3(Graph<long long> &g,triple &t);
    double scoreTriple(int ai, int bi, int ci);
    static void setPhiResolution(int pr) { phires = pr; }
    static void setThetaResolution(int tr) { theres = tr; }
    static int phiResolution() { return phires; }
    static int thetaResolution() { return theres; }

private:

    void circle(point&a, point&b, point&c, point&p, double&r);
    inline double dist2(double x, double y) { return (x*x+y*y); }
    double field(double z);
    inline int getIndex(int&li, double x);
    void loadAdjacent(const char *filename);
    
    // Data
    
    int eventnum;
    std::string workPath;
    double threshold1, threshold2, threshold3;
    static int phires, theres;
    std::map<std::pair<int,int>, Graph<long long> > &tgraph; // graph to represent path information
    std::vector<point> &hits; //hit position
    std::vector<point> &polar; //hit position in polar / cylindrical coordinates
    std::vector<point> &meta; //detector,layer,module
    std::vector<int> &metai; //layer id
    std::vector<int> &metaz; //classification of z for disc layers in [0,4)
    Layer (&layer)[LAYERS];
    double (&disc_z)[LAYERS][4];
    std::vector<std::pair<std::pair<int, int>, double> > (&hit_cells)[MAXDIM];
    point (&hit_dir)[MAXDIM][2];

    static const int Tube = 0, Disc = 1;
    static int next_layer[LAYERS][LAYERS]; //Number of directly adjacent hits in layers
    static const int adj_thres = 450; //Threshold in next_layer to consider layers adjacent
    int match[MAXDIM]; //reused global array for storing matching hits from PolarModule->getNear function
    const double Bfield = 1673.; //Empirical field strengh, to scale the momentum
    const double stretch = 0.02; //Decides how curved "approximately straight" helices are supposed to be
    std::vector<double> sorted_hits[LAYERS]; //List of coordinates for each layer, used to estimate outlier densities
    static const int crude_steps = 1<<10; //Acceleration look-up table for faster lookup in sorted_hits
    int crudeIndex[LAYERS][crude_steps]; //Scaling to get sorted_hits in range [0,crude_steps)
    std::pair<double, double> crudeIndex_a[LAYERS]; //Accumulated polynomial coefficients of second order polynomial for O(1) density look-up. Second dimension (20000) must be bigger than number of hits in the most populated layer
    double poly[LAYERS][20000][3];
    
    static TFile *file;
    static TNtuple *ntuple1,*ntuple2,*ntuple3;
};

#endif
