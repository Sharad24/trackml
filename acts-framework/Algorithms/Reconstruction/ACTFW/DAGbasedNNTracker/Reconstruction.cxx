// Neural Network based path reconstruction
// M.Kunze, Heidelberg University, 2018

#include <cmath>
#include <queue>
#include <sstream>
#include <iomanip>
#include "Parameters.h"
#include "Reconstruction.h"
#include "PolarModule.h"
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

extern PolarModule *mod[LAYERS];
int Reconstruction::next_layer[LAYERS][LAYERS];
TFile *Reconstruction::file(NULL);
TNtuple *Reconstruction::ntuple1(NULL);
TNtuple *Reconstruction::ntuple2(NULL);
TNtuple *Reconstruction::ntuple3(NULL);
int Reconstruction::phires(PHIRESOLUTION);
int Reconstruction::theres(THETARESOLUTION);

using namespace std;

long long Reconstruction::hitHash(int hitid) {
    long l=metai[hitid];
    long m=meta[hitid].z;
    auto h = graphHash(hitid);
    long long index = (((long)h.first)<<32) | (((long)h.second)<<24) | (l<<16) | m;
    //long long index = (l<<16) | m;
    //long long index = l*MODULES + m;
    return index;
}

Reconstruction::Reconstruction(int event,const char *workpath,std::map<std::pair<int,int>, Graph<long long> > &t,vector<point> &h,std::vector<point> &p,vector<point> &m,std::vector<int> &mi,std::vector<int> &mz,Layer (&l)[LAYERS],double (&dz)[LAYERS][4],vector<pair<pair<int, int>, double> > (&hc)[MAXDIM],point (&hd)[MAXDIM][2]) : eventnum(event), workPath(workpath), tgraph(t), hits(h), polar(p), meta(m), metai(mi), metaz(mz), layer(l), disc_z(dz), hit_cells(hc), hit_dir(hd)
{
    static bool initialized(false);
    
    if (!initialized) {
        string filename = workPath + "/adjacency";
        loadAdjacent(filename.c_str());
        initialized = true;
    }
    
    initDensity3();
    
    threshold1 = 0.1;
    threshold2 = 0.1;
    threshold3 = 0.3;
}

Reconstruction::~Reconstruction()
{
}

std::pair<int,int> Reconstruction::graphHash(int hitid) {
    
    point pol = polar[hitid];
    
    double rt = pol.x;
    double ph = pol.y;
    double z0 = pol.z;
    double th = asinh(z0/rt); // slope (-4...+4)
    int i1 = (int) (phires*0.15*(M_PI+ph));
    int i2 = (int) (theres*0.1*(5-th));
    
    return make_pair(i1,i2);
}


// Calculate a vertex for two points
double Reconstruction::xyVertex(int ai, int bi)
{
    point &a = hits[ai];
    point &b = hits[bi];
    point d = a - b;
    double ppxy = a.x*a.x+a.y*a.y;
    double pdxy = a.x*d.x+a.y*d.y;
    double ddxy = d.x*d.x+d.y*d.y;
    if (ddxy==0) return 1.E3;
    return sqrt(ppxy-pdxy*pdxy/ddxy);
}


double Reconstruction::zVertex(int ai, int bi)
{
    point &a = hits[ai];
    point &b = hits[bi];
    point d = a - b;
    double pp = a.x*a.x+a.y*a.y+a.z*a.z;
    double pd = a.x*d.x+a.y*d.y+a.z*d.z;
    double dd = d.x*d.x+d.y*d.y+d.z*d.z;
    if (dd==0) return 1.E3;
    return sqrt(pp-pd*pd/dd);
}


pair<double,double> Reconstruction::zVertexScore(int ai, int bi)
{
    double z(0);
    double dz = 10.0;
    double minscore = 1.E3;
    hits[0] = point(0,0,0);
    if (hits[ai].z<0) dz = -dz; // negative direction
    for (int i=0;i<16;i++) {
        hits[0].z = i*dz;
        double score = scoreTriple(0,ai,bi);
        if (score < minscore) {
            z = hits[0].z;
            minscore = score;
        }
    }
    return make_pair(z,minscore);
}


// Look for the point of closest approach to the origin
point Reconstruction::getVertex(int ai, int bi)
{
    point origin(0,0,0);
    point linePnt = hits[ai];
    point line = hits[ai] - hits[bi];
    double factor = 1./point::norm(line);
    point lineDir = line*factor;
    point v = origin - linePnt;
    double d = point::dot(v,lineDir);
    return linePnt + lineDir * d;
}

// Look for seeding pair/triple candidates by hit pair/triple combinations in consecutive layers
vector<triple> Reconstruction::findCandidatesGraph(Graph<long long> &g)
{
    static const int n(STARTLAYERS); // number of seeding layers
    static const int startlayer[48] = {0,11,4,18,1,5,12,13,6,2,3,19,20,7,14,21,24,36,15,8,22,9,16,38,40,42,26,28,30,25,37,10,17,23,32,34,44,46,27,39,29,41,31,43,33,45,35,47}; // seeding layers
    static const double vertexscore[48] = {470,590,430,430,140,100,80,210,100,160,260,430,560,110,110,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE,VERTEXSCORE}; // max helix score
    static const double offset[48] = {0.0,0.0,0.0,0.2,0.1,0.0,0.0,0.0,0.0,0.2,0.3,0.2,0.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // adaptation wrt. recall distribution per layer
    
    vector<triple> candidates;
    
    for (int i = 0; i < n; i++)
    {
        int layer1 = startlayer[i];
        
        for (auto start : g.hash(layer1)) { // all modules in first layer
            const auto &edgelist = g.edges(start);
            if (edgelist.size() == 0) continue;
            
            for (auto &edge : edgelist) {
                long nextindex = edge.first;
                double vz = 0.0;
#ifdef VERTEXCUT
                vz = edge.second; // get z vertex from graph
                if (abs(vz) > VERTEXCUT) continue;
#endif
                
                for (auto ai : g.data(start)) { // all hits in module
                    
                    auto &p = g.data(nextindex); // all hits in following modules
                    if (p.size() == 0) continue;
                    
                    for (auto bi:p)
                    {
                        int l1 = metai[ai];
                        int l2 = metai[bi];
                        if (l1 == l2) continue; // Same layer (Double hits)
                        
                        double xy0 = xyVertex(ai,bi); // check the radial distance from origin
                        if (i<3 && xy0 > VERTEXCUTXY) continue;
                        double z0 = zVertex(ai,bi); // check the z distance from origin
                        if (z0 > VERTEXCUTZ) continue;
                        //double zdist = zdist2(hits[ai],hits[bi]);
                        //if (zdist > VERTEXCUTZ) continue;
                        hits[0] = point(0,0,vz); // origin
                        double score = scoreTriple(0,ai,bi); // helix score wrt. origin
                        if (score > vertexscore[i]) continue;
                        double dir1 = dir_miss(ai,bi);
                        double dir2 = dir_miss(bi,ai);
                        double recall(1.0);
                        bool tube = (l1<4) || (l1>=18 && l1<=23);

                        if (!tube) {
                            recall = recall1(g,ai,bi,dir1,dir2,vz); // Search for hit pairs using coordinates and hit direction
                            if (recall < threshold1+offset[i]) continue;
                        }
                        else {
                            recall = recall2(g,ai,bi,dir1,dir2,score,vz); // Search for hit pairs using coordinates, hit direction, and helix score
                            if (recall < threshold2+offset[i]) continue;
                        }
                        
                        candidates.push_back(triple(ai,bi,0,recall,vz)); // Note the good pair candidates

                    }
                    
                }
                
            }
            
        }
    }
    
    return candidates;
}


// Generate tracklets of 3 points wrt. the first point in seed
std::vector<triple> Reconstruction::findTriplesGraph(Graph<long long> &g,vector<triple> &pairs) {
    
    static const int target(2);
    vector<triple> triples;
    
    int n = (int) pairs.size();
    for (int i=0;i<n;i++)
    {
        vector<triple> batchtriples;
        auto &pa = pairs[i];
        int ai = get<0>(pa);
        int bi = get<1>(pa);
        
#ifdef COMBINEDMETHOD
        pair<int,int> index = graphHash(bi);
        auto *gg = &tgraph[index];
#else
        Graph<long long> *gg = &g;
#endif
        
        // Search for triples
        long indexg = hitHash(bi);
        const auto &edgelist = gg->edges(indexg);
        if (edgelist.size() == 0) continue;
        
        for (auto &edge : edgelist) {
            long nextindex = edge.first;
            double vz = 0.0; //edge.second;

            auto &p = g.data(nextindex); // all hits in following modules
            if (p.size() == 0) continue;
            
            vector<triple> v;
            for (auto ci:p)
            {
                if (ci==ai || ci==bi) continue; // Same hit
                double score = scoreTriple(ai,bi,ci); // Check helix propagation
                if (score>SWEIGHT) continue; // bad helix hypothesis
                triple t(ai,bi,ci,score,vz);
                v.push_back(t);
            }
            
            // Sort the result by score using lambda function
            sort(begin(v), end(v),[](triple const &t1, triple const &t2) {return get<3>(t1) < get<3>(t2);} );
            
            if (v.size()>target) v.resize(target);
            
            for (auto &t : v) {
                double recall = recall3(*gg,t); // Check quality with ANN
                if (recall<threshold3) continue;
                triples.push_back(t); // Note the good triple candidates
            }
        }
    }
    
    return triples;
}


// Recall function for 2 points (cylinder coordinates)
double Reconstruction::recall1(Graph<long long> &g,int a, int b, double f1, double f2, double z0)
{
    point &p1 = polar[a];
    point &p2 = polar[b];
    
#ifdef USETMVA
    float *x1 = g.getX1();
#ifdef FOLDEDINPUT1
    x1[0] = p1.x; // rz1
    x1[1] = fabs(fabs(p1.y)-M_PI_2); // phi1
    x1[2] = fabs(p1.z); // z1
    x1[3] = p2.x; // rz2
    x1[4] = fabs(fabs(p2.y)-M_PI_2); // phi2
    x1[5] = fabs(p2.z); // z2
#else
    x1[0] = p1.x;
    x1[1] = p1.y;
    x1[2] = p1.z;
    x1[3] = p2.x;
    x1[4] = p2.y;
    x1[5] = p2.z;
#endif
    x1[6] = f1; // dirmiss1
    x1[7] = f2; // dirmiss2
    
#ifdef TMVAREADER
    TMVA::Reader *reader = g.getReader1();
    double recall = reader->EvaluateMVA("MLP method");
#else
    vector<double> x;
    for (int i=0;i<8;i++) x.push_back(x1[i]);
    double recall = g.getNet1()->GetMvaValue(x);
#endif
    
#else

    float x[10];
#ifdef FOLDEDINPUT1
    x[0]    = p1.x*0.001;               // rz1 [m]
    x[1]    = fabs(fabs(p1.y)-M_PI_2);  // phi1
    x[2]    = fabs(p1.z-z0)*0.001;      // z1 [m]
    x[3]    = p2.x*0.001;               // rz2 [m]
    x[4]    = fabs(fabs(p2.y)-M_PI_2);  // phi2
    x[5]    = fabs(p2.z-z0)*0.001;      // z2 [m]
#else
    x[0]    = p1.x*0.001;   // rz1 [m]
    x[1]    = p1.y;         // phi1
    x[2]    = p1.z*0.001;   // z1 [m]
    x[3]    = p2.x*0.001;   // rz2 [m]
    x[4]    = p2.y;         // phi2
    x[5]    = p2.z*0.001;   // z2 [m]
#endif
    x[6]    = f1;           // feature
    x[7]    = f2;           // feature

    double recall = g.net1()->Recallstep(x)[0];
#endif
    
    return recall;
}


// Recall function for 2 points (folded cylinder coordinates)
double Reconstruction::recall2(Graph<long long> &g,int a, int b, double f1, double f2, double f3, double z0)
{
    point &p1 = polar[a];
    point &p2 = polar[b];

#ifdef USETMVA
    float *x2 = g.getX2();
#ifdef FOLDEDINPUT2
    x2[0] = p1.x; // rz1
    x2[1] = fabs(fabs(p1.y)-M_PI_2); // phi1
    x2[2] = fabs(p1.z); // z1
    x2[3] = p2.x; // rz2
    x2[4] = fabs(fabs(p2.y)-M_PI_2); // phi2
    x2[5] = fabs(p2.z); // z2
#else
    x2[0] = p1.x;
    x2[1] = p1.y;
    x2[2] = p1.z;
    x2[3] = p2.x;
    x2[4] = p2.y;
    x2[5] = p2.z;
#endif
    x2[6] = f1; // dirmiss1
    x2[7] = f2; // dirmiss2
    x2[8] = log(f3); // score
    
#ifdef TMVAREADER
    TMVA::Reader *reader = g.getReader2();
    double recall = reader->EvaluateMVA("MLP method");
#else
    vector<double> x;
    for (int i=0;i<9;i++) x.push_back(x2[i]);
    double recall = g.getNet2()->GetMvaValue(x);
#endif
    
#else

    float x[10];
#ifdef FOLDEDINPUT2
    x[0]    = p1.x*0.001;               // rz1 [m]
    x[1]    = fabs(fabs(p1.y)-M_PI_2);  // phi1
    x[2]    = fabs(p1.z-z0)*0.001;      // z1 [m]
    x[3]    = p2.x*0.001;               // rz2 [m]
    x[4]    = fabs(fabs(p2.y)-M_PI_2);  // phi2
    x[5]    = fabs(p2.z-z0)*0.001;      // z2 [m]
#else
    x[0]    = p1.x*0.001;   // rz1 [m]
    x[1]    = p1.y;         // phi1
    x[2]    = p1.z*0.001;   // z1 [m]
    x[3]    = p2.x*0.001;   // rz2 [m]
    x[4]    = p2.y;         // phi2
    x[5]    = p2.z*0.001;   // z2 [m]
#endif
    x[6]    = f1;           // feature
    x[7]    = f2;           // feature
    x[8]    = 0.001*f3;     // feature
    
    double recall = g.net2()->Recallstep(x)[0];
#endif
    
    return recall;
}


// Recall function for 3 points
double Reconstruction::recall3(Graph<long long> &g,triple &t)
{
    int ai, bi, ci;
    double r3,z0;
    tie(ai,bi,ci,r3,z0) = t;
    r3 = recall3(g,ai,bi,ci,r3,z0);
    get<3>(t) = r3;
    return r3;
}


// Recall function for 3 points (cylinder coordinates)
double Reconstruction::recall3(Graph<long long> &g,int a, int b, int c, double f1, double z0)
{
    point &p1 = polar[a];
    point &p2 = polar[b];
    point &p3 = polar[c];
    
#ifdef USETMVA
#ifdef TMVAREADER
    TMVA::Reader *reader = g.getReader3();
    float *x = g.getX3();
    x[0] = p1.x; // rz1
    x[1] = p1.y; // phi1
    x[2] = p1.z; // z1
    x[3] = p2.x; // rz2
    x[4] = p2.y; // phi2
    x[5] = p2.z; // z2
    x[6] = p3.x; // rz3
    x[7] = p3.y; // phi23
    x[8] = p3.z; // z3
    x[9] = log(f1); // score
    double recall = reader->EvaluateMVA("MLP method");
#else
    vector<double> x;
    x.push_back(p1.x);
    x.push_back(p1.y);
    x.push_back(p1.z);
    x.push_back(p2.x);
    x.push_back(p2.y);
    x.push_back(p2.z);
    x.push_back(p3.x);
    x.push_back(p3.y);
    x.push_back(p3.z);
    x.push_back(log(f1));
    double recall = g.getNet3()->GetMvaValue(x);
#endif
    
#else

    float x[12];
#ifdef FOLDEDINPUT3
    x[0]    = p1.x*0.001;               // rz1 [m]
    x[1]    = fabs(fabs(p1.y)-M_PI_2);  // phi1
    x[2]    = fabs(p1.z-z0)*0.001;      // z1 [m]
    x[3]    = p2.x*0.001;               // rz2 [m]
    x[4]    = fabs(fabs(p2.y)-M_PI_2);  // phi2
    x[5]    = fabs(p2.z-z0)*0.001;      // z2 [m]
    x[6]    = p3.x*0.001;               // rz3 [m]
    x[7]    = fabs(fabs(p3.y)-M_PI_2);  // phi3
    x[8]    = fabs(p3.z-z0)*0.001;      // z3 [m]
#else
    x[0]    = p1.x*0.001;        // rz1 [m]
    x[1]    = p1.y;              // phi1
    x[2]    = (p1.z-z0)*0.001;   // z1 [m]
    x[3]    = p2.x*0.001;        // rz2 [m]
    x[4]    = p2.y;              // phi2
    x[5]    = (p2.z-z0)*0.001;   // z2 [m]
    x[6]    = p3.x*0.001;        // rz3 [m]
    x[7]    = p3.y;              // phi3
    x[8]    = (p3.z-z0)*0.001;   // z3 [m]
#endif
    x[9]    = f1;           // feature
    double recall = g.net3()->Recallstep(x)[0];
    
#endif
    
    return recall;
}


// The following code has been adopted from Johan Sokrates Wind's award winning trackml Kaggle contribution
// Parameters have been slightly optimized and at some places values are checked to make the code more robust
// Triples have been substituted by standard C++ tuples
// see https://www.kaggle.com/c/trackml-particle-identification/discussion/63249

//Find circle with center p, radius r, going through a, b, and c (in xy plane)
void Reconstruction::circle(point&a, point&b, point&c, point&p, double&r) {
    double ax = a.x-c.x, ay = a.y-c.y, bx = b.x-c.x, by = b.y-c.y;
    double aa = ax*ax+ay*ay, bb = bx*bx+by*by;
    double x = ax*by-ay*bx;
    if (x==0) {
        r = 0;
        return;
    }
    double idet = .5/x;
    p.x = (aa*by-bb*ay)*idet;
    p.y = (ax*bb-bx*aa)*idet;
    p.z = 0;
    r = dist(p.x, p.y);
    p.x += c.x;
    p.y += c.y;
}


//Score triple based on the deviation from a perfect helix, no prior that it should be straight
double Reconstruction::scoreTriple(int ai, int bi, int ci) {
    point center;
    double radius;
    if (ai==bi || ai==ci || bi==ci) return 1e3;
    circle(hits[ai], hits[bi], hits[ci], center, radius);
    if (radius==0) return 1e3;
    
    point cb = hits[ci]-hits[bi];
    point ba = hits[bi]-hits[ai];
    double ang_cb = asin(dist(cb.x, cb.y)*.5/radius)*2;
    double ang_ba = asin(dist(ba.x, ba.y)*.5/radius)*2;
    if (radius != radius || fabs(radius) > 1e50) {
        ang_cb = dist(cb.x, cb.y);
        ang_ba = dist(ba.x, ba.y);
    }
    if (ba.z*cb.z < 0) ang_ba *= -1;
    
    if (ang_cb-ba.z == 0) return 1.E3;
    if (ang_ba == 0) return 1.E3;
    
    double y = ang_cb ? (fabs(cb.z*ang_ba/ang_cb-ba.z)) : 1e3;
    double z = ang_ba ? (fabs(cb.z-ba.z*ang_cb/ang_ba)) : 1e3;
    double score = min(y, z);
    
    return score;
}


//Expand triples into paths, do this by expanding helices through the outermost 3 point in each direction repeatedly
vector<vector<int> > Reconstruction::findPaths(vector<triple>&triples) {
    
    vector<vector<int> > paths;
    static const int minpath = 3;
    
    int n = (int) triples.size();
    for (int i=0;i<n;i++)
    {
        vector<vector<int> > batchpaths;
        auto &t = triples[i];
        vector<int> v;
        int ai, bi, ci;
        double r3,z0;
        tie(ai,bi,ci,r3,z0) = t;
        int misses = 0;
        
        for (int li = metai[ai]-1; li >= 0; li--) {
            if (next_layer[li][metai[ai]] < adj_thres) continue;
            int di = extend3(ci, bi, ai, li);
            if (di > 0) {
                v.push_back(di);
                ci = bi;
                bi = ai;
                ai = di;
                misses = 0;
            } else if (di == -1) misses++;
            if (misses == 2) break;
        }
        reverse(v.begin(), v.end());
        
        tie(ai,bi,ci,r3,z0) = t;
        
        v.push_back(ai);
        v.push_back(bi);
        v.push_back(ci);
        
        misses = 0;
        for (int li = metai[ci]+1; li < LAYERS; li++) {
            if (next_layer[metai[ci]][li] < adj_thres) continue;
            int di = extend3(ai, bi, ci, li);
            if (di > 0) {
                v.push_back(di);
                ai = bi;
                bi = ci;
                ci = di;
                misses = 0;
            } else if (di == -1) misses++;
            if (misses == 2) break;
        }
        
        v.shrink_to_fit();
        
        if (v.size()>minpath) batchpaths.push_back(v);
        
        paths.insert(paths.end(),batchpaths.begin(),batchpaths.end()); // append the candidates
        
    }
    
    // Remove duplicate entries
    sort(paths.begin(),paths.end());
    paths.erase(unique(paths.begin(),paths.end() ),paths.end());
    
    return paths;
}


//Extend the helix going through hits with ids "ai", "bi", "ci", to layer "li". Do this by looking at the intersection with the layer, and expecting around "target" continuation for each outlier triple. "li" must be after metai[ci]
int Reconstruction::extend3(int ai, int bi, int ci, int li, double target) {
    point d, dp, xp, bap;
    if (prepareQuadrupleScore(ai,bi,ci,li,d,dp,xp,bap,1)) return -2;
    
    double tt = 400;//findDensity(dp, xp, target, li);
    double mins = 1e9;//target;//1e9;
    
    double fac = getDensity3(dp, xp, tt, li)/tt;
    tt = target/fac;
    int matches = mod[li]->getNear(dp, xp, bap, tt, match);
    int mini = -1;
    for (int i = 0; i < matches; i++) {
        int ti = match[i];
        if (ti >= (int)hits.size()){
            cout << "ERROR: " << ti << endl;
            continue;
        }
        //if (assignment[ti]) continue;
        double s = evaluateScore(ti,dp,xp,bap)*fac;
        if (s < mins) {
            mins = s;
            mini = ti;
        }
    }
    return mini;
}


//Default target is average position of layer
int Reconstruction::prepareQuadrupleScore(int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&bap, double sign = 1) {
    Layer l = layer[li];
    point target(l.avgr, 0, l.avgz);
    return prepareQuadrupleScore(ai, bi, ci, li, d, dp, xp, bap, target, sign);
}


//Similar to prepareTripleScore, but now we extend the helix that passes through hits with id "ai", "bi", "ci". Assumes li > metai[ci] if sign = 1, and metai[bi] < li < metai[ci] if sign = -1
int Reconstruction::prepareQuadrupleScore(int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point &dirp, point target, double sign = 1) {
    Layer l = layer[li];
    
    point p;
    double r, ir;
    
    point c = hits[ci];
    point cb = hits[ci]-hits[bi];//c-b;
    double ang_cb;
    
    if (1) { //Take into account varying magnetic field strength
        point a = hits[ai], b = hits[bi], c = hits[ci];
        if (l.type == Disc && (c.z-b.z)*(target.z-c.z)*sign < 0) return -1;
        double B1 = field((a.z+b.z)*.5), B2 = field((b.z+c.z)*.5), B3;
        if (l.type == Disc) B3 = field((c.z+target.z)*.5);
        else B3 = field(c.z);
        //B1 = B2 = B3 = 1;
        double ax = b.x-a.x, ay = b.y-a.y, bx = c.x-b.x, by = c.y-b.y;
        double aa = ax*ax+ay*ay, dot = ax*bx+ay*by, cross = ax*by-ay*bx;
        double alpha = B2/(2*B3), beta = (-B1*aa-B2*dot)/(2*cross*B3);
        //alpha *= -1;
        double rx = alpha*bx-beta*by, ry = alpha*by+beta*bx;
        p = point(c.x-rx, c.y-ry, 0);
        r = dist(rx, ry);
        ir = 1./r;
        ang_cb = B3/B2*asin(dist(cb.x, cb.y)*.5*ir*B2/B3)*2;
    }
    
    //exit(0);
    
    const double slack = 1.00;
    
    xp = point(0,0,0); //Circle for now
    
    if (l.type == Tube) {
        double RR = target.x*target.x, pp = dist2(p.x, p.y);
        double s = .5+(RR-r*r)/(2*pp);
        double sq = RR/pp-s*s;
        
        if (sq < 0) return -1;
        
        double t = sqrt(sq);
        if (p.y*c.x-p.x*c.y < 0) t *= -1;
        d.x = p.x*s+p.y*t;
        d.y = p.y*s-p.x*t;
        
        point dc = d-c;
        double A = dist(dc.x, dc.y);
        double B = A*.5*ir;
        double ang_dc = asin(B)*2;
        if (dc.x*cb.x+dc.y*cb.y < 0) ang_cb *= -1;
        
        d.z = c.z+cb.z*ang_dc/ang_cb;
        
        if (!(d.z > l.minz*slack && d.z < l.maxz*slack)) return -1;
        
        point dir;
        double s_ = target.x/pp, t_ = s_*(1-s)/t;
        dir.x = p.x*s_+p.y*t_;
        dir.y = p.y*s_-p.x*t_;
        dir.z = (dc.x*dir.x+dc.y*dir.y)*ir*cb.z/(ang_cb*A*sqrt(1-B*B));
        
        dp = point(dist(d.x,d.y), atan2(d.y, d.x), d.z);
        dirp = point(d.x*dir.x+d.y*dir.y, d.x*dir.y-d.y*dir.x, dir.z*dp.x);
        //cout << dirp << endl; //dirp.x = l.avgr
        
        dirp = dirp*(1./dirp.x);
        
    } else if (l.type == Disc) {
        d.z = target.z;
        double fac = ang_cb/cb.z;
        double ang_dc = (d.z-c.z)*fac;
        
        double sa = sin(ang_dc), ca = cos(ang_dc);
        
        double rx = c.x-p.x, ry = c.y-p.y;
        double cross = rx*cb.y-ry*cb.x;
        if (cross < 0) sa *= -1;
        
        d.x = ca*rx-sa*ry+p.x;
        d.y = sa*rx+ca*ry+p.y;
        
        
        point dir;
        dir.x =-fac*(rx*sa+ry*ca);
        dir.y = fac*(rx*ca-ry*sa);
        dir.z = cross < 0 ? -1 : 1;
        
        
        dp = point(dist(d.x,d.y), atan2(d.y, d.x), d.z);
        
        if (!(dp.x > l.minr*(1./slack) && dp.x < l.maxr*slack)) return -1;
        
        dirp = point(d.x*dir.x+d.y*dir.y, d.x*dir.y-d.y*dir.x, dir.z*dp.x);
        //cout << dirp << endl; //dirp.x = l.avgr
        
        dirp = dirp*(1./dirp.z);
    }
    return 0;
}


//Init everything needed for fast density calculations, includes most global variables above
void Reconstruction::initDensity3() {
    vector<int>*tube = new vector<int>[48]();
    
    for (int i = 1; i < (int)hits.size(); i++) {
        //if (!assignment[i])
        tube[metai[i]].push_back(i);
    }
    
    for (int li = 0; li < 48; li++) {
        sorted_hits[li].clear();
        if (layer[li].type == Tube) {
            for (int i : tube[li])
                sorted_hits[li].push_back(hits[i].z);
        } else {
            for (int i : tube[li])
                sorted_hits[li].push_back(polar[i].x);
        }
        sorted_hits[li].push_back(-1e50);
        sorted_hits[li].push_back(1e50);
        sort(sorted_hits[li].begin(), sorted_hits[li].end());
        
        double minx = *next(sorted_hits[li].begin())-1e-8;
        double maxx = *next(sorted_hits[li].rbegin())+1e-8;
        //cout << maxx << ' ' << minx << endl;
        double f = crude_steps/(maxx-minx);
        crudeIndex_a[li] = make_pair(f, -minx*f);
        
        for (int i = 0; i < crude_steps; i++) {
            double x = (i+.5)/f+minx;
            crudeIndex[li][i] = (int) (upper_bound(sorted_hits[li].begin(), sorted_hits[li].end(), x)-sorted_hits[li].begin());
        }
        double acc[3] = {};
        for (int i = 1; i < (int)sorted_hits[li].size(); i++) {
            for (int j = 0; j < 3; j++) poly[li][i][j] = acc[j];
            double x = sorted_hits[li][i];
            for (int j = 0; j < 3; j++)
                acc[j] += pow(x, j);
        }
    }
    delete[]tube;
}


//Get expected number of hits on layer "li" in area (in polar/cylindrical coordnates) spanned by (p-dp)^2-dot(p-dp, xp)^2 < tt
double Reconstruction::getDensity3(point&dp, point&xp, double tt, int li) {
    Layer l = layer[li];
    
    double x0, dx, dy;
    if (l.type == Tube) {
        x0 = dp.z;
        dx = xp.z;
        dy = xp.y;
    } else {
        x0 = dp.x;
        dx = xp.x;
        dy = xp.y;
    }
    double b = tt*(1-dy*dy)/(1-dx*dx-dy*dy);
    double a = sqrt(1-dx*dx-dy*dy)/((1-dy*dy)*M_PI*dp.x);
    double rx = sqrt(b);
    
    int ai = getIndex(li, x0-rx);
    int bi = getIndex(li, x0+rx);
    if (bi-ai > 10) {//Approximate integration by 2. order polynomial approximation to half disc
        //cout << ai << ' ' << bi << endl;
        const double A = 21*M_PI/64., B = -15*M_PI/64.;
        double ib = 1./b;
        double c0 = A+B*x0*x0*ib, c1 = -2*B*ib*x0, c2 = B*ib;
        double ret =
        ((poly[li][bi][0]-poly[li][ai][0])*c0+
         (poly[li][bi][1]-poly[li][ai][1])*c1+
         (poly[li][bi][2]-poly[li][ai][2])*c2)*a*rx;
        return max(ret,0.);
    } else { //Exact integration, uses half disc
        double density = 0;
        for(int i = ai; i < bi; i++) {
            double x = sorted_hits[li][i]-x0;
            double h = a*sqrt(b-x*x);/// *it;
            //cout << h << endl;
            density += h;
        }
        return density;
    }
}


//Use the prepared "dp", "xp", "bap" and return the area that is closer to the collision line (taking into account xp for elliptic behaviour) compared to the hit with id "ci"
double Reconstruction::evaluateScore(int ci, point&dp, point&xp, point&bap) {
    point&r = polar[ci];
    point err = r-dp;
    if (err.y > M_PI) err.y -= M_PI*2;
    if (err.y <-M_PI) err.y += M_PI*2;
    err.y *= dp.x;
    
    err = err-bap*(layer[metai[ci]].type == Disc ? err.z : err.x);
    double r2 = err*err-pow(err*xp, 2);
    return r2;
}


//Approximate magnetic field strengh as a function of z coordinate, decays drastically near the ends
double Reconstruction::field(double z) {
    z *= 1./2750;
    double z2 = z*z;
    return 1.002-z*3e-2-z2*(0.55-0.3*(1-z2));
}


// O(1) indexing in sorted_hits
// Functionally similar to "upper_bound(sorted_hits[li].begin(), sorted_hits[li].end(), x)-sorted_hits[li].begin();"
inline int Reconstruction::getIndex(int&li, double x) {
    int ci = x*crudeIndex_a[li].first+crudeIndex_a[li].second;
    ci = min(crude_steps-1, max(0, ci));
    int i = crudeIndex[li][ci];
    
    //Might segfault sometimes :)
    while (x >= sorted_hits[li][i]) i++;
    while (x < sorted_hits[li][i-1]) i--;
    
    return i;//max(0,min(i,int(sorted_hits[li].size())-1));
}


//Initialize next_layer
void Reconstruction::loadAdjacent(const char *filename) {
    FILE*fp = fopen(filename, "r");
    if (!fp) {
        cout << "Could not open " << filename << endl;
        exit(0);
    }
    for (int i = 0; i < LAYERS; i++)
        for (int j = 0; j < LAYERS; j++)
            if (!fscanf(fp, "%d", &next_layer[i][j])) cout << "Error reading adjacent file" << endl;
    fclose(fp);
}


//Attempt to add all duplicate hits (hits on same layer as an already added hit in the path) to the paths
vector<vector<int> > Reconstruction::addDuplicates(vector<vector<int> >&paths) {
    vector<vector<int> > extended;
    
    int n = (int) paths.size();
    for (int i=0;i<n;i++)
    {
        auto &path = paths[i];
        if (path.size() < 3) continue;
        
        vector<int> ext;
        for (int i = 0; i < (int)path.size(); i++) {
            ext.push_back(path[i]);
            int ai, bi, ci;
            if (i < (int)path.size()-2) {
                ai = path[i];
                bi = path[i+1];
                ci = path[i+2];
            }
            else if (i == (int)path.size()-1) {
                ai = path[i-2];
                bi = path[i-1];
                ci = path[i];
            }
            else {
                ai = path[i-1];
                bi = path[i];
                ci = path[i+1];
            }
            
            int li = metai[path[i]];
            point d, dp, xp, bap;
            
            //Add on average about "target" outliers to each path
            const double target = 0.1;
            
            if (prepareDuplicateScore(ai, bi, ci, li, d, dp, xp, bap)) continue;
            
            double tt = findDensity(dp, xp, target, li);
            
            int pi = path[i];
            int matches = mod[li]->getNear(dp, xp, bap, tt, match);//target/fac
            
            map<int, pair<double, int> > mins;
            for (int i = 0; i < matches; i++) {
                int di = match[i];
                //if (di == bi) continue; // same hit
                //if (assignment[di]!=0) continue;
                double s = evaluateScore(di, dp, xp, bap);//*fac; scoreTriple(ai,bi,di);
                //if (s > 1.0) continue;
                if (meta[di].z != meta[pi].z) {
                    int zi = metaz[di];
                    if (!mins.count(zi) || s < mins[zi].first) {
                        mins[zi] = make_pair(s, di);
                    }
                }
            }
            
            for (auto &p : mins)
                ext.push_back(p.second.second);
        }
        path.clear();
        path.shrink_to_fit();
        extended.push_back(ext);
    }
    
    // Remove duplicate hits in a path
    // by moving the hits through a set
    for (auto &path : extended) {
        set<int> s;
        for(auto &it : path) s.insert(it);
        path.assign( s.begin(), s.end() );
        path.shrink_to_fit();
    }
    
    // Remove duplicate paths
    // by sorting and pruning
    sort(extended.begin(),extended.end());
    extended.erase(unique(extended.begin(),extended.end() ),extended.end());
    
    return extended;
}


//Similar to the other prepareXScore functions, but now we try to find duplicates on layer "li" to one of the input hits. This means looking for hits that are close to the helix fitted through "ai", "bi", "ci"
int Reconstruction::prepareDuplicateScore(int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&dirp) {
    Layer&l = layer[li];
    point a = hits[ai], b = hits[bi], c = hits[ci];
    //Find circle with center p, radius r, going through a, b, and c (in xy plane)
    double ax = a.x-c.x, ay = a.y-c.y, bx = b.x-c.x, by = b.y-c.y;
    double aa = ax*ax+ay*ay, bb = bx*bx+by*by;
    double idet = .5/(ax*by-ay*bx);
    point p;
    p.x = (aa*by-bb*ay)*idet;
    p.y = (ax*bb-bx*aa)*idet;
    p.z = 0;
    double r = dist(p.x, p.y), ir = 1./r;
    p.x += c.x;
    p.y += c.y;
    
    int di = -1;
    if (metai[ai] == li) di = ai;
    else if (metai[bi] == li) di = bi;
    else if (metai[ci] == li) di = ci;
    else {
        cout << "prepareDuplocateScore given layeri not corresponding to ai, bi or ci" << endl;
        return -1;
    }
    d = hits[di];
    dp = polar[di];
    
    double rx = hits[di].x-p.x, ry = hits[di].y-p.y;
    
    //TODO: do with respect to nearest circle arc, not ca
    point ca = hits[ci]-hits[ai];
    double ang_ca = asin(dist(ca.x, ca.y)*.5*ir)*2;
    double cross = rx*ca.y-ry*ca.x;
    
    point dir;
    if (ir) {
        dir.x =-ry*ang_ca;
        dir.y = rx*ang_ca;
        dir.z = ca.z;
        if (cross < 0) dir.z *= -1;
    } else {
        dir = ca;
    }
    
    xp = point(0,0,0);
    dirp = point(d.x*dir.x+d.y*dir.y, d.x*dir.y-d.y*dir.x, dir.z*dp.x);
    //cout << dirp << endl; //dirp.x = l.avgr
    if (l.type == Tube)
        dirp = dirp.x ? dirp*(1./dirp.x) : point(1,0,0);
    else
        dirp = dirp.z ? dirp*(1./dirp.z) : point(0,0,1);
    
    return 0;
}


//Find density by binary search
//This means we want to find (and return) tt such that getDensity3(dp, xp, tt, li) = target
double Reconstruction::findDensity(point&dp, point&xp, double target, int li) {
    double Ad = 0, A = 0, B = 1, Bd;
    while (1) {
        Bd = getDensity3(dp, xp, B, li);
        //cout << B << ' ' << Bd << endl;
        if (B > 1e20) {
            cout << "No density?" << endl;
            //cout << dp << ' ' << xp << ' ' << li << endl;
            //exit(0);
            return 1e20;
        }
        if (Bd > target) break;
        B *= 10;
        if (target/Bd < 1e8) B = max(B, target/Bd);
    }
    double mid = B/2;
    int cc = 0;
    while (1) {
        double density = getDensity3(dp, xp, mid, li);
        if (density > target) {
            B = mid;
            Bd = density;
        }
        else {
            A = mid;
            Ad = density;
        }
        
        //cout << A << ' ' << mid << ' ' << B << ' ' << density << endl;
        if ((B-A) < A*1e-3 || (density > target*0.9 && density < target*1.1) || cc >= 100) break;
        mid = max(A*0.9+B*0.1, min(B*0.9+A*0.1, (target-Ad)*(B-A)/(Bd-Ad)+A));
        if (++cc == 100) { //Should never happen
            cout << "Warning: Infinite loop in findDensity" << endl;
            /*cout << dp << endl;
             cout << xp << endl;
             cout << mid << endl;
             cout << target << endl;
             cout << li << endl;
             exit(0);*/
        }
    }
    return mid;
}


//Prune away paths based on path score, also take away the start of the path if it doesn't fit well. Also remove duplicate paths using hashing
vector<vector<int> > Reconstruction::prunePaths(vector<vector<int> >&paths) {
    double const static thres = 1.8;
    double const static density = 1.0;
    
    vector<vector<int> > pruned;
    vector<vector<int> > hash_list;
    hash_list.resize(1<<17);
    
    for (auto &path : paths) {
        while (path.size() >= 4 && scoreQuadrupleDensity(path[3], path[2], path[1], path[0]) > thres) {
            path.erase(path.begin());
        }
        
        double s = -scorepathDensity(path);
        if (s < density && path.size() >= 3) {
            int h = 0;
            for (int i : path) h ^= i;
            h &= (1<<17)-1;
            int found = 0;
            for (int j : hash_list[h]) {
                if (pruned[j] == path) {
                    found = 1;
                    break;
                }
            }
            if (!found) {
                hash_list[h].push_back((int)pruned.size());
                path.shrink_to_fit();
                pruned.push_back(path);
            }
        }
        path.clear();
        path.shrink_to_fit();
    }
    
    pruned.shrink_to_fit();
    return pruned;
}


//Data structure to extract index of path with lowers score and dynamically updating, supports:
// - Adding new (index, score) pair
// - Updating (index, score) pair
// - Extracting pair with lowest score
//Note that std::priority_queue is vastly faster than std::set (which was previously used, and makes for easier implementation)
class myMap2 {
public:
    priority_queue<pair<double, int> > pq;
    double*b;
    int realsize;
    
    myMap2(int n) {
        b = new double[n];
        realsize = 0;
    }
    
    ~myMap2() {
        delete[]b;
    }
    
    void add(int i, double score) {
        pq.push(make_pair(score, i));
        b[i] = score;
        realsize++;
    }
    
    void update(int i, double score) {
        add(i, score);
        realsize--;
        //Change *4 to a lower constant > 1 for (slightly) lower memory usage, this is currently the bottleneck for memory I think
        if ((int)pq.size() >= realsize*4 && pq.size() > 1e6) {
            priority_queue<pair<double, int> > clean;
            pair<double, int> last(1e9, 1e9);
            while (pq.size()) {
                auto p = pq.top();
                pq.pop();
                if (b[p.second] == p.first && p != last) clean.push(p);
                last = p;
            }
            swap(pq, clean);
        }
    }
    
    pair<int, double> pop() {
        while (pq.size() && b[pq.top().second] != pq.top().first) pq.pop();
        int i = pq.top().second;
        double score = pq.top().first;
        b[i] = -1e9;
        pq.pop();
        realsize--;
        while (pq.size() && b[pq.top().second] != pq.top().first) pq.pop();
        return make_pair(i, score);
    }
    
    int notempty() {
        return !pq.empty();
    }
    
};


//Find assignment of hits to paths
//Do this by iteratively:
// 1. Take path with highest score
// 2. Assign all hits in that path to the path's index
// 3. Remove all hits in the path from all other paths
// 4. Repeat from step 1 until all paths are empty
vector<vector<int> > Reconstruction::findAssignment(vector<vector<int> >&paths, PolarModule *mod[LAYERS], int use_trash = 1) {
    paths.insert(paths.begin(), vector<int>());
    map<int, int> map_assignment;
    myMap2 path_score((int)paths.size()+(int)hits.size());
    vector<vector<pair<int, int> > > used_by;
    long n = hits.size();
    used_by.resize(n);
    
    for (int i = 1; i < (int)paths.size(); i++) {
        double score = scorepathDensity(paths[i]);
        path_score.add(i, score);
        for (int j = 0; j < (int)paths[i].size(); j++) {
            paths[i][j] = abs(paths[i][j]); // to enable multiple calls
            int hitid = paths[i][j];
            if (hitid>=n) {
                cout << i << " hitid " << hitid << " >= " << n << endl;
                continue;
            }
            used_by[hitid].push_back(make_pair(i, j));
        }
    }
    
    int total = (int) hits.size()-1;
    for (int i = 1; i < (int)hits.size(); i++) {
        used_by[i].shrink_to_fit();
        if (/*!assignment[i] &&*/ used_by[i].empty()) total--;
    }
    
    vector<vector<int> > solution_paths;
    solution_paths.push_back(vector<int>());
    
    while (path_score.notempty()) {
        pair<int, double> pop = path_score.pop();
        int i = pop.first;
        
        if (paths[i].empty()) continue;
        
        int r[2] = {};
        for (int p : paths[i]) r[p<0]++;
        
        int trash = 0;
        if (use_trash) {
            int c = 0;
            for (int hit_id : paths[i]) if (hit_id > 0) c++;
            if (c <= 2) trash = 1;// || -pop.second >= 1e-2) trash = 1;
        }
        for (int hit_id : paths[i]) {
            if (hit_id < 0 || map_assignment.count(hit_id)) {
                if (hit_id > 0)
                    cout << "ERROR: " << hit_id << endl;
                continue;
            }
            map_assignment[hit_id] = trash ? 0 : (int) solution_paths.size();
            //if (trash) lost += truth_weight[hit_id];
        }
        if (!trash)
            solution_paths.push_back(paths[i]);
        
        for (int hit_id : paths[i]) {
            if (hit_id < 0 || hit_id >=n) continue; // out of bounds
            
            vector<pair<int, int> >&used = used_by[hit_id];
            if (used.size()==0) continue;
            for (int k = 0; k < (int)used.size(); k++) {
                if (used[k].first == i) continue;
                vector<int>&pi = paths[used[k].first];
                if (pi.empty()) continue;
                //cout << assignment[pi[used[k].second]] << ' '<< c << endl;
                pi[used[k].second] *= -1;
                double score = scorepathDensity(pi);
                path_score.update(used[k].first, score);
            }
        }
        
        paths[i].clear();
        paths[i].shrink_to_fit();
    }
    
    return solution_paths;
}


// Score a full path by looking at the probability that it could happen by outliers alone
// Multiply average number of outliers gotten from
// - The first triple
// - All consecutive quadruples
// - Duplicates
// This function is very important for score, and full of tuning opportunities
double Reconstruction::scorepathDensity(vector<int>&path) {
    vector<int> u(48);
    u.resize(0);
    //cout << u.capacity() << endl;
    
    int last_metai = -1;
    for (int i : path) {
        if (i <= 0) continue;
        int mi = metai[i];
        if (mi != last_metai) {
            //if (last_metai != -1 && next_layer[last_metai][mi] < adj_thres) return -1e4;
            u.push_back(i);
            last_metai = mi;
        }
    }
    if (u.size() < 3) return -0.1;
    double prod = 1;
    double s = scoreTripleDensity(u[0], u[1], u[2]);//+1e-4;
    s = min(s, 1e4);
    prod *= s;
    
    //Hugely important tuning parameters
    const static double quad_off = 1.0;
    const static double dup_off = 11.0;
    const static double quad_max1 = 1.0;
    const static double dup_max = 1.0;
    
    for (int i = 3; i < (int)u.size(); i++) {
        double s = min(scoreQuadrupleDensity(u[i-3], u[i-2], u[i-1], u[i]), quad_max1);
        s *= min(scoreQuadrupleDensity(u[i], u[i-1], u[i-2], u[i-3]), quad_max1)*quad_off;
        prod *= s;
    }
    
    int j = 0;
    for (int i = 0; i < (int)path.size(); i++) {
        if (path[i] <= 0) continue;
        while (metai[path[i]] != metai[u[j]]) j++;
        if (path[i] == u[j]) continue;
        int k = max(j, 2);
        double s = pow(scoreDuplicateDensity(u[k-2], u[k-1], u[k], path[i]), 0.90)*dup_off; // Really important magic constant of 0.92, I have no clue why
        s = min(s, dup_max);
        prod *= s;
    }
    
    return -prod;
}


//Return this if no points were found, somewhat tunable parameter
const double density_eps = 1e-6;


//How many outliers do we expect to fit better than "di" to the triple "ai", "bi", "ci"
double Reconstruction::scoreDuplicateDensity(int ai, int bi, int ci, int di) {
    //cout << metai[ai] << ' ' << metai[bi] << ' ' << metai[ci] << ' ' << metai[di] << endl;
    point d, dp, xp, bap;
    if (prepareDuplicateScore(ai, bi, ci, metai[di], d, dp, xp, bap)) return 1e9;
    double s = evaluateScore(di, dp, xp, bap);
    //cout << "S0: " << s << endl;
    s = getDensity3(dp, xp, s, metai[di]);
    //if (!s) cout << "What?" << endl;
    return s+density_eps;
}


//How many outliers do we expect to fit better than "ci" in the triple "ai", "bi", "ci"?
double Reconstruction::scoreTripleDensity(int ai, int bi, int ci) {
    point d, dp, xp, bap;
    if (prepareTripleScore(ai, bi, metai[ci], d, dp, xp, bap, polar[ci])) return 1e9;
    double s = evaluateScore(ci, dp, xp, bap);
    s = getDensity3(dp, xp, s, metai[ci]);
    return s+density_eps;
}


//How many outliers do we expect to fit better than "di" in the triple "ai", "bi", "ci", "di"?
double Reconstruction::scoreQuadrupleDensity(int ai, int bi, int ci, int di) {
    point d, dp, xp, bap;
    if (prepareQuadrupleScore(ai, bi, ci, metai[di], d, dp, xp, bap, polar[di])) return 1e9;
    double s = evaluateScore(di, dp, xp, bap);
    //cout << "S0: " << s << endl;
    s = getDensity3(dp, xp, s, metai[di]);
    //if (!s) cout << "What?" << endl;
    return s+density_eps;
}


//Prepare ellipse equation of collision between line extrapolated through hits with id "ai" and "bi" and layer "li". Return collision coordinate "d", in polar coordinates "dp", ellipse stretching "xp", and direction of hit in polar coordnates "bap". "target" describes the layer, possibly corrected for a single point we are evaluating a helix quadruple
int Reconstruction::prepareTripleScore(int ai, int bi, int li, point&d, point&dp, point&xp, point&bap, point target) {
    const double slack = 1.00; //No slack
    
    Layer&l = layer[li];
    point&a = hits[ai], &b = hits[bi];
    
    point ba = b-a;
    if (l.type == Tube) {
        double vv = ba.x*ba.x+ba.y*ba.y;
        double pv = ba.x*a.x+ba.y*a.y;
        double pp = a.x*a.x+a.y*a.y;
        double RR = target.x*target.x;
        double sq = pv*pv-vv*(pp-RR);
        if (sq < 0) return -1;
        
        double t = (-pv+sqrt(sq))/vv;
        if (t < 0 || !vv) return -1;
        d.x = a.x+ba.x*t;
        d.y = a.y+ba.y*t;
        d.z = a.z+ba.z*t;
        
        if (d.z < l.minz*slack || d.z > l.maxz*slack) return -1;
        
        dp = point(dist(d.x,d.y),atan2(d.y,d.x),d.z);
        
        xp = point(0, -dp.x*(ba.x*ba.x+ba.y*ba.y), ba.z);
        bap = point(ba.x*d.x+ba.y*d.y, d.x*ba.y-d.y*ba.x, ba.z*dp.x);
        
        bap = bap*(1./bap.x);
    } else if (l.type == Disc) {
        double t = (target.z-a.z)/ba.z;
        if (t < 0 || !ba.z) return -1;
        d.x = a.x+ba.x*t;
        d.y = a.y+ba.y*t;
        d.z = a.z+ba.z*t;
        
        dp = point(dist(d.x,d.y),atan2(d.y,d.x),d.z);
        
        if (dp.x < l.minr*(1./slack) || dp.x > l.maxr*slack) return -1;
        
        xp = point(ba.x*d.y-ba.y*d.x, d.x*ba.x+d.y*ba.y, 0);
        bap = point(ba.x*d.x+ba.y*d.y, d.x*ba.y-d.y*ba.x, ba.z*dp.x);
        
        bap = bap*(1./bap.z);
    }
    double xp2 = xp.x*xp.x+xp.y*xp.y+xp.z*xp.z;
    if (xp2)
        xp = xp*sqrt((1-stretch)/xp2);
    return 0;
}


//Default is using average position in the layer
int Reconstruction::prepareTripleScore(int ai, int bi, int li, point&d, point&dp, point&xp, point&bap) {
    Layer&l = layer[li];
    point target(l.avgr, 0, l.avgz);
    return prepareTripleScore(ai, bi, li, d, dp, xp, bap, target);
}


//Get some features for the logistic regression model
bool Reconstruction::getFeatures(int ai, int bi,  float* feature, point v) {
    point &a = hits[ai], &b = hits[bi];
    point d = a-b;
    
    //  double dr2 = dist(d.x, d.y);
    double r1 = dist(a.x, a.y);
    double r2 = dist(b.x, b.y);
    double dr = r2-r1;
    
    feature[0] = dir_miss(ai, bi);//Cell's data of ai
    feature[1] = dir_miss(bi, ai);//Cell's data of bi
    //Different distances from origin (like how far does the line through ai-bi pass from the origin)
    feature[2] = wdist(a, d, 0);
    feature[3] = zdist2(a, b);
    feature[4] = wdistr(r1, dr, a.z, d.z, 1);
    feature[5] = wdist(a, d, 1);
    hits[0] = v;
    feature[6] = scoreTriple(0,ai,bi);
    
    for (int i=0;i<7;i++) {
        if (isnan(feature[i])) {
            cout << "NAN " << i << endl;
            return false;
        }
        if (feature[i]==0.0) {
            return false;
        }
    }
    
    return true;
}


//Angle between line through hits ai-bi and cell's data direction of at hit id ai
double Reconstruction::dir_miss(int ai, int bi) {
    point d = hits[ai]-hits[bi];
    double r = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
    if (r==0.0) return M_PI;
    return acos(max(fabs(hit_dir[ai][0]*d),
                    fabs(hit_dir[ai][1]*d))/r);
}


//Different distances from origin (like how far does the line through ai-bi pass from the origin)
double Reconstruction::wdistr(double r1, double dr, double az, double dz, double w) {
    double pp = r1*r1+az*az*w;
    double pd = r1*dr+az*dz*w;
    double dd = dr*dr+dz*dz*w;
    if (dd==0) return 1.E3;
    return sqrt(pp-pd*pd/dd);
}


double Reconstruction::wdist(point&a, point&d, double w) {
    double pp = a.x*a.x+a.y*a.y+a.z*a.z*w;
    double pd = a.x*d.x+a.y*d.y+a.z*d.z*w;
    double dd = d.x*d.x+d.y*d.y+d.z*d.z*w;
    if (dd==0) return 1.E3;
    return sqrt(pp-pd*pd/dd);
}


double Reconstruction::zdist(point&a, point&b) {
    static point origin(0,0,0);
    point p;
    double r;
    circle(origin, a, b, p, r);
    if (r<=0) return 1.E3;
    double ang_ab = 2*asin(dist(a.x-b.x, a.y-b.y)*.5/r);
    double ang_a = 2*asin(dist(a.x, a.y)*.5/r);
    if (ang_ab==0) return 1.E3;
    return fabs(a.z-(b.z-a.z)*ang_a/ang_ab);
}


double Reconstruction::zdist2(point&a, point&b) {
    static point origin(0,0,0);
    point p;
    double r;
    circle(origin, a, b, p, r);
    if (r<=0) return 1.E3;
    double ang_ab = 2*asin(dist(a.x-b.x, a.y-b.y)*.5/r);
    double ang_a = 2*asin(dist(a.x, a.y)*.5/r);
    if (ang_ab==0) return 1.E3;
    return fabs(b.z-a.z-a.z*ang_ab/ang_a);
}
