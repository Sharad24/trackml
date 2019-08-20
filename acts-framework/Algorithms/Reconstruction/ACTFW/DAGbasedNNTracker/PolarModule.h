#ifndef _POLARMOD_H_
#define _POLARMOD_H_
// Neural Network based tracker
// The code has been adopted from Johan Sokrates Wind's award winning trackml Kaggle contribution
// M.Kunze, Heidelberg University, 2018

class Tracker;

//Levels of detail
static const int lods = 8;

class PolarModuleInternal {
public:
    Tracker *tracker = NULL;
    int *mem = NULL;
    int *ind[lods];
    int *num[lods];
    int layeri, zi, aspectx, aspecty;
    PolarModuleInternal(int li, int zi, Tracker *t);
    void recIndex(int&i, int ix, int iy, int lod);
    ~PolarModuleInternal();
    inline double calcX(point&p);
    inline double calcX(double r);
    int getNear(point&dp0, point&xp, point&bap, double tt, int*match);
};

class PolarModule {
    PolarModuleInternal *internal;
    int internals;
public:
    PolarModule() {internal = NULL;internals = 0;}
    PolarModule(int li, Tracker *t);
    int getNear(point&dp, point&xp, point&bap, double tt, int*match);
    ~PolarModule() {
        if (internal) {
            delete[]internal;
            internal = NULL;
        }
    }
};

#endif

