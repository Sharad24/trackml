// Neural Network based tracker
// The code has been adopted from Johan Sokrates Wind's award winning trackml Kaggle contribution
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"
#include "PolarModule.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

//Wrapper of PolarModuleInternal to keep 4 separate modules, one for each z value of disc layers

//Acceleration data structure to quickly query all hits close to a helix intersection
//Some features are:
// - Quadtree structure (with fully allocated layers)
// - Logarithmic scaling in radial dimension / aspect ratios to keep cells square
// - Bounding the elliptic cylinder intersection to guarantee all appropriate hits are found
PolarModuleInternal::PolarModuleInternal(int li, int zi = 0, Tracker *t = NULL) : tracker(t), zi(zi)  {
    layeri = li;
    Layer &l = t->layer[layeri];
    if (l.type == t->Disc) {
        aspectx = 1;
        aspecty = 2*M_PI/log(l.maxr/l.minr)+.5;
    } else {
        double a = 2*M_PI*l.avgr/(l.maxz-l.minz);
        aspectx = aspecty = 1;
        if (a > 1)
            aspecty = ceil(a);
        else
            aspectx = ceil(1./a);
    }
    for (int lod = 0; lod < lods; lod++) {
        ind[lod] = new int[aspectx*aspecty<<lod*2]{};
        num[lod] = new int[aspectx*aspecty<<lod*2]{};
    }
    for (int i = 1; i < (int)tracker->hits.size(); i++) {
        if (tracker->assignment[i]) continue;
        if (t->metai[i] == layeri && (l.type == tracker->Tube || t->metaz[i] == zi)) {
            double x = calcX(tracker->polar[i]), y = tracker->polar[i].y*.5/M_PI+.5;
            if (!(x >= 0 && y >= 0 && x < 1 && y < 1)) {
                std::cout << "ERROR: hit " << i << " outside of bounds of detector " << layeri << std::endl;
                std::cout << x << ' ' << y << std::endl;
                continue;
            }
            for (int lod = 0; lod < lods; lod++) {
                int ix = x*(aspectx<<lod);
                int iy = y*(aspecty<<lod);
                num[lod][ix+iy*(aspectx<<lod)]++;
            }
        }
    }
    int tot = 0;
    for (int j = 0; j < aspecty; j++)
        for (int i = 0; i < aspectx; i++)
            tot += num[0][i+aspectx*j];
    int k = 0;
    mem = new int[tot];
    for (int i = 0; i < aspectx; i++)
        for (int j = 0; j < aspecty; j++)
            recIndex(k, i, j, 0);
    for (int i = 1; i < (int)tracker->hits.size(); i++) {
        if (tracker->assignment[i]) continue;
        if (t->metai[i] == layeri && (l.type == t->Tube || t->metaz[i] == zi)) {
            double x = calcX(tracker->polar[i]), y = tracker->polar[i].y*.5/M_PI+.5;
            int ix = x*(aspectx<<(lods-1));
            int iy = y*(aspecty<<(lods-1));
            mem[--ind[lods-1][ix+iy*(aspectx<<(lods-1))]] = i;
        }
    }
    
    /*if (l.type == Tube)
     cout << "Tube: " << l.avgr*2*M_PI << ' ' << l.maxz-l.minz << endl;
     else
     cout << "Disc: " << l.minr << ' ' << l.maxr << ' ' << M_PI*2*l.maxr << endl;*/
}

void PolarModuleInternal::recIndex(int&i, int ix = 0, int iy = 0, int lod = 0) {
    if (lod == lods-1) {
        i += num[lod][ix+iy*(aspectx<<lod)];
        ind[lod][ix+iy*(aspectx<<lod)] = i;
    } else {
        ind[lod][ix+iy*(aspectx<<lod)] = i;
        for (int y = 0; y < 2; y++)
            for (int x = 0; x < 2; x++)
                recIndex(i, ix*2+x, iy*2+y, lod+1);
    }
}

PolarModuleInternal::~PolarModuleInternal() {
    if (mem) {
        for (int lod = 0; lod < lods; lod++) {
            delete[]ind[lod];
            delete[]num[lod];
        }
        delete[]mem;
        mem = NULL;
    }
}

inline double PolarModuleInternal::calcX(point&p) {
    Layer&l = tracker->layer[layeri];
    if (l.type == tracker->Disc) //Logarithmic scaling for polar coordinates
        return p.x > l.minr ? log(p.x/l.minr)/log(l.maxr/l.minr) : 0.;//(p.x-l.minr)/(l.maxr-l.minr);
    else
        return (p.z-l.minz)/(l.maxz-l.minz);
}

inline double PolarModuleInternal::calcX(double r) {
    Layer&l = tracker->layer[layeri];
    if (l.type == tracker->Disc)
        return r > l.minr ? log(r/l.minr)/log(l.maxr/l.minr) : 0.;//(r-l.minr)/(l.maxr-l.minr);
    else
        return (r-l.minz)/(l.maxz-l.minz);
}

//Find all hits di with "evaluateScore(di, dp, xp, bap) < tt" in array match, return number of matches
int PolarModuleInternal::getNear(point&dp0, point&xp, point&bap, double tt, int*match) {
    Layer&l = tracker->layer[layeri];
    
    //cout << layeri << ' ' << l.minr << ' ' << l.maxr << endl;
    
    double yscale = 1./(dp0.x*M_PI*2);
    
    double ext_xm = 0, ext_xp = 0, ext_ym = 0, ext_yp = 0;
    point dp = dp0;
    if (l.type == Tracker::Disc) {
        point off = bap*(tracker->disc_z[layeri][zi]-dp0.z);
        dp.x = dp0.x+off.x;
        dp.y = dp0.y+off.y/dp0.x;
        dp.z = dp0.z+off.z;
    } else {
        double out = l.maxr-dp.x;
        double in  = l.minr-dp.x;
        ext_xm = std::max(0.,-std::min(out*bap.z, in*bap.z));
        ext_xp = std::max(0., std::max(out*bap.z, in*bap.z));
        ext_ym = std::max(0.,-std::min(out*bap.y, in*bap.y)*yscale);
        ext_yp = std::max(0., std::max(out*bap.y, in*bap.y)*yscale);
        //cout << ext_xm << ' ' << ext_xp << ' ' << ext_ym << ' ' << ext_yp << endl;
    }
    double x, y = dp.y*.5/M_PI+.5;
    double dx, dy = xp.y;
    if (l.type == Tracker::Disc) {
        x = dp.x;
        dx = xp.x;
    } else {
        x = dp.z;
        dx = xp.z;
    }
    
    int lod = -log2(tt*pow(yscale*aspecty, 2))*.5;
    if (lod > lods-1) lod = lods-1;//{cout << "Too large lod: " << lod << endl; lod = lods-1;}
    if (lod < 0) lod = 0;
    int stride = aspectx<<lod;
    
    double inv = 1./(1-dx*dx-dy*dy);
    double rx = sqrt(tt*(1-dy*dy)*inv);
    double ry = yscale*sqrt(tt*(1-dx*dx)*inv);
    if (rx != rx || ry != ry || x != x || y != y) {
        /*            if (DEBUG) {
         std::cout << "NaN in PolarModule" << std::endl;
         std::cout << layeri << ' ' << dp.x << ' ' << x << ' ' << y << ' ' << rx << ' ' << ry << std::endl;
         exit(0);
         } else
         */
        return 0;
    }
    
    int sx =  calcX(x-rx-ext_xm)*(aspectx<<lod);
    int sy = floor((y-ry-ext_ym)*(aspecty<<lod));
    int ex =  calcX(x+rx+ext_xp)*(aspectx<<lod)+1;
    int ey = floor((y+ry+ext_yp)*(aspecty<<lod))+1;
    
    //cout << sx << ' ' << ex << ' ' << sy << ' ' << ey << ' ' << (1<<lod) << endl;
    sx = std::max(sx, 0);
    ex = std::min(ex, stride);
    
    ey = std::min(ey, sy+(aspecty<<lod));
    
    //cout << (ex-sx)*(ey-sy) << endl;
    int considered = 0;
    int matches = 0;
    
    int ii = sy;
    int i = ii%(aspecty<<lod);
    if (i < 0) i += aspecty<<lod;
    do {
        for (int j = sx; j < ex; j++) {
            int&n = num[lod][i*stride+j];
            if (!n) continue;
            //if (rx*rx+ry*ry-dot*dot < t2) //TODO
            considered += n;
            int k0 = ind[lod][i*stride+j];
            for (int k = 0; k < n; k++) {
                int hit_id = mem[k0+k];
                
                point&r = tracker->polar[hit_id];
                //if (r.x < l.minr || r.x > l.maxr) cout << "What?" << endl;
                point err = r-dp0;
                if (err.y > M_PI) err.y -= M_PI*2;
                if (err.y <-M_PI) err.y += M_PI*2;
                err.y *= dp0.x; //Why doesn't dp.x work better than r.x?
                
                err = err-bap*(l.type == Tracker::Disc ? err.z : err.x);
                double dot = err*xp;
                double r2 = err*err-dot*dot;
                if (r2 < tt) match[matches++] = hit_id;
            }
        }
        ii++;
        if (++i == (aspecty<<lod)) i = 0;
    } while (ii < ey);
    //cout << considered << endl;
    //cout << endl;
    return matches;
}


PolarModule::PolarModule(int li, Tracker *t) {
    Layer&l = t->layer[li];
    if (l.type == t->Disc) {
        internals = 4;
        internal = new PolarModuleInternal[4]{{li,0,t}, {li,1,t}, {li,2,t}, {li,3,t}};
    } else {
        internals = 1;
        internal = new PolarModuleInternal[1]{{li,0,t}};
    }
}

int PolarModule::getNear(point&dp, point&xp, point&bap, double tt, int*match) {
    int matches = 0;
    for (int i = 0; i < internals; i++) {
        matches += internal[i].getNear(dp, xp, bap, tt, match+matches);
    }
    return matches;
}
