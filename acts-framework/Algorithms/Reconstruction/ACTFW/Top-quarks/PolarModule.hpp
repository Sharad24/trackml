#pragma once

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/AlgorithmContext.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "main.hpp"
#include "Point.hpp"
#include "Triple.hpp"
#include "Layer.hpp"
#include <array>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cmath>
#include <array>
namespace FW {

static const int lods = 8;
// static const std::array<Layer, 48>* layer = nullptr;

class Topquarks;

class PolarModuleInternal {
public:
    static const int Tube = 0, Disc = 1;
    // AlgorithmContext ctx
    // std::array<int, 200000> match;
    // const std::vector<int>* metai = nullptr;
    // const std::array<std::array<double, 4>, 48>* disc_z = nullptr;
    // const std::vector<int>* metaz = nullptr;
    // const std::vector<point>* hits = nullptr;
    // const std::array<int, 200000>* assignment = nullptr;
    // const std::vector<point>* polar = nullptr;
    int *mem = nullptr;
    // std::array<int, lods>* mem;
    // int *ind[lods];
    std::array<std::vector<int>, lods> num;
    std::array<std::vector<int>, lods> ind;
    // std::vector<int> data;
    // int *num[lods];
    int layeri, zi, aspectx, aspecty;
    PolarModuleInternal(int li, AlgorithmContext ctx, int zi);
    void recIndex(AlgorithmContext ctx, int&i, std::array<std::vector<int>, lods>&ind, std::array<std::vector<int>, lods>&num, int ix, int iy, int lod);
    ~PolarModuleInternal();
    inline double calcX(AlgorithmContext ctx, point&p, int li, int zi) const;
    inline double calcX(AlgorithmContext ctx, double r, int li, int zi) const;
    int getNear(AlgorithmContext ctx, point&dp0, point&xp, point&bap, double tt, std::array<int, 200000>& match, int li, int zi, int pointer) const;
};

class PolarModule {
    // PolarModuleInternal *internal;
    std::vector<PolarModuleInternal> internal_data;
    std::vector<PolarModuleInternal>* internal;
    int layeri;
    int internals;
public:

    // PolarModule() {internal = nullptr ;internals = 0;}
    PolarModule() {internals = 0; layeri=0;}// internal_data{{{0, FW::AlgorithmContext ctx, 0}, {0, ctx, 1}, {0, ctx, 2}, {0, ctx, 3}}};}
    PolarModule(int li, AlgorithmContext ctx);
    int getNear(AlgorithmContext ctx, point&dp, point&xp, point&bap, double tt, std::array<int, 200000>& match, int li) const;
    // ~PolarModule() {
    //     if (internal) {
    //         delete[]internal;
    //         internal = NULL;
    //     }
    // }
    ~PolarModule();
    // std::array<PolarModuleInternal, 4> internal_data;

};
}
