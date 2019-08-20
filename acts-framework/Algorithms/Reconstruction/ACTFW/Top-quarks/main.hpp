// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


/// @author Sharad Chitlangia <sharad.chitlangia@cern.ch>

#pragma once

#include <stdexcept>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cmath>
#include <array>
#include <queue>

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "Point.hpp"
#include "PolarModule.hpp"
#include "Triple.hpp"
namespace FW {

/// @brief An empty reconstruction algorithm to be use as a template.

// static const int lods = 8;
static const int crude_steps = 1<<10;
static int done = 0;
// extern std::array<PolarModule, 48> mod;

// class PolarModule;

class Topquarks : public BareAlgorithm
{
public:
  struct Config
  {
    /// Input space point collection
    std::string spacePointCollection;
    /// Some config option that changes the reconstruction
    double aConfigOption = 1.23;
  };

  /*struct point {
    double x, y, z;
    point() {}
    point(double x, double y, double z) : x(x),y(y),z(z) {}
    inline point operator-(const point&p) {
      return point(x-p.x, y-p.y, z-p.z);
    }
    inline point operator+(const point&p) {
      return point(x+p.x, y+p.y, z+p.z);
    }
    inline double operator*(const point&p) {
      return x*p.x+y*p.y+z*p.z;
    }
    inline point operator*(double f) {
      return point(x*f, y*f, z*f);
    }
  };*/

  struct Layer {
    double minr, avgr, maxr;
    double minz, avgz, maxz;
    int count;
    int type;

    double var0, var1;
  };

  struct Detector {
    int volume_id, layer_id, module_id;
    point c;
    point rx, ry, rz;
    double d, minw, maxw, h, cell_w, cell_h;
  };

  struct myMap2 {
    std::priority_queue<std::pair<double, int> > pq;
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
      pq.push(std::make_pair(score, i));
      b[i] = score;
      realsize++;
    }
    void update(int i, double score) {
      add(i, score);
      realsize--;
      //Change *4 to a lower constant > 1 for (slightly) lower memory usage, this is currently the bottleneck for memory I think
      if (pq.size() >= realsize*4 && pq.size() > 1e6) {
        std::priority_queue<std::pair<double, int> > clean;
        std::pair<double, int> last(1e9, 1e9);
        while (pq.size()) {
  	auto p = pq.top();
  	pq.pop();
  	if (b[p.second] == p.first && p != last) clean.push(p);
  	last = p;
        }
        std::swap(pq, clean);
        //cout << pq.size() << ' ' << realsize << endl;
      }
      //static int cc = 0;
      //if (++cc%100000 == 0)
        //  cout << pq.size() * 1. / realsize << endl;
    }
    std::pair<int, double> pop() {
      while (pq.size() && b[pq.top().second] != pq.top().first) pq.pop();
      int i = pq.top().second;
      double score = pq.top().first;
      b[i] = -1e9;
      pq.pop();
      realsize--;
      while (pq.size() && b[pq.top().second] != pq.top().first) pq.pop();
      return std::make_pair(i, score);
    }
    int notempty() {
      return !pq.empty();
    }
  };

  /// Construct the algorithm.
  ///
  /// @param [in] cfg is the configuration struct
  /// @param [in] loglevel is the logging level
  Topquarks(const Config&        cfg,
                               Acts::Logging::Level logLevel);

  /// Execute the algorithm and reconstruct nothing.
  ///
  /// @param [in] ctx is the algorithm context for event consistency
  /// @return is a process code indicating success or not
  FW::ProcessCode
  execute(AlgorithmContext ctx) const final override;

  void
  readEventData(AlgorithmContext ctx) const;

  void
  readDetectorInfo(AlgorithmContext ctx) const;

  // std::vector<std::pair<int, int> >
  void
  findPairs(AlgorithmContext ctx, int modeli) const;

  std::vector<triple>
  findTriples(AlgorithmContext ctx, std::array<PolarModule, 48> &mod, std::array<int, 200000>&match, int method, double target) const;

  // FW::ProcessCode
  // prunePaths(AlgorithmContext ctx);
  //
  // FW::ProcessCode
  // addDuplicates(AlgorithmContext ctx);

  FW::ProcessCode
  hitAssignment(AlgorithmContext ctx);

  FW::ProcessCode
  writeSubmission(AlgorithmContext ctx);

  void
  readTubes(AlgorithmContext ctx) const;

  bool
  z_cmp(int &a, int&b, AlgorithmContext ctx) const;

  bool
  r_cmp(int&a, int&b, AlgorithmContext ctx) const;

  double
  wdist(point&a, point&d, double w) const;

  double
  wdistr(double r1, double dr, double az, double dz, double w) const;

  double
  zdist(point&a, point&b);

  double
  zdist2(point&a, point&b) const;

  double
  dist(double x, double y) const;

  double
  dist2(double x, double y) const;

  double dist0(point&p) const;

  void
  circle(const point&a, const point&b, const point&c, point&p, double&r) const;

  double
  dir_miss(AlgorithmContext ctx, int ai, int bi) const;

  point
  normalize(point &a) const;

  double
  findDensity(AlgorithmContext ctx, point&dp, point&xp, double target, int li) const;

  double
  getDensity3(AlgorithmContext ctx, point&dp, point&xp, double tt, int li) const;

  int
  getIndex(FW::AlgorithmContext ctx, int&li, double x) const;

  double
  extendTripleLine(AlgorithmContext ctx, std::vector<triple>&triples, int ai, int bi, const point&a, const point&b, int li, PolarModule&mod, std::array<int, 200000> &match, int rev) const;

  double
  extendTripleOrigin(AlgorithmContext ctx, std::vector<triple>&triples, int ai, int bi, const point&a, const point&b, int li, PolarModule&mod, std::array<int, 200000> &match, double target , int rev) const;

  double
  scoreTriple(AlgorithmContext ctx, int ai, int bi, int ci) const;

  int
  acceptTriple(AlgorithmContext ctx, triple&t) const;

  double
  scoreTripleDensity(AlgorithmContext ctx, int ai, int bi, int ci) const;

  int
  prepareTripleScore(AlgorithmContext ctx, int ai, int bi, int li, point&d, point&dp, point&xp, point&bap, point target) const;

  int
  prepareTripleScoreDefault(AlgorithmContext ctx, int ai, int bi, int li, point&d, point&dp, point&xp, point&bap) const;

  int
  prepareQuadrupleScore(AlgorithmContext ctx, int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&dirp, point target, double sign) const;

  double
  field(double z) const;

  double
  evaluateScore(AlgorithmContext ctx, int ci, point&dp, point&xp, point&bap) const;

  int
  prepareQuadrupleScoreDefault(AlgorithmContext ctx, int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&bap, double sign) const;

  double
  scoreTripleLogRadius_and_HitDir(AlgorithmContext ctx, int ai, int bi, int ci, double*L) const;

  void
  addDuplicateTriples(AlgorithmContext ctx, std::vector<triple>&triples, std::array<PolarModule, 48> &mod, std::array<int, 200000> &match, int ai, int bi) const;

  point
  topolar(point&dir, point&ref, point&refp) const;

  std::vector<std::pair<int, int> >
  findDuplicates(AlgorithmContext ctx, std::array<PolarModule, 48> &mod, std::array<int, 200000> &match) const;

  std::vector<triple>
  findTriplesDuplicates(AlgorithmContext ctx, std::array<PolarModule, 48> &mod, std::array<int, 200000> &match) const;

  std::vector<triple>
  pruneTriples(AlgorithmContext ctx, std::vector<triple> &triples) const;

  point
  scoreTripleHitDir(triple t) const;

  void
  storeInEventStore(AlgorithmContext ctx) const;

  std::vector<std::vector<int> >
  findPaths(FW::AlgorithmContext ctx, std::vector<triple>&triples, std::array<PolarModule, 48> *mod, std::array<int, 200000> &match) const;

  int
  extend3(AlgorithmContext ctx, std::array<PolarModule, 48> *mod, std::array<int, 200000> &match, int ai, int bi, int ci, int li, double target) const;

  std::vector<std::vector<int> >
  addDuplicates(AlgorithmContext ctx, std::vector<std::vector<int> >&paths, std::array<PolarModule, 48> &mod, std::array<int, 200000> &match) const;

  int
  prepareDuplicateScore(AlgorithmContext ctx, int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&dirp, std::array<int, 200000> &match) const;

  std::vector<std::vector<int> >
  prunePaths(AlgorithmContext ctx, std::vector<std::vector<int> >&paths, std::array<int, 200000> &match) const;

  double
  scoreQuadrupleDensity(AlgorithmContext ctx, int ai, int bi, int ci, int di) const;

  double
  scorepathDensity(AlgorithmContext ctx, std::vector<int>&path, std::array<int, 200000> &match) const;

  double
  scoreDuplicateDensity(AlgorithmContext ctx, int ai, int bi, int ci, int di, std::array<int, 200000> &match) const;

  std::vector<std::vector<int> >
  findAssignment(AlgorithmContext ctx, std::vector<std::vector<int> >&paths, std::array<PolarModule, 48> *mod, std::array<int, 200000> &match, int use_trash = 1) const;

  void
  initDensity3(AlgorithmContext ctx) const;

  std::vector<std::vector<int> >
  extendPaths(AlgorithmContext ctx, std::vector<std::vector<int> >&paths, std::array<PolarModule, 48> *mod, std::array<int, 200000>& match) const;

  void
  loadAdjacent(AlgorithmContext ctx) const;

  void
  writeSubmission(AlgorithmContext ctx, std::map<int, int>&assignment) const;

  // static const int Tube = 0, Disc = 1;
  // static const int adj_thres = 1000;
  static constexpr double density_eps = 0.000001;
  static constexpr double stretch = 0.02;

private:
  Config m_cfg;

protected:
  static const int Tube = 0, Disc = 1;
  std::vector<int>* check;
  std::array<std::vector<int>, 48>* tube;
  std::array<Layer, 48>* layer;
  std::array<std::array<point, 2>, 200000>* hit_dir;
  int adj_thres = 1000;
  std::array<int, 48>* itopo;
  std::array<int, 48>* topo;
  int next_layer[48][48];
  std::array<std::vector<double>, 48>* sorted_hits;
  std::vector<point>* hits;
  std::array<std::vector<std::pair<std::pair<int, int>, double> >, 200000>* hit_cells;
  std::vector<std::pair<int, int> > origin_pairs, *pairs;
  std::vector<triple> *triples;
  std::array<int, 200000>* match; //reused global array for storing matching hits from PolarModule->getNear function
  std::array<std::array<int, 48>, crude_steps>* crudeIndex;
  std::array<std::pair<double, double>, 48>* crudeIndex_a;
  std::array<std::array<std::array<double, 3>, 20000>, 48>* poly;
  std::vector<point>* polar;
  std::vector<point>* meta;
  std::vector<int>* metai, *metaz;
  std::map<std::pair<int, int>, double> *v_duplicates;
  std::array<std::array<double, 4>, 48> *disc_z;
  std::vector<std::vector<int> > paths;
  // std::vector<std::pair<double, int> >* v;
  friend class PolarModule;
  // std::array<PolarModule, 48> *mod;
  friend class PolarModuleInternal;
};

}  // namespace FW
