// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "main.hpp"

#include <stdexcept>
#include <cmath>


#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/GeometryID.hpp>

#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

// extern std::array<PolarModule, 48> mod;

FW::Topquarks::Topquarks(
    const FW::Topquarks::Config& cfg,
    Acts::Logging::Level                            logLevel)
  : BareAlgorithm("Topquarks", logLevel), m_cfg(cfg)
{
  if (m_cfg.spacePointCollection.empty()) {
    throw std::invalid_argument("Missing input space points collection");
  }
}

FW::ProcessCode
FW::Topquarks::execute(FW::AlgorithmContext ctx) const
{
  ACTS_INFO("empty reconstruction on event " << ctx.eventNumber);

  // if (FW::Topquarks::readDetectorInfo() == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Error in reading detector info");
  //   return FW::ProcessCode::ABORT;
  // }
  // ACTS_INFO("Read Detector info");
  //
  // if (FW::Topquarks::readEventData() == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Error in reading data");
  //   return FW::ProcessCode::ABORT;
  // }

  // PolarModule *mod[48];

  // int check[3] = {0, 1, 2};
  //
  // if (ctx.eventStore.add("check", std::move(check)) == FW::ProcessCode::ABORT){
  //   return FW::ProcessCode::ABORT;
  // }
  // std::array<int, 48>*topo = new std::array<int, 48>();
  // std::array<int, 48>*itopo = new std::array<int, 48>();
  readDetectorInfo(ctx);
  ACTS_INFO("Read Detector Info");

  readEventData(ctx);
  ACTS_INFO("Read Event files");
  int Tube = 0, Disc = 1;
  ACTS_INFO("Check " << Tube << Disc);

  // const std::array<std::vector<int>, 48>* tube = nullptr;

  const std::vector<FW::point>* hits = nullptr;
  const std::vector<FW::point>* polar = nullptr;
  const std::vector<FW::point>* meta = nullptr;
  const std::vector<int>* metai = nullptr;
  const std::vector<int>* metaz = nullptr;
  const std::array<Layer, 48>* layer = nullptr;

  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load hits in execute");
  }
  if (ctx.eventStore.get("polar", polar) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load polar in execute");
  }
  if (ctx.eventStore.get("meta", meta) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load meta in execute");
  }
  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load hits in execute");
  }
  if (ctx.eventStore.get("metaz", metaz) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load metaz in execute");
  }
  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load kayer in execute");
  }

  // ACTS_INFO(hits->size());
  // ACTS_INFO(metai->size());
  // ACTS_INFO(meta->size());
  // ACTS_INFO(metaz->size());
  // ACTS_INFO(polar->size());
  //ACTS_INFO(layer->size());


  std::array<PolarModule, 48> mod;
  // std::cout << "\n\n\n\n\n\n\n";
  for (int i = 0; i < 48; i++){
    ACTS_VERBOSE("Creating polarmodule " << i);
    mod[i] = PolarModule(i, ctx);
  }
  ACTS_INFO("Created PolarModule array");

  // if (ctx.eventStore.add("mod", std::move(mod)) == FW::ProcessCode::ABORT){
  //   ACTS_WARNING("Error in adding mod to eventStore in execute");
  // }
  // ACTS_INFO("Added PolarModule array to eventStore");


  // std::cout << hits->at(0).x;

  readTubes(ctx);

  // if (ctx.eventStore.add("tube", std::move(tube)) == FW::ProcessCode::ABORT){
  //   ACTS_INFO("Not able to store tube in eventStore");
  // }
  ACTS_INFO("Read and Stored Tubes");

  ACTS_VERBOSE("initializing density objects")
  initDensity3(ctx);

  // findPairs(ctx, 3);
  std::vector<std::pair<int, int> >* pairs = new std::vector<std::pair<int, int> >();
  FILE*fp6 = fopen("/Users/sharadchitlangia/Desktop/trackml/test_sub/pairs", "r");
  int n;
  int tmp = fscanf(fp6, "%d", &n);
  ACTS_INFO(n << " pairs");
  // pairs->resize(n);
  for (int i = 0; i < n; i++){
    int a, b;
    // ACTS_INFO("Read " << i);
    tmp = fscanf(fp6, "%d%d", &a, &b);
    pairs->push_back(std::make_pair(a,b));
  }
  fclose(fp6);
  if (ctx.eventStore.add("pairs", std::move(*pairs)) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to add pairs from eventStore in execute");
  }

  // const std::vector<std::pair<int, int> >* pairs = nullptr;
  //
  // if (ctx.eventStore.get("pairs", pairs) == FW::ProcessCode::ABORT){
  //   ACTS_WARNING("Not able to get pairs from eventStore in execute");
  // }
  //
  // ACTS_INFO("Found pair " << pairs->size() << " candidates");

  loadAdjacent(ctx);
  const std::array<std::pair<double, double>, 48>* crudeIndex_a = nullptr;
  const std::array<std::array<int, 48>, crude_steps>* crudeIndex = nullptr;
  const std::array<std::vector<double>, 48>* sorted_hits = nullptr;

  if (ctx.eventStore.get("crudeIndex", crudeIndex) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get crudeIndex from eventStore in getIndex");
  }
  if (ctx.eventStore.get("crudeIndex_a", crudeIndex_a) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get crudeIndex_a from eventStore in getIndex");
  }
  if (ctx.eventStore.get("sorted_hits", sorted_hits) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get sorted_hits from eventStore in getIndex");
  }

  FILE*fp14 = fopen("sorted_hits", "w");
  // fprintf(fp14, "%d\n", int(pairs->size()));
  for (int i = 0; i < 48; i++){
    for (auto a = sorted_hits->at(i).begin(); a!=sorted_hits->at(i).end(); a++){
      fprintf(fp14, "%lf\n", sorted_hits->at(i)[*a]);
    }
  }
  fclose(fp14);

  FILE*fp15 = fopen("crudeIndex_a", "w");
  // fprintf(fp14, "%d\n", int(pairs->size()));
  for (int i = 0; i < 48; i++)
    fprintf(fp15, "%lf %lf\n", crudeIndex_a->at(i).first, crudeIndex_a->at(i).second);
  fclose(fp15);

  FILE*fp16 = fopen("crudeIndex", "w");
  // fprintf(fp14, "%d\n", int(pairs->size()));
  for (int i = 0; i < crude_steps; i++){
    for (int j = 0; j < 48; j++){
      fprintf(fp16, "%d %d\n", crudeIndex->at(i)[j], crudeIndex->at(i)[j]);
    }
  }
  fclose(fp16);


  // if (FW::Topquarks::findTriples(ctx) == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Error in finding triplets");
  //   return FW::ProcessCode::ABORT;
  // }
  std::array<int, 200000> match;
  std::vector<FW::triple> triples = findTriples(ctx, mod, match, 0, 0.5);

  FILE*fp13 = fopen("triples", "w");
  fprintf(fp13, "%d\n", int(triples.size()));
  for (int i = 0; i < triples.size(); i++)
    fprintf(fp13, "%d %d %d\n", triples[i].x, triples[i].y, triples[i].z);
  fclose(fp13);
	// for (auto&t : straight_triples) triples.push_back(t);
  ACTS_INFO("Found " << triples.size() << " candidate triples");
  ACTS_INFO("Finding duplicates");
  std::vector<FW::triple> duplicate_triples = findTriplesDuplicates(ctx, mod, match);//findTriplesDuplicates(pairs, mod);
	//scoreTriples(triples2);

  // const std::vector<triple>* triples = nullptr;

  // if (ctx.eventStore.get("triples", triples) == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Not able to get triples from eventStore in execute");
  // }

	for (auto&t : duplicate_triples) triples.push_back(t);

  // if (ctx.eventStore.add("final_triples", std::move(*triples)) == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Not able to store final triples in eventStore in execute");
  // }

  ACTS_INFO("Found duplicate triples and stored in eventStore");
      //scoreOutsideTriples(triples);
      //scoreOriginTriples(triples);
      //scoreTriples(triples);
  std::vector<FW::triple> final_triples = pruneTriples(ctx, triples);
  ACTS_INFO("Pruned triples");

  std::vector<std::vector<int> > paths = findPaths(ctx, triples, &mod, match);
  ACTS_INFO("Found paths");

  paths = prunePaths(ctx, paths, match);
  ACTS_INFO("Pruned paths")

  paths = addDuplicates(ctx, paths, mod, match);
  ACTS_INFO("Added duplicates");

  paths = findAssignment(ctx, paths, &mod, match);
  ACTS_INFO("Assigned tracks");

  paths = extendPaths(ctx, paths, &mod, match);
  ACTS_INFO("Extended paths");

  ACTS_INFO("Assigning hits to tracks");
  std::map<int, int> assignment_part;
  for (int i = 1; i < paths.size(); i++)
    for (int j : paths[i])
      if (j > 0)
        assignment_part[j] = i;
  std::array<int, 200000> assignment;
  for (std::pair<int, int> p : assignment_part)
    assignment[p.first] = p.second;

  std::map<int, int> map_assignment;
  for (int i = 1; i < hits->size(); i++)
    if (assignment[i]) map_assignment[i] = assignment[i];

  // Write out solution to /tmp/submission0.csv
  writeSubmission(ctx, map_assignment);
  ACTS_INFO("CSV with tracks written");
}

void
FW::Topquarks::readTubes(FW::AlgorithmContext ctx) const
{
  // std::vector<int>*tube = new std::vector<int>[48](); // List of hits in each layer
  std::map<long long,int> added[48];
  const std::vector<int>* metai = nullptr;
  const std::vector<FW::point>* hits = nullptr;
  const std::array<Layer, 48>* layer = nullptr;
  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting metai from eventstore in readTubes");
  }
  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting hits from eventstore in readTubes");
  }
  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting layer from eventstore in readTubes");
  }

  ACTS_VERBOSE("Got metaz, layer and hits for readtubes");

  std::array<std::vector<int>, 48>* tube = new std::array<std::vector<int>, 48>();
  const int Tube = 0, Disc = 1;
  for (int i = 1; i < hits->size(); i++) {
    int j = metai->at(i);
    tube->at(j).push_back(i);
  }
  for (int i = 0; i < 48; i++) {
    ACTS_VERBOSE("Sorting tube vector " << i << " of type " << layer->at(i).type << " having " << tube->at(i).size() << " hits");
    if (layer->at(i).type == Tube)
      sort(tube->at(i).begin(), tube->at(i).end(), [this, ctx](int a, int b) {return z_cmp(a, b, ctx); });
    else
      sort(tube->at(i).begin(), tube->at(i).end(), [this, ctx](int a, int b) {return r_cmp(a, b, ctx); });
  }
  // ACTS_VERBOSE("Sorted tubes. Storing in eventStore");
  FILE*fp17 = fopen("tube", "w");
  // fprintf(fp14, "%d\n", int(pairs->size()));
  for (int i = 0; i < 48; i++){
    ACTS_INFO("Writing " << i << " out of " << 48);
    ACTS_INFO("Vector size " << tube->at(i).size())
    for (auto b = tube->at(i).begin(); b!=tube->at(i).end(); b++){
      fprintf(fp17, "%d\n", tube->at(i)[*b]);
      ACTS_INFO("Writing");
    }
  }
  fclose(fp17);
  if (ctx.eventStore.add("tube", std::move(*tube)) == FW::ProcessCode::ABORT){
    ACTS_INFO("Error in storing tube in eventStore in readTubes");
  }
}

double
FW::Topquarks::wdist(point&a, point&d, double w) const
{
  double pp = a.x*a.x+a.y*a.y+a.z*a.z*w;
  double pd = a.x*d.x+a.y*d.y+a.z*d.z*w;
  double dd = d.x*d.x+d.y*d.y+d.z*d.z*w;
  return sqrt(pp-pd*pd/dd);
}

double
FW::Topquarks::wdistr(double r1, double dr, double az, double dz, double w) const
{
  double pp = r1*r1+az*az*w;
  double pd = r1*dr+az*dz*w;
  double dd = dr*dr+dz*dz*w;
  return sqrt(pp-pd*pd/dd);
}

double
FW::Topquarks::zdist(point&a, point&b)
{
  FW::point origin(0,0,0);
  FW::point p;
  double r;
  circle(origin, a, b, p, r);
  double ang_ab = 2*asin(dist(a.x-b.x, a.y-b.y)*.5/r);
  double ang_a = 2*asin(dist(a.x, a.y)*.5/r);
  return fabs(a.z-(b.z-a.z)*ang_a/ang_ab);
}

double
FW::Topquarks::zdist2(point&a, point&b) const
{
  static point origin(0,0,0);
  point p;
  double r;
  circle(origin, a, b, p, r);
  double ang_ab = 2*asin(dist(a.x-b.x, a.y-b.y)*.5/r);
  double ang_a = 2*asin(dist(a.x, a.y)*.5/r);

  return fabs(b.z-a.z-a.z*ang_ab/ang_a);
}

double
FW::Topquarks::dist(double x, double y) const
{
  return sqrt(x*x + y*y);
}

double
FW::Topquarks::dist2(double x, double y) const
{
  return x*x + y*y;
}

void
FW::Topquarks::circle(const FW::point&a, const FW::point&b, const FW::point&c, FW::point&p, double&r) const
{
  double ax = a.x-c.x, ay = a.y-c.y, bx = b.x-c.x, by = b.y-c.y;
  double aa = ax*ax+ay*ay, bb = bx*bx+by*by;
  double idet = .5/(ax*by-ay*bx);
  p.x = (aa*by-bb*ay)*idet;
  p.y = (ax*bb-bx*aa)*idet;
  p.z = 0;
  r = dist(p.x, p.y);
  p.x += c.x;
  p.y += c.y;
}

bool
FW::Topquarks::r_cmp(int&a, int&b, FW::AlgorithmContext ctx) const
{
  const std::vector<point>* hits = nullptr;
  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get hits from eventStore in z_cmp")
  }
  return hits->at(a).x*hits->at(a).x+hits->at(a).y*hits->at(a).y < hits->at(b).x*hits->at(b).x+hits->at(b).y*hits->at(b).y;
}

bool
FW::Topquarks::z_cmp(int &a, int&b, FW::AlgorithmContext ctx) const
{
  const std::vector<point>* hits = nullptr;
  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get hits from eventStore in z_cmp")
  }
  return hits->at(a).z < hits->at(b).z;
}

double
FW::Topquarks::findDensity(FW::AlgorithmContext ctx, point&dp, point&xp, double target, int li) const
{
  double Ad = 0, A = 0, B = 1, Bd;
  while (1) {
    Bd = getDensity3(ctx, dp, xp, B, li);
    // ACTS_VERBOSE("Came back from getDensity3 in findDensity0 " << Bd << " Target value is " << target);
    if (B > 1e20) {
      ACTS_WARNING("No density?");
      ACTS_VERBOSE("dp.x is " << dp.x << "dp.y is " << dp.y << "dp.z is " << dp.z << "\n xp.x is " << xp.x << "xp.y is " << xp.y << "xp.z is " << xp.z << "\n Value of li is " << li);
      exit(0);
      // return 0;
    }
    if (Bd > target) break;
    B *= 10;
    if (target/Bd < 1e8) B = std::max(B, target/Bd);
  }
  double mid = B/2;
  int cc = 0;
  while (1) {
    double density = getDensity3(ctx, dp, xp, mid, li);
    // ACTS_VERBOSE("Came back from getDensity3 in findDensity1 " << density << " Target value is " << target);
    if (density > target) B = mid, Bd = density;
    else A = mid, Ad = density;

    //cout << A << ' ' << mid << ' ' << B << ' ' << density << endl;
    if ((B-A) < A*1e-3 || density > target*0.9 && density < target*1.1 || cc >= 100) break;
    mid = std::max(A*0.9+B*0.1, std::min(B*0.9+A*0.1, (target-Ad)*(B-A)/(Bd-Ad)+A));
  }

  return mid;
}

double
FW::Topquarks::getDensity3(FW::AlgorithmContext ctx, point&dp, point&xp, double tt, int li) const
{
  const std::array<Layer, 48>* layer = nullptr;
  const std::array<std::vector<double>, 48>* sorted_hits = nullptr;
  const std::array<std::array<std::array<double, 3>, 20000>, 48>* poly = nullptr;

  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting layer from eventStore in getDensity3");
  }
  if (ctx.eventStore.get("sorted_hits", sorted_hits) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting sorted_hits from eventStore in getDensity3");
  }
  if (ctx.eventStore.get("poly", poly) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting poly from eventStore in getDensity3");
  }
  // ACTS_VERBOSE("got layer, poly and sorted_hits from eventStore in getDensity3 for itr", li);
  const Layer&l = layer->at(li);
  int Tube = 0, Disc = 1;
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
  // ACTS_VERBOSE(b << "\n" << tt << "\n" << dx << "\n" << dy << "\n");
  double a = sqrt(1-dx*dx-dy*dy)/((1-dy*dy)*M_PI*dp.x);
  double rx = sqrt(b);
  // ACTS_VERBOSE("Going to getIndex")
  int ai = getIndex(ctx, li, x0-rx);
  int bi = getIndex(ctx, li, x0+rx);
  if (bi-ai > 10) {//Approximate integration by 2. order polynomial approximation to half disc
    //cout << ai << ' ' << bi << endl;
    const double A = 21*M_PI/64., B = -15*M_PI/64.;
    double ib = 1./b;
    double c0 = A+B*x0*x0*ib, c1 = -2*B*ib*x0, c2 = B*ib;
    // ACTS_VERBOSE(A << " " << B << " " << x0 << " " <<  ib);
    double ret =
      ((poly->at(li)[bi][0]-poly->at(li)[ai][0])*c0+
       (poly->at(li)[bi][1]-poly->at(li)[ai][1])*c1+
       (poly->at(li)[bi][2]-poly->at(li)[ai][2])*c2)*a*rx;
    // ACTS_VERBOSE("Check\n");
    // ACTS_VERBOSE(poly->at(li)[bi][0] << "\n" << poly->at(li)[ai][0] << "\n" << c0 << "\n\n");
    // ACTS_VERBOSE(a << "\n" << rx);
    // ACTS_VERBOSE("Returning " << std::max(ret,0.) << " from getDensity3");
    return std::max(ret,0.);
  }
  else { //Exact integration, uses half disc
    double density = 0;
    for(int i = ai; i < bi; i++) {
      // ACTS_VERBOSE("Running else loop in getdensity3. Iterations " << i << " out of " << bi);
      double x = sorted_hits->at(li)[i]-x0;
      double h = a*sqrt(b-x*x);/// *it;
      //cout << h << endl;
      density += h;
    }
    // ACTS_VERBOSE("Returning " << density << " from getDensity3");
    return density;
  }
}


//Angle between line through hits ai-bi and cell's data direction of at hit id ai
double
FW::Topquarks::dir_miss(FW::AlgorithmContext ctx, int ai, int bi) const
{
  const std::array<std::array<point, 2>, 200000>* hit_dir = nullptr;
  const std::vector<point>* hits = nullptr;

  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to get hits from eventStore in dir_miss");
  }
  if (ctx.eventStore.get("hit_dir", hit_dir) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to get hit_dir from eventStore in hit_dir");
  }
  FW::point d = hits->at(ai)-hits->at(bi);
  return acos(std::max(fabs(hit_dir->at(ai).at(0)*d),
		  fabs(hit_dir->at(ai).at(1)*d))/dist0(d));
}

FW::point
FW::Topquarks::normalize(point &a) const
{
  FW::point ret = a*(1./sqrt(a.x*a.x + a.y*a.y + a.z*a.z));
  if (ret.z < 0) ret = ret*-1;
  return ret;
}

void
FW::Topquarks::readEventData(FW::AlgorithmContext ctx) const
{
  // std::vector<int>*tube = new std::vector<int>[48](); // List of hits in each layer
  // const Layer* layer = nullptr;
  // const Detector* detectors = nullptr;
  //
  // if (ctx.eventStore.get("detectors", detectors) == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Error in loading detectors from eventStore");
  //   return FW::ProcessCode::ABORT;
  // }

  // if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Error in loading layer from eventStore");
  //   return FW::ProcessCode::ABORT;
  // }

  // const int* itopo = nullptr; const int* topo = nullptr;
  // if (ctx.eventStore.get("topo", topo) == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Error in loading itopo from eventStore");
  //   return FW::ProcessCode::ABORT;
  // }
  // if (ctx.eventStore.get("itopo", itopo) == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Error in loading topo from eventStore");
  //   return FW::ProcessCode::ABORT;
  // }

  // initialize character array and file name
  char file[1000];
  sprintf(file, "/Users/sharadchitlangia/Desktop/trackml/train_100_events/event000001000-hits.csv");

  FILE*fp = fopen(file, "r");
  if (!fp) {
    ACTS_INFO("couldn't open hits\n");
    exit(1);
  }
  char tmpstr[1000];
  int tmp = fscanf(fp, "%s", tmpstr);

  std::vector<FW::point>* hits = new std::vector<FW::point>();
  std::vector<FW::point>* polar = new std::vector<FW::point>();
  std::vector<FW::point>* meta = new std::vector<FW::point>();
  std::vector<int>* metai = new std::vector<int>();
  std::vector<int>* metaz = new std::vector<int>();
  hits->push_back(point(0,0,0));
  metai->push_back(0);
  meta->push_back(point(0,0,0));
  polar->push_back(point(0,0,0));

  int layers_temp[9] = {7,4,7,6,4,6,6,2,6};
  int metai_list[9][7];
  int c = 0;
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < layers_temp[i]; j++)
      metai_list[i][j] = c++;
  }
  const std::array<int, 48>*itopo;
  std::array<Layer, 48>*layer = new std::array<Layer, 48>();
  const std::array<Layer, 48>*layer_temp = nullptr;

  if (ctx.eventStore.get("layer_temp", layer_temp) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting layer_temp from eventStore in readEventData");
  }

  for (int i = 0; i < 48; i++) {
    layer->at(i).minr = layer->at(i).minz = 1e9;
    layer->at(i).maxr = layer->at(i).maxz =-1e9;
  }

  if (ctx.eventStore.get("itopo", itopo) == FW::ProcessCode::ABORT){
    ACTS_INFO("Error in loading itopo from eventStore in readEventData");
  }
  ACTS_INFO("Reading data")
  for (int i=0; i<48; i++){
    layer->at(i).minr = layer_temp->at(i).minr;
    layer->at(i).maxr = layer_temp->at(i).maxr;
    layer->at(i).avgr = layer_temp->at(i).avgr;
    layer->at(i).minz = layer_temp->at(i).minz;
    layer->at(i).maxz = layer_temp->at(i).maxz;
    layer->at(i).avgz = layer_temp->at(i).avgz;
    layer->at(i).var0 = layer_temp->at(i).var0;
    layer->at(i).var1 = layer_temp->at(i).var1;
    layer->at(i).count = layer_temp->at(i).count;
    layer->at(i).type = layer_temp->at(i).type;
  }
  while (1) {
    // ACTS_INFO("running...")
    // char out_str[50];
    long long hit_id;
    double tx, ty, tz;
    int volume_id, layer_id, module_id;
    if (fscanf(fp, "%lld,%lf,%lf,%lf,%d,%d,%d", &hit_id, &tx, &ty, &tz, &volume_id, &layer_id, &module_id) == -1) break;
    meta->push_back(point(volume_id, layer_id, module_id));
    polar->push_back(point(sqrt(tx*tx+ty*ty), atan2(ty,tx), tz));
    // std::cout << tx << ty << tx << "\n";
    hits->push_back(point(tx, ty, tz));
    // sprintf(out_str, "%lld, %lf, %lf, %lf", hit_id, tx, ty, tz);
    // ACTS_INFO(out_str);
    if (volume_id <= 9)
      volume_id -= 7;
    else if (volume_id <= 14)
      volume_id -= 9;
    else
      volume_id -= 10;
    int mi = itopo->at(metai_list[volume_id][layer_id/2-1]);
    metai->push_back(mi);

    double r = sqrt(tx*tx+ty*ty);
    Layer&l = layer->at(metai_list[volume_id][layer_id/2-1]);
    l.minr = std::min(l.minr, r);
    l.avgr += r;
    l.maxr = std::max(l.maxr, r);
    l.minz = std::min(l.minz, tz);
    l.avgz += tz;
    l.maxz = std::max(l.maxz, tz);
    l.count++;

    // if (hit_id>20000){
    //   break;
    // }
  }
  ACTS_INFO("Read hits data");
  // std::cout << hits->size();

  // if (ctx.eventStore.add("tube", std::move(*tube)) == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Error in storing tube into eventStore");
  //   // return FW::ProcessCode::ABORT;
  // }
  fclose(fp);

  // for (int i=0; i<48; i++){
  //   layer->at(i).minr = layer_temp->at(i).minr;
  //   layer->at(i).maxr = layer_temp->at(i).maxr;
  //   layer->at(i).avgr = layer_temp->at(i).avgr;
  //   layer->at(i).minz = layer_temp->at(i).minz;
  //   layer->at(i).maxz = layer_temp->at(i).maxz;
  //   layer->at(i).avgz = layer_temp->at(i).avgz;
  //   layer->at(i).var0 = layer_temp->at(i).var0;
  //   layer->at(i).var1 = layer_temp->at(i).var1;
  //   layer->at(i).count = layer_temp->at(i).count;
  //   layer->at(i).type = layer_temp->at(i).type;
  // }

  // if (ctx.eventStore.add("layer", std::move(*layer)) == FW::ProcessCode::ABORT){
  //   ACTS_WARNING("Error in adding layer to eventStore");
  //   // return FW::ProcessCode::ABORT;
  // }

  ACTS_INFO("Made final layer structure");

  // initialize structure for cells data
  // pair<pair<ch0, ch1>, value>

  //initialize filename and file pointers
  char file2[1000];
  sprintf(file2, "/Users/sharadchitlangia/Desktop/trackml/train_100_events/event000001000-cells.csv");

  FILE*fp2 = fopen(file2, "r");
  if (!fp2) {
    printf("couldn't open cells\n");
    exit(1);
  }
  std::array<std::vector<std::pair<std::pair<int, int>, double> > ,200000>* hit_cells = new std::array<std::vector<std::pair<std::pair<int, int>, double> >, 200000>();
  ACTS_INFO("Reading cells data");
  char tmpstr2[1000];
  int tmp2 = fscanf(fp2, "%s", tmpstr);
  int hit_id, ch0, ch1;
  double value;
  while (fscanf(fp2, "%d,%d,%d,%lf", &hit_id, &ch0, &ch1, &value) == 4) {
    hit_cells->at(hit_id).push_back(std::make_pair(std::make_pair(ch0, ch1), value));
    // if (hit_id>20000){
    //   break;
    // }
  }
  fclose(fp2);
  ACTS_INFO("Read cells data");
  // Read the input space points
  // const DetectorData<geo_id_value, Acts::Vector3D>* spacePoints;
  // if (ctx.eventStore.get(m_cfg.spacePointCollection, spacePoints)
  //     == FW::ProcessCode::ABORT) {
  //   ACTS_WARNING("missing space point input");
  //   return FW::ProcessCode::ABORT;
  // }
  // ACTS_INFO("Loaded input points for event " << ctx.eventNumber);

  std::map<int, Detector> detectors;

  char file3[1000];
  sprintf(file3, "/Users/sharadchitlangia/Desktop/trackml/detectors.csv");
  FILE*fp3 = fopen(file3, "r");
  if (!fp3) {
    printf("couldn't open detectors\n");
    exit(1);
  }

  char tmpstr3[1000];
  int tmp3 = fscanf(fp3, "%s", tmpstr3);
  //cout << tmpstr << endl;

  Detector d;
  while (fscanf(fp3, "%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &d.volume_id, &d.layer_id, &d.module_id, &d.c.x, &d.c.y, &d.c.z,
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
  fclose(fp3);
  ACTS_INFO("Read Detectors");

  std::array<std::array<FW::point, 2>, 200000>* hit_dir = new std::array<std::array<FW::point, 2>, 200000>();
  for (int hit_id = 1; hit_id < hits->size(); hit_id++) {
    if (hit_id%100 == 0){
      ACTS_VERBOSE(hit_id << " done out of " << hits->size());
    }
      //ACTS_VERBOSE("Running " << hit_id << " out of " << hits->size());
    point m = meta->at(hit_id);
    Detector&d = detectors[int(m.x)*10000000+int(m.y)*10000+int(m.z)];

    //if (!hit_cells[hit_id].size()) cout << "Hit with zero cells" << endl;
    //if (metai[hit_id] < 18) continue;

    //Use linear regression for direction
    double mx = 0, my = 0, mw = 0;
    auto&cells = hit_cells->at(hit_id);
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
    }
    //Find eigenvector with minimum eigenvalue
    double a = mxx-myy, b = 2*mxy;
    double x = a+b+sqrt(a*a+b*b);
    double y =-a+b+sqrt(a*a+b*b);
    if (0) {
      double lambda = (mxx+myy+sqrt(a*a+b*b))/2;
    }

    //Analytical formula for z
    double z = 2*d.d*(fabs(x)/d.cell_w+fabs(y)/d.cell_h+1e-8);
    x *= (cells.size()*1.-1.3);//1.3 != 1 was adjusted empirically
    y *= (cells.size()*1.-1.3);
    FW::point d1(x,y,z), d2(x,y,-z);
    d1 = d.rx*d1.x+d.ry*d1.y+d.rz*d1.z;
    d2 = d.rx*d2.x+d.ry*d2.y+d.rz*d2.z;
    hit_dir->at(hit_id).at(0) = normalize(d1);
    hit_dir->at(hit_id).at(1) = normalize(d2);

    for (int k = 0; k < 2; k++)
      if (hit_dir->at(hit_id).at(k)*hits->at(hit_id) < 0)
  	hit_dir->at(hit_id).at(k) = hit_dir->at(hit_id).at(k)*-1;
  }

  // for (int hit_id = 1; hit_id < hits->size(); hit_id++) {
  //   metai_weight[truth_part[hit_id]][metai[hit_id]] += truth_weight[hit_id];
  // }
  // const std::array<Layer, 48>* layer = nullptr;
  // if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
  //   ACTS_INFO("Error in getting layer from event store in readEventData");
  // }
  // std::cout<< check->at(0)<<"\n";
  std::array<std::array<double, 4>, 48> *disc_z = new std::array<std::array<double, 4>, 48>();
  std::map<double, double> mir[48], mar[48];
  for (int i = 1; i < hits->size(); i++) {
    int mi = metai->at(i);
    if (layer->at(mi).type != Disc) continue;
    double&mir_ = mir[mi][polar->at(i).z];
    double&mar_ = mar[mi][polar->at(i).z];
    if (!mir_) mir_ = 1e9;
    mir_ = std::min(mir_, polar->at(i).x);
    mar_ = std::max(mar_, polar->at(i).x);
  }
  // std::cout<< FW::Topquarks::layer->at(0).type <<"\n";
  double z_minr[48][4], z_maxr[48][4];
  std::map<double, int> zi[48];
  for (int mi = 0; mi < 48; mi++) {
    if (layer->at(mi).type != Disc) continue;
    int k = 0;
    for (auto p : mir[mi]) {
      double mir_ = mir[mi][p.first]-1e-5;
      double mar_ = mar[mi][p.first]+1e-5;
      z_minr[mi][k] = mir_;
      z_maxr[mi][k] = mar_;
      const double temp = p.first;
      disc_z->at(mi)[k] = temp;
      zi[mi][p.first] = k++;
    }
  }

  metaz->resize(hits->size());
  metaz->at(0) = 0;
  // ACTS_VERBOSE("Initiating metaz loop");
  for (int i = 1; i < hits->size(); i++) {
    int mi = metai->at(i);
    metaz->at(i) = meta->at(i).z;
    if (layer->at(mi).type == Disc){
      metaz->at(i) = zi[mi][hits->at(i).z];
    }
  }
  ACTS_INFO("Adding read objects to eventStore");
  if (ctx.eventStore.add("detectors", std::move(detectors)) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in adding detectors to eventStore");
    // return FW::ProcessCode::ABORT;
  }
  if (ctx.eventStore.add("hits", std::move(*hits)) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in storing hits");
    // return FW::ProcessCode::ABORT;
  }
  if (ctx.eventStore.add("hit_cells", std::move(*hit_cells)) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in addong cells to eventStore in readEventData");
    // return FW::ProcessCode::ABORT;
  }
  if (ctx.eventStore.add("metai", std::move(*metai)) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in adding metai to eventStore in readEventData");
    // return FW::ProcessCode::ABORT;
  }
  if (ctx.eventStore.add("meta", std::move(*meta)) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in adding meta to eventStore in readEventData");
  }
  if (ctx.eventStore.add("polar", std::move(*polar)) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in adding polar to eventStore in readEventData");
  }
  if (ctx.eventStore.add("metaz", std::move(*metaz)) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in adding metaz to eventStore in readEventData");
  }
  if (ctx.eventStore.add("disc_z", std::move(*disc_z)) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in adding dics_z to eventStore in readEventData");
  }
  if (ctx.eventStore.add("hit_dir", std::move(*hit_dir)) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in adding hit_dir to eventStore in readEventData");
  }
  if (ctx.eventStore.add("layer", std::move(*layer)) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in adding layer to eventStore in readEventData");
  }
}

void
FW::Topquarks::readDetectorInfo(FW::AlgorithmContext ctx) const
{
  // ACTS_INFO("No errors");
  const int Tube = 0, Disc = 1;
  // Layer layer[48];
  double z_minr[48][4], z_maxr[48][4];
  // int topo[48], itopo[48];
  // ACTS_INFO("No errors1");
  int handpicked[48] = {};
  int c = 0;
  for (int i = 0; i < 4; i++) handpicked[c++] = 7+i;
  for (int i = 0; i < 7; i++) handpicked[c++] = 7-1-i;
  for (int i = 0; i < 7; i++) handpicked[c++] = 11+i;
  // ACTS_INFO("No errors2");
  for (int i = 0; i < 4; i++) handpicked[c++] = 24+i;
  for (int i = 0; i < 2; i++) handpicked[c++] = 40+i;
  for (int i = 0; i < 6; i++) {
    handpicked[c++] = 24-1-i;
    handpicked[c++] = 40-1-i;
  }
  // ACTS_INFO("No errors3");
  for (int i = 0; i < 6; i++) {
    handpicked[c++] = 28+i;
    handpicked[c++] = 42+i;
  }
  // topo->clear();
  // topo->reserve(48);
  // itopo->clear();
  // itopo->reserve(48);

  std::array<int, 48>*topo = new std::array<int, 48>();
  std::array<int, 48>*itopo = new std::array<int, 48>();

  for (int i = 0; i < 48; i++) {
    topo->at(i) = handpicked[i];
    itopo->at(topo->at(i)) = i;
  }

  // int* p_itopo = itopo;
  // int* p_topo = topo;
  //
  // // topo and itopo exist. Put them in AlgorithmContext
  //

  if (ctx.eventStore.add("itopo", std::move(*itopo)) == FW::ProcessCode::ABORT){
   ACTS_VERBOSE("Error in storing itopo");
   // return FW::ProcessCode::ABORT;
  }

  // std::string topo_name = "topo";
  if (ctx.eventStore.add("topo", std::move(*topo)) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in storing topo");
    // return FW::ProcessCode::ABORT;
  }


  // ACTS_INFO("Stored topo and itopo");
  double avgz1[2][7] = {{-1500,-1300,-1100,-960,-820,-700,-600},
			{ 600, 700, 820, 960, 1100, 1300, 1500}};

  std::array<Layer, 48>* layer = new std::array<Layer, 48>();

  // for (int i = 0; i < 48; i++) {
  //   layer[i].minr = layer[i].minz = 1e9;
  //   layer[i].maxr = layer[i].maxz =-1e9;
  // }

  for (int k = 0; k < 2; k++)
    for (int i = 0; i < 7; i++) {
      layer->at(k*11+i).minr = 30;
      layer->at(k*11+i).maxr = 176.5;
      layer->at(k*11+i).avgz = avgz1[k][i];
      layer->at(k*11+i).type = Disc;
    }
  double avgz2[2][6] = {{-2950,-2550,-2150,-1800,-1500,-1220},
			{ 1220, 1500, 1800, 2150, 2550, 2950}};
  for (int k = 0; k < 2; k++)
    for (int i = 0; i < 6; i++) {
      layer->at(k*10+i+18).minr = 240;
      layer->at(k*10+i+18).maxr = 701;
      layer->at(k*10+i+18).avgz = avgz2[k][i];
      layer->at(k*10+i+18).type = Disc;

      layer->at(k*8+i+34).minr = 755;
      layer->at(k*8+i+34).maxr = 1018;
      layer->at(k*8+i+34).avgz = avgz2[k][i];
      layer->at(k*8+i+34).type = Disc;
    }

  double avgr1[4] = {32.3, 72.1, 116.1, 172.1};
  double avgr2[4] = {260.3, 360.2, 500.2, 660.2};
  double avgr3[2] = {820.2, 1020.2};

  for (int i = 0; i < 4; i++) {
    layer->at(i+7).minz =-491;
    layer->at(i+7).maxz = 491;
    layer->at(i+7).avgr = avgr1[i];
    layer->at(i+7).type = Tube;
  }
  for (int i = 0; i < 4; i++) {
    layer->at(i+24).minz =-1084;
    layer->at(i+24).maxz = 1084;
    layer->at(i+24).avgr = avgr2[i];
    layer->at(i+24).type = Tube;
  }
  for (int i = 0; i < 2; i++) {
    layer->at(i+40).minz =-1084;
    layer->at(i+40).maxz = 1084;
    layer->at(i+40).avgr = avgr3[i];
    layer->at(i+40).type = Tube;
  }
  Layer layer2[48];
  for (int i = 0; i < 48; i++) layer2[i] = layer->at(i);
  for (int i = 0; i < 48; i++) layer->at(i) = layer2[topo->at(i)];

  layer->at(0).var0 = 1e-3;
  layer->at(1).var0 = 5e-4;
  for (int i = 2; i < 18; i++) layer->at(i).var0 = 3e-4;
  for (int i = 18; i < 22; i++) layer->at(i).var0 = 5e-2;
  for (int i = 22; i < 48; i++) layer->at(i).var0 = i%2 || i == 22 ? 9 : 0.1;

  for (int i = 0; i < 4; i++) layer->at(i).var1 = 0.5;
  for (int i = 4; i < 18; i++) layer->at(i).var1 = 5;
  for (int i = 18; i < 24; i++) layer->at(i).var1 = 7;
  for (int i = 24; i < 48; i++) layer->at(i).var1 = i%2 ? 19 : 11;
  // delete layer2;
  // Layer* p_layer = layer;
  //
  if (ctx.eventStore.add("layer_temp", std::move(*layer)) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in adding layer_temp to eventStore");
    // return FW::ProcessCode::ABORT;
  }

  ACTS_INFO("Stored layer_temp in eventStore");
}

// void
// FW::Topquarks::initLayers(FW::AlgorithmContext ctx) const
// {
//   const std::array<int, 48>*topo = nullptr;
//   const std::vector<FW::point>* hits = nullptr;
//
//   if (ctx.eventStore.get("topo", topo) == FW::ProcessCode::ABORT){
//     ACTS_WARNING("Error in getting topo from eventStore in initLayers");
//   }
//   if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
//     ACTS_WARNING("Error in getting hits from eventStore in initLayers");
//   }
//   double avgz1[2][7] = {{-1500,-1300,-1100,-960,-820,-700,-600},
// 			{ 600, 700, 820, 960, 1100, 1300, 1500}};
//
//   std::array<Layer, 48>* layer = new std::array<Layer, 48>();
//
//
//
//   for (int i = 0; i < 48; i++) {
//     layer->at(i).minr = layer->at(i).minz = 1e9;
//     layer->at(i).maxr = layer->at(i).maxz =-1e9;
//   }
//
//   for (int k = 0; k < 2; k++)
//     for (int i = 0; i < 7; i++) {
//       layer->at(k*11+i).minr = 30;
//       layer->at(k*11+i).maxr = 176.5;
//       layer->at(k*11+i).avgz = avgz1[k][i];
//       layer->at(k*11+i).type = Disc;
//     }
//   double avgz2[2][6] = {{-2950,-2550,-2150,-1800,-1500,-1220},
// 			{ 1220, 1500, 1800, 2150, 2550, 2950}};
//   for (int k = 0; k < 2; k++)
//     for (int i = 0; i < 6; i++) {
//       layer->at(k*10+i+18).minr = 240;
//       layer->at(k*10+i+18).maxr = 701;
//       layer->at(k*10+i+18).avgz = avgz2[k][i];
//       layer->at(k*10+i+18).type = Disc;
//
//       layer->at(k*8+i+34).minr = 755;
//       layer->at(k*8+i+34).maxr = 1018;
//       layer->at(k*8+i+34).avgz = avgz2[k][i];
//       layer->at(k*8+i+34).type = Disc;
//     }
//
//   double avgr1[4] = {32.3, 72.1, 116.1, 172.1};
//   double avgr2[4] = {260.3, 360.2, 500.2, 660.2};
//   double avgr3[2] = {820.2, 1020.2};
//
//   for (int i = 0; i < 4; i++) {
//     layer->at(i+7).minz =-491;
//     layer->at(i+7).maxz = 491;
//     layer->at(i+7).avgr = avgr1[i];
//     layer->at(i+7).type = Tube;
//   }
//   for (int i = 0; i < 4; i++) {
//     layer->at(i+24).minz =-1084;
//     layer->at(i+24).maxz = 1084;
//     layer->at(i+24).avgr = avgr2[i];
//     layer->at(i+24).type = Tube;
//   }
//   for (int i = 0; i < 2; i++) {
//     layer->at(i+40).minz =-1084;
//     layer->at(i+40).maxz = 1084;
//     layer->at(i+40).avgr = avgr3[i];
//     layer->at(i+40).type = Tube;
//   }
//   Layer layer2[48];
//   for (int i = 0; i < 48; i++) layer2[i] = layer->at(i);
//   for (int i = 0; i < 48; i++) layer->at(i) = layer2[topo->at(i)];
//
//   layer->at(0).var0 = 1e-3;
//   layer->at(1).var0 = 5e-4;
//   for (int i = 2; i < 18; i++) layer->at(i).var0 = 3e-4;
//   for (int i = 18; i < 22; i++) layer->at(i).var0 = 5e-2;
//   for (int i = 22; i < 48; i++) layer->at(i).var0 = i%2 || i == 22 ? 9 : 0.1;
//
//   for (int i = 0; i < 4; i++) layer->at(i).var1 = 0.5;
//   for (int i = 4; i < 18; i++) layer->at(i).var1 = 5;
//   for (int i = 18; i < 24; i++) layer->at(i).var1 = 7;
//   for (int i = 24; i < 48; i++) layer->at(i).var1 = i%2 ? 19 : 11;
//   // delete layer2;
//   // Layer* p_layer = layer;
//   //
//   if (ctx.eventStore.add("layer", std::move(*layer)) == FW::ProcessCode::ABORT){
//     ACTS_VERBOSE("Error in adding 'layer' to eventStore");
//     // return FW::ProcessCode::ABORT;
//   }
//
//   ACTS_INFO("Stored layer in eventStore");
// }


// std::vector<std::pair<int, int> >
void
FW::Topquarks::findPairs(FW::AlgorithmContext ctx, int modeli = 3) const
{
  // ctx.eventStore.get("hit_cells", hit_cells);
  //ctx.eventStore.get("layer", layer);
  //ctx.eventStore.get("metai", metai);
  const int features = 6;
  double feature[features];
  const std::vector<FW::point>* hits = nullptr;
  const std::vector<FW::point>* polar = nullptr;
  const std::vector<FW::point>* meta = nullptr;
  const std::vector<int>* metai = nullptr;
  const std::vector<int>* metaz = nullptr;
  const std::array<std::vector<int>, 48>* tube = nullptr;
  const std::array<Layer, 48>* layer = nullptr;
  // const std::vector<int> hits;

  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load hits in findPairs");
  }
  if (ctx.eventStore.get("polar", polar) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load polar in findPairs");
  }
  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load metai in findPairs");
  }
  if (ctx.eventStore.get("metaz", metaz) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load metaz in findPairs");
  }
  if (ctx.eventStore.get("meta", meta) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load meta in findPairs");
  }
  if (ctx.eventStore.get("tube", tube) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load tube in findPairs");
  }
  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    ACTS_INFO("Not able to load layer in findPairs");
  }



  // std::array<std::vector<int>, 48>* tube = new std::array<std::vector<int>, 48>();
  // for (int i = 1; i < hits->size(); i++) {
  //   // ACTS_VERBOSE("Assigning the tubes");
  //   int j = metai->at(i);
  //   tube->at(j).push_back(i);
  // }
  // ACTS_VERBOSE("Sorting tube vectors")
  // for (int i = 0; i < 48; i++) {
  //   ACTS_VERBOSE("Sorted " << i << " out of " << 48);
  //   if (layer->at(i).type == Tube)
  //     sort(tube->at(i).begin(), tube->at(i).end(), [this, ctx](int a, int b) {return z_cmp(a, b, ctx); });
  //   else
  //     sort(tube->at(i).begin(), tube->at(i).end(), [this, ctx](int a, int b) {return r_cmp(a, b, ctx); });
  // }

  ACTS_VERBOSE("Got objects from eventStore in find in findPairs");

  const int n = 10;//How many pairs of layers to consider. Roughly proportional to run-time, and setting this to 30 gave practically the same score (less than 0.0002 reduction)
  std::pair<int, int> start_list[100] = {{0, 1}, {11, 12}, {4, 5}, {0, 4}, {0, 11}, {18, 19}, {1, 2}, {5, 6}, {12, 13}, {13, 14}, {6, 7}, {2, 3}, {3, 18}, {19, 20}, {0, 2}, {20, 21}, {1, 4}, {7, 8}, {11, 18}, {1, 11}, {14, 15}, {4, 18}, {2, 18}, {21, 22}, {0, 18}, {1, 18}, {24, 26}, {36, 38}, {15, 16}, {8, 9}, {22, 23}, {9, 10}, {16, 17}, {38, 40}, {5, 18}, {18, 24}, {18, 36}, {12, 18}, {40, 42}, {28, 30}, {26, 28}, {0, 12}, {18, 20}, {6, 18}, {2, 11}, {13, 18}, {2, 4}, {0, 5}, {19, 36}, {19, 24}, {4, 6}, {19, 22}, {20, 22}, {11, 13}, {3, 19}, {7, 18}, {14, 18}, {3, 4}, {22, 25}, {1, 3}, {20, 24}, {15, 18}, {3, 11}, {22, 37}, {30, 32}, {42, 44}, {8, 18}, {9, 18}, {8, 26}, {15, 38}, {20, 36}, {14, 36}, {7, 24}, {1, 5}, {16, 18}, {22, 24}, {18, 22}, {25, 27}, {16, 40}, {10, 30}, {25, 26}, {17, 40}, {36, 39}, {1, 12}, {10, 28}, {7, 26}, {17, 42}, {24, 27}, {21, 24}, {23, 37}, {13, 36}, {15, 36}, {22, 36}, {14, 38}, {8, 28}, {19, 21}, {6, 24}, {9, 28}, {16, 38}, {0, 3}};
  std::vector<std::pair<int, int> > *pairs = new std::vector<std::pair<int, int> >();
  for (int i = 0; i < n; i++) {
    //load(n, "Finding pairs");
    ACTS_VERBOSE("Finding pairs for layer " << i);
    double weight[23];
    {
      char filename[100];
      int modeli = 3; // Config
      sprintf(filename, "/Users/sharadchitlangia/Desktop/trackml/code/acts-framework/Algorithms/Reconstruction/ACTFW/Top-quarks/trained/pair_logreg%d/model%d", modeli, i);
      FILE*fp = fopen(filename, "r");
      if (!fp) {
        ACTS_VERBOSE("Could not open trained model number " << i);
        continue;
      }
      for (int i = 0; i < 23; i++) int tmp = fscanf(fp, "%lf", &weight[i]);
      fclose(fp);
      ACTS_VERBOSE("Loaded model " << i);
    }
    // ACTS_INFO(tube->at(10).size());
    // std::cout << std::endl << start_list[0].first << std::endl;
    // std::cout << std::endl << tube->at(start_list[0].first)[0] << std::endl;
    // int a = 0;
    // int b = 0;
    for (auto a = tube->at(start_list[i].first).begin(); a!=tube->at(start_list[i].first).end(); a++) {
      // ACTS_VERBOSE("Running loop " << *a << " out of " << tube->at(start_list[i].first).size());
      for (auto b = tube->at(start_list[i].second).begin(); b!=tube->at(start_list[i].second).end(); b++) {
        // ACTS_VERBOSE(*a << " " << *b);
        //ACTS_VERBOSE("Running inner loop " << *b << " out of " << tube->at(start_list[i].second).size());
        // if ((*a==16938) && (*b == 13850)){
        //   ACTS_INFO("Got in " << *a);}
        // if ((*a==77208) && (*b == 85109)){
        //   ACTS_INFO("Got in " << *a);}
        // if ((*a==17181) && (*b == 14175)){
        //   ACTS_INFO("Got in " << *a);}
        double dot = hits->at(*a).x * hits->at(*b).x+ hits->at(*a).y * hits->at(*b).y;
        double alen = hits->at(*a).x*hits->at(*a).x + hits->at(*a).y*hits->at(*a).y;
        double blen = hits->at(*b).x*hits->at(*b).x + hits->at(*b).y*hits->at(*b).y;
        if (dot < 0 || dot*dot < alen*blen*(.7*.7)) continue;
        // if (((*a==16938) && (*b == 13850)) || ((*a==77208) && (*b == 85109)) || ((*a==17181) && (*b == 14175))){
        //   ACTS_INFO("Surpassed first if " << *a);
        //   ACTS_INFO(dot << " " << alen << " " << blen);}

        dot += hits->at(*a).z*hits->at(*b).z;
        alen += hits->at(*a).z*hits->at(*a).z;
        blen += hits->at(*b).z*hits->at(*b).z;
        // if (((*a==16938) && (*b == 13850)) || ((*a==77208) && (*b == 85109)) || ((*a==17181) && (*b == 14175))){
        //   ACTS_INFO(dot << " " << alen << " " << blen);}
        if (dot < 0 || dot*dot < alen*blen*(.7*.7)) continue;
        // if (((*a==16938) && (*b == 13850)) || ((*a==77208) && (*b == 85109)) || ((*a==17181) && (*b == 14175))){
        //   ACTS_INFO("Surpassed second if " << *a);}

        //getFeatures3(a, b, feature);

        int x = metai->at(*a), y = metai->at(*b);
        FW::point p_a = hits->at(*a), p_b = hits->at(*b);
        FW::point d = p_a-p_b;

        double dr2 = sqrt(d.x*d.x + d.y*d.y);
        double r1 = sqrt(p_a.x*p_a.x + p_a.y*p_a.y);
        double r2 = sqrt(p_b.x*p_b.x + p_b.y*p_b.y);
        double dr = r2-r1;

        feature[0] = dir_miss(ctx, *a, *b);//Cell's data of ai
        feature[1] = dir_miss(ctx, *b, *a);//Cell's data of bi
        //Different distances from origin (like how far does the line through ai-bi pass from the origin)
        feature[2] = wdist(p_a, d, 0);
        feature[3] = zdist2(p_a, p_b);
        feature[4] = wdistr(r1, dr, p_a.z, d.z, 1);
        feature[5] = wdist(p_a, d, 1);
        double&A = feature[0], &B = feature[1], &C = feature[2], &D = feature[3], &E = feature[4], &F = feature[5];
        double score = weight[0]+
            A*(weight[1]+
               B*weight[7]+
               C*weight[8]+
               D*weight[9]+
               E*weight[10]+
               F*weight[11])+
            B*(weight[2]+
               C*weight[12]+
               D*weight[13]+
               E*weight[14]+
               F*weight[15])+
            C*(weight[3]+
               D*weight[16]+
               E*weight[17]+
               F*weight[18])+
            D*(weight[4]+
               E*weight[19]+
               F*weight[20])+
            E*(weight[5]+
               F*weight[21])+
            F*weight[6];
        if (score > weight[22]){
          pairs->push_back(std::make_pair(*a, *b));
        }
        // if ((*a==16938) && (*b == 13850)){
        //   ACTS_INFO(A << " " << B << " " << C << " " << D << " " << E << " " << F);}
        // if ((*a==77208) && (*b == 85109)){
        //   ACTS_INFO(A << " " << B << " " << C << " " << D << " " << E << " " << F);}
        // if ((*a==17181) && (*b == 14175)){
        //   ACTS_INFO(A << " " << B << " " << C << " " << D << " " << E << " " << F);}
        }
      }
      ACTS_VERBOSE("Found pairs for model " << i);
    }

    // FILE*fp5 = fopen("pairs", "w");
    // fprintf(fp5, "%d\n", int(pairs->size()));
    // for (int i = 0; i < pairs->size(); i++)
    //   fprintf(fp5, "%d %d\n", pairs->at(i).first, pairs->at(i).second);
    // fclose(fp5);

    // FILE*fp6 = fopen("hits", "w");
    // fprintf(fp5, "%d\n", int(pairs->size()));
    // for (int i = 0; i < pairs->size(); i++)
    //   fprintf(fp5, "%d %d\n", pairs->at(i).first, pairs->at(i).second);
    // fclose(fp5);
    //
    const std::array<std::array<FW::point, 2>, 200000>* hit_dir = nullptr;
    if (ctx.eventStore.get("hit_dir", hit_dir) == FW::ProcessCode::ABORT){
      ACTS_VERBOSE("Error in getting hit_dir from eventStore in scoreTripleLogRadius_and_HitDir");
    }
    FILE*fp7 = fopen("hit_dir", "w");
    fprintf(fp7, "%d\n", 200000);
    for (int i = 0; i < 200000; i++)
      fprintf(fp7, "%lf %lf %lf %lf %lf %lf\n", hit_dir->at(i)[0].x, hit_dir->at(i)[0].y, hit_dir->at(i)[0].z, hit_dir->at(i)[1].x, hit_dir->at(i)[1].y, hit_dir->at(i)[1].z);
    fclose(fp7);
    //
    // FILE*fp8 = fopen("polar", "w");
    // fprintf(fp5, "%d\n", int(pairs->size()));
    // for (int i = 0; i < pairs->size(); i++)
    //   fprintf(fp5, "%d %d\n", pairs->at(i).first, pairs->at(i).second);
    // fclose(fp5);
    //
    // FILE*fp9 = fopen("meta", "w");
    // fprintf(fp5, "%d\n", int(pairs->size()));
    // for (int i = 0; i < pairs->size(); i++)
    //   fprintf(fp5, "%d %d\n", pairs->at(i).first, pairs->at(i).second);
    // fclose(fp5);
    //
    // FILE*fp10 = fopen("metai", "w");
    // fprintf(fp5, "%d\n", int(pairs->size()));
    // for (int i = 0; i < pairs->size(); i++)
    //   fprintf(fp5, "%d %d\n", pairs->at(i).first, pairs->at(i).second);
    // fclose(fp5);
    //
    // FILE*fp11 = fopen("metaz", "w");
    // fprintf(fp5, "%d\n", int(pairs->size()));
    // for (int i = 0; i < pairs->size(); i++)
    //   fprintf(fp5, "%d %d\n", pairs->at(i).first, pairs->at(i).second);
    // fclose(fp5);
    //
    // FILE*fp12 = fopen("layer", "w");
    // fprintf(fp5, "%d\n", int(pairs->size()));
    // for (int i = 0; i < pairs->size(); i++)
    //   fprintf(fp5, "%d %d\n", pairs->at(i).first, pairs->at(i).second);
    // fclose(fp5);


    if (ctx.eventStore.add("pairs", std::move(*pairs)) == FW::ProcessCode::ABORT){
      ACTS_VERBOSE("Error in storing hits");
      // return FW::ProcessCode::ABORT;
    }
  // return pairs;
}

std::vector<FW::triple>
FW::Topquarks::findTriples(FW::AlgorithmContext ctx, std::array<PolarModule, 48>& mod, std::array<int, 200000> &match, int method = 0, double target = 0.5) const
{
  std::vector<FW::triple> triples;
  // std::vector<double> v_temp[48];
  // const std::vector<std::pair<int, int> >* pairs = nullptr;
  const std::vector<FW::point>* hits = nullptr;
  const std::vector<int>* metai = nullptr;
  const std::array<std::array<int, 48>, 48>* next_layer = nullptr;
  //
  if (ctx.eventStore.get("next_layer", next_layer) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in loading next_layer from eventStore");
    // return FW::ProcessCode::ABORT;
  }
  //
  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in loading hits from eventStore");
    // return FW::ProcessCode::ABORT;
  }
  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in loading metai from the eventStore");
    // return FW::ProcessCode::ABORT;
  }
  const std::vector<std::pair<int, int> >* pairs = nullptr;

  if (ctx.eventStore.get("pairs", pairs) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get pairs from eventStore in findTriples");
  }
  ACTS_VERBOSE("Starting iteration over pairs to find triples with method " << method);
  int i = 0;
  for (auto p : *pairs) {
    i++;
    if (i%10000 == 9999){
      ACTS_VERBOSE("Iteration " << ++i << " out of " << pairs->size());
    }
    int ai = p.first, bi = p.second;
    if ((ai == 16880) && (bi == 25841)){
      ACTS_INFO("Got into " << ai << " " << bi);
    }
    //if (!samepart(ai, bi)) continue;
    const FW::point&a = hits->at(ai), &b = hits->at(bi);
    // ACTS_VERBOSE("1");
    int added = 0;
    for (int li = metai->at(bi)+1; added < 100*method+1 && li < 48; li++) {
      // ACTS_VERBOSE("Running inner loop iteration " << li << " out of " << 48 << ". Outer Loop iter: " << i << " out of " << pairs->size());
      if (next_layer->at(metai->at(bi))[li] < adj_thres){
        // ACTS_VERBOSE("Continuing");
        continue;
      }
      if ((ai == 16880) && (bi == 25841)) {ACTS_INFO("Surpassed initial condition");}
      if (method == 1 && extendTripleLine(ctx, triples, ai, bi, a, b, li, mod[li], match, 0)) added++;
      if (method == 0 && extendTripleOrigin(ctx, triples, ai, bi, a, b, li, mod[li], match, target, 0)) added++;
    }
  }
  return triples;
}

void
FW::Topquarks::loadAdjacent(FW::AlgorithmContext ctx) const
{
  FILE*fp = fopen("/Users/sharadchitlangia/Desktop/trackml/code/top-quarks/trained/adjacency", "r");
  if (!fp) {
    ACTS_INFO("Could not open adjacency");
    exit(0);
  }
  std::array<std::array<int, 48>, 48> next_layer;
  for(int i=0; i<48; i++){
    for(int j=0; j<48; j++){
      int tmp = fscanf(fp, "%d", &next_layer[i][j]);
    }
  }

  if (ctx.eventStore.add("next_layer", std::move(next_layer)) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in storing next_layer in eventStore");
  }
  ACTS_VERBOSE("Stored next_layer in eventStore");
  fclose(fp);
}

double
FW::Topquarks::extendTripleLine(FW::AlgorithmContext ctx, std::vector<FW::triple>&triples, int ai, int bi, const point&a, const point&b, int li, PolarModule&mod, std::array<int, 200000> &match, int rev = 0) const
{
  // const std:vector<FW::triple>* triples = nullptr;
  //
  // if (ctx.eventStore.get("triples", triples) == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Not able to get triples from eventStore in extendTripleLine");
  // }

  FW::point d, dp, xp, bap;// Where are these declared?
  // ACTS_VERBOSE("Reached extendTripleLine")
  if (prepareTripleScoreDefault(ctx, ai, bi, li, d, dp, xp, bap)) return 0;

  //xp = normalize(xp)*0.999;
  //xp = point(0,0,0.5);
  const double target0 = 0.1, target = 10;//0.5;
  // ACTS_VERBOSE("Not returning 0");
  //double mid0 = -findDensity(dp, xp, target0, li);
  double mid = findDensity(ctx, dp, xp, target, li);
  // ACTS_VERBOSE("Found density" << mid);
  int matches = mod.getNear(ctx, dp, xp, bap, mid, match, li);
  // ACTS_VERBOSE("Came back from get Near")
  int mini;
  double best = target;
  std::vector<std::pair<double, int> > v;
  for (int i = 0; i < matches; i++) {
    int ci = match.at(i);
    //if (ci == ai || ci == bi) continue;
    double s = scoreTriple(ctx, ai, bi, ci);//evaluateScore(ci, dp, xp, bap);//
    v.push_back(std::make_pair(s, ci));
  }
  //Take only best ones
  std::sort(v.begin(), v.end());
  for (int i = 0; i < v.size(); i++) {
    if (i >= target) break;
    FW::triple t(ai, bi, v[i].second);
    if (rev) std::swap(t.x,t.z);
    if (acceptTriple(ctx, t))
    { //Prune most triples
      // ACTS_VERBOSE("Found a triplet");
      triples.push_back(t);
    }
  }

  // if (ctx.eventStore.add("triples", std::move(*triples)) == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Not able to store triples into eventStore in extendTripleLine");
  // }

  return 1;
}

double
FW::Topquarks::extendTripleOrigin(FW::AlgorithmContext ctx, std::vector<FW::triple>&triples, int ai, int bi, const point&a, const point&b, int li, PolarModule&mod, std::array<int, 200000> &match, double target = 0.5, int rev = 0) const
{
  // polar->at(0) = point(0,0,0);


  FW::point d, dp, xp, bap;
  if (prepareQuadrupleScoreDefault(ctx, 0, ai, bi, li, d, dp, xp, bap, 1)) return 0;
  if ((ai == 16880) && (bi == 25841)) {ACTS_INFO("Surpassed second condition for " << ai << " " << bi);}

  // const std:vector<FW::triple>* triples = nullptr;
  //
  // if (ctx.eventStore.get("triples", triples) == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Not able to get triples from eventStore in extendTripleLine");
  // }

  //const double target0 = 0.1, target = 0.5;

  //double mid0 = findDensity(dp, xp, target0, li);
  double mid = findDensity(ctx, dp, xp, target, li);
  if ((ai == 16880) && (bi == 25841)) ACTS_INFO(mid << " for " << ai << " " << bi);

  int matches = mod.getNear(ctx, dp, xp, bap, mid, match, li);
  std::vector<std::pair<double, int> > v;
  int mini;
  double best = target;
  for (int i = 0; i < matches; i++) {
    int ci = match.at(i);
    if ((ci == 40350) || (ci==73133) || (ci==73138)){ ACTS_INFO("Found a correct pair");}
    else{ exit(0);}
    //double s = scoreTriple(ai, bi, ci);
    double s = evaluateScore(ctx, ci, dp, xp, bap);
    v.push_back(std::make_pair(s, ci));
  }
  if (!matches) return 0;
  std::sort(v.begin(), v.end());
  double thres = v[0].first*4;
  for (int i = 0; i < v.size(); i++) {
    if (i >= 2 || v[i].first > thres) break;// && v[i].first > mid0) break; //i >= 1 and added < 1 gives better than old
    if (rev)
      triples.push_back(triple(v[i].second, bi, ai));
    else
      triples.push_back(triple(ai, bi, v[i].second));
  }

  // if (ctx.eventStore.add("triples", std::move(*triples)) == FW::ProcessCode::ABORT){
  //   ACTS_VERBOSE("Not able to store triples into eventStore in extendTripleLine");
  // }

  return 1;
}

int
FW::Topquarks::getIndex(FW::AlgorithmContext ctx, int&li, double x) const
{
  const int crude_steps = 1<<10;
  const std::array<std::pair<double, double>, 48>* crudeIndex_a = nullptr;
  const std::array<std::array<int, 48>, crude_steps>* crudeIndex = nullptr;
  const std::array<std::vector<double>, 48>* sorted_hits = nullptr;

  if (ctx.eventStore.get("crudeIndex", crudeIndex) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get crudeIndex from eventStore in getIndex");
  }
  if (ctx.eventStore.get("crudeIndex_a", crudeIndex_a) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get crudeIndex_a from eventStore in getIndex");
  }
  if (ctx.eventStore.get("sorted_hits", sorted_hits) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get sorted_hits from eventStore in getIndex");
  }
  // ACTS_VERBOSE("got objects from eventStore in getIndex ")

  int ci = x*crudeIndex_a->at(li).first+crudeIndex_a->at(li).second;
  // ACTS_VERBOSE("Value of ci in getIndex " << ci << "\n");
  ci = std::min(crude_steps-1, std::max(0, ci));
  int i = crudeIndex->at(li)[ci];
  // ACTS_VERBOSE("Value of i in getIndex " << i << "\n");
  // ACTS_VERBOSE("Value of x in getIndex " << x << "\n");
  while (x >= sorted_hits->at(li)[i]) i++;
  while (x < sorted_hits->at(li)[i-1]) i--;
  // ACTS_VERBOSE("Returning " << i);
  // for(double j: sorted_hits->at(li)) ACTS_VERBOSE(j << "\n");
  // return 0;
  return i;//max(0,min(i,int(sorted_hits[li].size())-1));
}

double
FW::Topquarks::scoreTriple(FW::AlgorithmContext ctx, int ai, int bi, int ci) const
{
  FW::point center;
  const std::vector<FW::point>* hits;
  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in getting hits from eventStore in scoreTriple");
  }
  // ACTS_VERBOSE("Got object in scoreTriple from eventStore");
  double radius;
  circle(hits->at(ai), hits->at(bi), hits->at(ci), center, radius);

  FW::point cb = hits->at(ci)-hits->at(bi);
  FW::point ba = hits->at(bi)-hits->at(ai);
  double ang_cb = asin(dist(cb.x, cb.y)*.5/radius)*2;
  double ang_ba = asin(dist(ba.x, ba.y)*.5/radius)*2;
  if (radius != radius || fabs(radius) > 1e50) {
    ang_cb = dist(cb.x, cb.y);
    ang_ba = dist(ba.x, ba.y);
  }
  if (ba.z*cb.z < 0) ang_ba *= -1;

  double x = ba.z ? (fabs(cb.z*ang_ba/ba.z-ang_cb))*radius : 1e9;
  double y = ang_cb ? (fabs(cb.z*ang_ba/ang_cb-ba.z)) : 1e9;
  double score = std::min(x, y);//, fabs(cb.z-ba.z*ang_cb/ang_ba)));
  // ACTS_VERBOSE("Scored triples. Returning Score " << score);
  return score;
}

int
FW::Topquarks::acceptTriple(FW::AlgorithmContext ctx, triple&t) const
{
  double A = log(scoreTriple(ctx, t.x, t.y, t.z)+1e-8);
  double B = log(scoreTripleDensity(ctx, t.x, t.y, t.z));
  double C = log(scoreTripleDensity(ctx, t.z, t.y, t.x));
  double L[3];
  double D = scoreTripleLogRadius_and_HitDir(ctx, t.x,t.y,t.z,L);
  double w[121] = {-13.215482291516638, -0.519368174205195, -0.6019168737814719, -0.400773825827796, -3.0689189279504614, -8.21987444638849, -1.773083608093787, -3.7271459966647913, -0.18753136282696767, -0.1700350202416788, -0.13325020734065293, -0.0712787103124509, 2.2365502305889295, -0.38264699613950004, -1.5996361946235698, 0.02602607302842127, -0.04090477074387659, -0.12800207114108786, -2.0616314706262706, 0.9350417490331662, -0.6313327964001432, 0.00830034532077729, -0.1021716887039019, 0.3719980432444666, 0.43818671158350325, 0.0338130884608543, 0.19225191422472998, -0.33226136371623366, -1.0631299954951279, -1.3290857109832128, 8.50340840417112, 4.489339976769724, -3.6866359902477703, -1.530677908510526, -0.3660375991432235, -0.2832850515900752, -0.003067393550643885, -0.06185860378584967, -0.004472073355177509, -0.034047188478205974, 0.056232208303619684, -0.09251101374546467, -0.3186456107148592, -0.011497406815609599, 0.0040898730087192275, 0.04166475101451824, 0.5313081554181062, 0.05691563704023761, 0.004054188315119864, 0.009440976068230187, 0.015452389083207108, 0.02857025202533131, -0.01788978369714811, -0.014820867725080077, -0.0032975179225221054, -0.2739810756530984, -0.209895536224461, -0.05555070596049059, -3.8393234795148445, -0.39189992715019867, 0.5302884318217037, -1.0560724732243318, 0.5808249742500916, 0.2085127159157602, -0.002879796716268462, -0.008289453513497825, -0.013327308424637882, 0.034516052559319284, 0.05612738574267425, -0.04698958101602463, 0.0007407605230615924, -0.015547995524776616, 0.06280040184070336, -0.056422842974113374, -0.02553695075115984, -0.030162351232030156, -0.216209409546151, 0.03852063554031595, -0.0693834129966963, -1.0570960495204662, 0.6811799884934827, 0.3386224510850844, -0.10244400357684635, -0.17437169642288441, 0.527447777429105, -0.0009197072806774356, -0.004512546596535816, -0.026048615023175962, -0.016165328699534447, -0.007957908851240184, -0.01677913671380496, 0.00448514125782629, -0.0164129789525374, -0.04792927265651915, 0.3459064488723725, 0.08305188504334206, -0.4177214300084773, -0.09227473501037928, 0.04508615512899353, -0.03988016215006392, 0.029600325479286028, -0.2533468783999991, -0.1438693183062194, -0.17942937900359165, -1.277174888294048, -0.12050721012197445, -1.306910361564254, -0.056617003726385146, -1.1681337555296898, 0.06259298498866638, -6.501290522349262, -10.841239719611016, 2.156866020752887, 1.3871764445557901, 4.945464722802966, -4.26463890108575, -1.510051189434741, -3.140021739172429, -5.693045331329942, 1.1610084032964856, 2.2204604560570425};
  double x[8] = {1,A,B,C,D,L[0],L[1],L[2]};
  int c = 0;
  double score = 0;
  for (int i = 0; i < 8; i++) {
    for (int j = i; j < 8; j++) {
      double a = x[i]*x[j];
      for (int k = j; k < 8; k++)
	score += a*x[k]*w[c++];
    }
  }
  return score > w[120];
}

// How many outliers do we expect to fit better than "ci" in the triple "ai", "bi" and "li"
double
FW::Topquarks::scoreTripleDensity(FW::AlgorithmContext ctx, int ai, int bi, int ci) const
{
  FW::point d, dp, xp, bap;
  const double density_eps = 1e-6;
  const std::vector<int>* metai = nullptr;
  const std::vector<FW::point>* polar = nullptr;
  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in getting metai from eventStore in scoreTripleDensity");
  }
  if (ctx.eventStore.get("polar", polar) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in getting polar from eventStore in scoreTripleDensity");
  }
  if (prepareTripleScore(ctx, ai, bi, metai->at(ci), d, dp, xp, bap, polar->at(ci))) return 1e9;
  double s = evaluateScore(ctx, ci, dp, xp, bap);
  // ACTS_VERBOSE("Got score" << s)
  s = getDensity3(ctx, dp, xp, s, metai->at(ci));
  // ACTS_VERBOSE("Got density" << s);
  return s+density_eps;
}

//Prepare ellipse equation of collision between line extrapolated through hits with id "ai" and "bi" and layer "li". Return collision coordinate "d", in polar coordinates "dp", ellipse stretching "xp", and direction of hit in polar coordnates "bap". "target" describes the layer, possibly corrected for a single point we are evaluating a helix quadruple
int
FW::Topquarks::prepareTripleScore(FW::AlgorithmContext ctx, int ai, int bi, int li, point&d, point&dp, point&xp, point&bap, point target) const
{
  double slack = 1.00; //No slack
  const std::array<Layer, 48>* layer = nullptr;
  const std::vector<FW::point>* hits = nullptr;

  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting hits from eventStore in prepareTripleScore");
  }
  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting layer from eventStore in prepareTripleScore");
  }
  // ACTS_VERBOSE("Got layer and hits from eventStore in prepareTripleScore");

  const Layer&l = layer->at(li);
  const FW::point&a = hits->at(ai), &b = hits->at(bi);
  const int Tube = 0, Disc = 1;
  double stretch = 0.02;

  // ACTS_VERBOSE("no error0");

  const FW::point ba = b-a;
  if (l.type == Tube) {
    double vv = ba.x*ba.x+ba.y*ba.y;
    double pv = ba.x*a.x+ba.y*a.y;
    double pp = a.x*a.x+a.y*a.y;
    double RR = target.x*target.x;
    double sq = pv*pv-vv*(pp-RR);
    if (sq < 0) return -1;
    // ACTS_VERBOSE("no error1");

    double t = (-pv+sqrt(sq))/vv;
    if (t < 0 || !vv) return -1;
    d.x = a.x+ba.x*t;
    d.y = a.y+ba.y*t;
    d.z = a.z+ba.z*t;
    // ACTS_VERBOSE("no error2");

    if (d.z < l.minz*slack || d.z > l.maxz*slack) return -1;
    // ACTS_VERBOSE("no error3");

    dp = FW::point(dist(d.x,d.y),atan2(d.y,d.x),d.z);
    // ACTS_VERBOSE("no error4");

    xp = FW::point(0, -dp.x*(ba.x*ba.x+ba.y*ba.y), ba.z);
    // ACTS_VERBOSE("xp.x is 0 and xp.z is " << ba.z);
    bap = FW::point(ba.x*d.x+ba.y*d.y, d.x*ba.y-d.y*ba.x, ba.z*dp.x);
    // ACTS_VERBOSE("no error5");

    bap = bap*(1./bap.x);
  }
  else if (l.type == Disc) {
    // ACTS_VERBOSE("no error1");
    double t = (target.z-a.z)/ba.z;
    if (t < 0 || !ba.z) return -1;
    d.x = a.x+ba.x*t;
    d.y = a.y+ba.y*t;
    d.z = a.z+ba.z*t;
    // ACTS_VERBOSE("no error2");

    dp = FW::point(dist(d.x,d.y),atan2(d.y,d.x),d.z);
    // ACTS_VERBOSE("no error3");

    if (dp.x < l.minr*(1./slack) || dp.x > l.maxr*slack) return -1;
    // ACTS_VERBOSE("no error4");

    xp = FW::point(ba.x*d.y-ba.y*d.x, d.x*ba.x+d.y*ba.y, 0);
    bap = FW::point(ba.x*d.x+ba.y*d.y, d.x*ba.y-d.y*ba.x, ba.z*dp.x);
    // ACTS_VERBOSE("no error5");

    bap = bap*(1./bap.z);
  }
  double xp2 = xp.x*xp.x+xp.y*xp.y+xp.z*xp.z;
  // ACTS_VERBOSE("no error6");
  if (xp2)
    xp = xp*sqrt((1-stretch)/xp2);

  return 0;
}

int
FW::Topquarks::prepareTripleScoreDefault(FW::AlgorithmContext ctx, int ai, int bi, int li, point&d, point&dp, point&xp, point&bap) const
{
  const std::array<Layer, 48>* layer = nullptr;
  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting layer from eventStore in prepareTripleScoreDefault");
  }
  // ACTS_VERBOSE("Got layer from eventStore in prepareTripleScoreDefault");
  const Layer&l = layer->at(li);
  FW::point target(l.avgr, 0, l.avgz);
  return prepareTripleScore(ctx, ai, bi, li, d, dp, xp, bap, target);
}

int
FW::Topquarks::prepareQuadrupleScore(FW::AlgorithmContext ctx, int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&dirp, point target, double sign = 1) const
{
  const int Tube = 0, Disc = 1;
  const std::array<Layer, 48>* layer = nullptr;
  const std::vector<FW::point>* hits = nullptr;
  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting hits from eventStore in prepareQuadrupleScore");
  }
  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in getting layer from eventStore in prepareQuadrupleScore");
  }
  const Layer&l = layer->at(li);

  FW::point p;
  double r, ir;

  const FW::point c = hits->at(ci);
  const FW::point cb = hits->at(ci)-hits->at(bi);//c-b;
  double ang_cb;

  if (0) {
    const FW::point a = hits->at(ai), b = hits->at(bi), c = hits->at(ci);

    //TODO: test if has bieffects
    if (l.type == Disc && (c.z-b.z)*(target.z-c.z) < 0) return -1;

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

    ang_cb = asin(dist(cb.x, cb.y)*.5*ir)*2;
  }

  if (1) { //Take into account varying magnetic field strength
    const FW::point a = hits->at(ai), b = hits->at(bi), c = hits->at(ci);
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

  double slack = 1.00;

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

    FW::point dc = d-c;
    double A = dist(dc.x, dc.y);
    double B = A*.5*ir;
    double ang_dc = asin(B)*2;
    if (dc.x*cb.x+dc.y*cb.y < 0) ang_cb *= -1;

    d.z = c.z+cb.z*ang_dc/ang_cb;

    if (!(d.z > l.minz*slack && d.z < l.maxz*slack)) return -1;

    FW::point dir;
    double s_ = target.x/pp, t_ = s_*(1-s)/t;
    dir.x = p.x*s_+p.y*t_;
    dir.y = p.y*s_-p.x*t_;
    dir.z = (dc.x*dir.x+dc.y*dir.y)*ir*cb.z/(ang_cb*A*sqrt(1-B*B));

    dp = point(dist(d.x,d.y), atan2(d.y, d.x), d.z);
    dirp = point(d.x*dir.x+d.y*dir.y, d.x*dir.y-d.y*dir.x, dir.z*dp.x);

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

    dirp = dirp*(1./dirp.z);
  }
  return 0;
}

double
FW::Topquarks::field(double z) const
{
  z *= 1./2750;
  double z2 = z*z;
  return 1.002-z*3e-2-z2*(0.55-0.3*(1-z2));
}

double
FW::Topquarks::evaluateScore(FW::AlgorithmContext ctx, int ci, point&dp, point&xp, point&bap) const
{
  const std::array<Layer, 48>* layer = nullptr;
  const std::vector<int>* metai = nullptr;
  const std::vector<FW::point>* polar = nullptr;
  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get layer from eventStore in evaluateScore");
  }
  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get metai from eventStore in evaluateScore");
  }
  if (ctx.eventStore.get("polar", polar) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get polar from eventStore in evaluateScore");
  }
  const int Disc = 1;
  const FW::point&r = polar->at(ci);
  // ACTS_VERBOSE("Check0");
  FW::point err = r-dp;
  // ACTS_VERBOSE("Check0");
  if (err.y > M_PI) err.y -= M_PI*2;
  if (err.y <-M_PI) err.y += M_PI*2;
  err.y *= dp.x;
  // ACTS_VERBOSE("Check2");
  err = err-bap*(layer->at(metai->at(ci)).type == Disc ? err.z : err.x);
  // ACTS_VERBOSE("Check3");
  double r2 = err*err-std::pow(err*xp, 2);
  // ACTS_VERBOSE("Check4");
  return r2;
}

int
FW::Topquarks::prepareQuadrupleScoreDefault(FW::AlgorithmContext ctx, int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&bap, double sign = 1) const
{
  const std::array<Layer, 48>* layer = nullptr;
  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in getting layer from eventStore in prepareQuadrupleScoreDefault");
  }
  const Layer&l = layer->at(li);
  FW::point target(l.avgr, 0, l.avgz);
  return prepareQuadrupleScore(ctx, ai, bi, ci, li, d, dp, xp, bap, target, sign);
}

//Return features for logistic regression of a triple. L has cell's data angle errors, return logarithm of inverse radius of helix
double
FW::Topquarks::scoreTripleLogRadius_and_HitDir(FW::AlgorithmContext ctx, int ai, int bi, int ci, double*L) const
{
  const std::vector<FW::point>* hits = nullptr;
  const std::array<std::array<FW::point, 2>, 200000>* hit_dir = nullptr;
  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in getting hits from eventStore in scoreTripleLogRadius_and_HitDir");
  }
  if (ctx.eventStore.get("hit_dir", hit_dir) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in getting hit_dir from eventStore in scoreTripleLogRadius_and_HitDir");
  }
  point a = hits->at(ai), b = hits->at(bi), c = hits->at(ci);
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

  for (int k = 0; k < 3; k++) {
    int di = k ? k==2 ? ci : bi : ai;
    double rx = hits->at(di).x-p.x, ry = hits->at(di).y-p.y;

    point ca = hits->at(ci)-hits->at(ai);
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
    L[k] = acos(std::max(fabs(hit_dir->at(di).at(0)*dir),
		    fabs(hit_dir->at(di).at(1)*dir))/dist0(dir));
  }
  return log(ir);
}

double
FW::Topquarks::dist0(point&p) const
{
  return sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
}

std::vector<FW::triple>
FW::Topquarks::findTriplesDuplicates(FW::AlgorithmContext ctx, std::array<PolarModule, 48> &mod, std::array<int, 200000> &match) const{
  ACTS_VERBOSE("In findTriplesDuplicates");
  std::vector<std::pair<int, int> > pairs_temp = findDuplicates(ctx, mod, match);
  std::vector<FW::triple> triples_temp;
  ACTS_INFO("addDuplicateTriples loop size " << pairs_temp.size());
  for (auto&p : pairs_temp) {
    addDuplicateTriples(ctx, triples_temp, mod, match, p.first, p.second);
  }
  ACTS_INFO("Found " << triples_temp.size() << " duplicate triple candidates");
  return triples_temp;
}

//(Try to) find all duplicates in the dataset by looking for point close together given cells' data velocity direction
std::vector<std::pair<int, int> >
FW::Topquarks::findDuplicates(FW::AlgorithmContext ctx, std::array<PolarModule, 48> &mod, std::array<int, 200000> &match) const
{
  ACTS_VERBOSE("In findDuplicates");
  std::vector<std::pair<int, int> > pairs_temp;
  const std::array<std::vector<int>, 48>* tube = nullptr;
  const std::vector<int>* metaz = nullptr;
  const std::vector<FW::point>* hits = nullptr;
  const std::vector<FW::point>* polar = nullptr;
  const std::array<std::array<point, 2>, 200000>* hit_dir = nullptr;
  const std::array<Layer, 48>* layer = nullptr;
  std::map<std::pair<int, int>, double>* v_duplicates = new std::map<std::pair<int, int>, double>();

  if (ctx.eventStore.get("tube", tube) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Not able to load tube from eventStore in findDuplicates");
  }
  if (ctx.eventStore.get("metaz", metaz) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Not able to load metaz from eventStore in findDuplicates");
  }
  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Not able to load hits from eventStore in findDuplicates");
  }
  if (ctx.eventStore.get("polar", polar) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Not able to load tube from eventStore in findDuplicates");
  }
  if (ctx.eventStore.get("hit_dir", hit_dir) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Not able to load hit_dir from eventStore in findDuplicates");
  }
  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Not able to load layer from eventStore in findDuplicates");
  }
  // ACTS_VERBOSE("Got objects from eventStore in findDuplicates");
  const int Tube=0, Disc = 1;
  for (int li = 0; li < 48; li++) {
    // ACTS_VERBOSE("Running " << li << " out of 48");
    FW::point d, dp, xp, bap;
    double target = 0.05;

    int fails = 0;
    std::set<long long> found;
    ACTS_VERBOSE("Inner loop size " << tube->at(li).size());
    for (int hit_id : tube->at(li)) {
      // ACTS_VERBOSE("Running inner loop " << hit_id << " out of " << tube->at(li).size());
      for (int k = 0; k < 2; k++) {
        d = hits->at(hit_id);
        dp = polar->at(hit_id);
        FW::point dir = hit_dir->at(hit_id).at(k);
        xp = FW::point(0,0,0);
        FW::point bap = topolar(dir, d, dp);
        // ACTS_VERBOSE("Made bap");
        bap = bap*(1./(layer->at(li).type == Disc ? bap.z : bap.x));
        double tt = findDensity(ctx, dp, xp, target, li);
        int matches = mod[li].getNear(ctx, dp, xp, bap, tt, match, li);
        for (int i = 0; i < matches; i++) {
        	int di = match[i];
        	if (metaz->at(di) == metaz->at(hit_id)) continue;
        	auto p = std::make_pair(hit_id, di);
        	double s = getDensity3(ctx, dp, xp, evaluateScore(ctx, di, dp, xp, bap), li);
          // ACTS_INFO()
        	if (!v_duplicates->count(p) || s < (*v_duplicates)[p])
        	  (*v_duplicates)[p] = s;
        }
      }
    }
    ACTS_VERBOSE("Second loop size " << v_duplicates->size());
    for (auto&i : *v_duplicates) {
      int a = i.first.first, b = i.first.second;
      if (a < b || !v_duplicates->count(std::make_pair(b,a)))
        pairs_temp.push_back(std::make_pair(a, b));
    }
  }
  return pairs_temp;
}

void
FW::Topquarks::addDuplicateTriples(FW::AlgorithmContext ctx, std::vector<triple>&triples, std::array<PolarModule, 48> &mod, std::array<int, 200000> &match,int ai, int bi) const
{
  const std::vector<int>* metai = nullptr;
  const std::array<std::array<int, 48>, 48>* next_layer = nullptr;
  const std::vector<FW::point>* hits = nullptr;
  // ACTS_INFO("In addDuplicateTriples");
  if (ctx.eventStore.get("next_layer", next_layer) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in loading next_layer from eventStore");
    // return FW::ProcessCode::ABORT;
  }
  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in loading metai from eventStore");
    // return FW::ProcessCode::ABORT;
  }
  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_VERBOSE("Error in loading hits from eventStore");
    // return FW::ProcessCode::ABORT;
  }
  std::vector<int> s[48];
  for (int li = 0; li < 48; li++) {
    int itr = li + 1;
    // ACTS_VERBOSE("First Iteration " << itr << " out of 48 in addDuplicateTriples");
    if (next_layer->at(li)[metai->at(ai)] < adj_thres/2 &&
	next_layer->at(metai->at(ai))[li] < adj_thres/2) continue;
    if (li == metai->at(ai)) continue;

    point d, dp, xp, bap;
    if (prepareTripleScoreDefault(ctx, ai, bi, li, d, dp, xp, bap) &&
	prepareTripleScoreDefault(ctx, bi, ai, li, d, dp, xp, bap)) continue;

    xp = normalize(xp)*0.99;
    //xp = point(0,0,0);
    const double target = 100;
    double mid = findDensity(ctx, dp, xp, target, li);
    int matches = mod[li].getNear(ctx, dp, xp, bap, mid, match, li);

    double best = target;
    std::vector<std::pair<double, int> > v;
    for (int i = 0; i < matches; i++) {
      int ci = match.at(i);
      double s = evaluateScore(ctx, ci, dp, xp, bap);
      v.push_back(std::make_pair(s, ci));
    }
    std::sort(v.begin(), v.end());
    for (int i = 0; i < v.size(); i++) {
      if (i >= target) break;
      s[li].push_back(v[i].second);
      //triples.push_back(triple(ai, bi, v[i].second));
    }
    //break;
  }

  for (int la = 0; la < 48; la++) {
    int itr2 = la + 1;
    // ACTS_INFO("Second Iteration " << itr2 << " out of 48 in addDuplicateTriples");

    for (int ci : s[la]) {
      int xi = ai, yi = bi, zi = ci;
      if ((hits->at(zi).z-hits->at(yi).z)*(hits->at(yi).z-hits->at(xi).z) < 0) std::swap(xi, yi);
      if (metai->at(zi) < metai->at(xi)) std::swap(xi, zi);
      for (int li = 0; li < 48; li++) {
      	if (li >= metai->at(xi) && li <= metai->at(zi)) continue;
      	if (next_layer->at(li)[metai->at(zi)] < adj_thres/2 &&
      	    next_layer->at(metai->at(xi))[li] < adj_thres/2) continue;

      	int xi_ = xi, zi_ = zi;
      	if (li < metai->at(xi_)) std::swap(xi_, zi_);

      	point d, dp, xp, bap;
      	if (prepareQuadrupleScoreDefault(ctx, xi_, yi, zi_, li, d, dp, xp, bap)) continue;

      	double target2 = 0.5;
      	double tt = findDensity(ctx, dp, xp, target2, li);
      	int matches = mod[li].getNear(ctx, dp, xp, bap, tt, match, li);
      	std::vector<std::pair<double, int> > v;
      	for (int i = 0; i < matches; i++) {
      	  int di = match[i];
      	  double s = evaluateScore(ctx, di, dp, xp, bap);
      	  v.push_back(std::make_pair(s, di));
      	}
      	std::sort(v.begin(), v.end());
      	for (int i = 0; i < v.size() && i < target2; i++) {
      	  triple t(ai, ci, v[i].second);
      	  if (metai->at(t.x) > metai->at(t.y)) std::swap(t.x,t.y);
      	  if (metai->at(t.y) > metai->at(t.z)) std::swap(t.y,t.z);
      	  if (metai->at(t.x) > metai->at(t.y)) std::swap(t.x,t.y);
      	  // assert(metai->at(t.x) < metai->at(t.y) && metai->at(t.y) < metai->at(t.z), "Wrong order!");
      	  triples.push_back(t);
      	}
      }
    }
  }
}

//Convert "dir" to polar (really cylindrical coordinates) coordinates, suffix p  (as in "refp") usually means polar coodinates throughout the code
FW::point
FW::Topquarks::topolar(point&dir, point&ref, point&refp) const
{
  return FW::point(ref.x*dir.x+ref.y*dir.y, ref.x*dir.y-ref.y*dir.x, dir.z*refp.x);
}

std::vector<FW::triple>
FW::Topquarks::pruneTriples(FW::AlgorithmContext ctx, std::vector<FW::triple>&triples) const
{
  int start_list_len = 100;
  triple start_list[100] = {{0, 1, 2}, {11, 12, 13}, {4, 5, 6}, {0, 1, 11}, {0, 1, 4}, {0, 4, 5}, {0, 11, 12}, {18, 19, 20}, {1, 2, 3}, {5, 6, 7}, {12, 13, 14}, {13, 14, 15}, {6, 7, 8}, {2, 3, 18}, {3, 18, 19}, {19, 20, 21}, {0, 1, 3}, {0, 2, 3}, {20, 21, 22}, {0, 1, 5}, {0, 1, 12}, {1, 4, 5}, {1, 2, 18}, {11, 18, 19}, {1, 11, 12}, {7, 8, 9}, {4, 18, 19}, {14, 15, 16}, {0, 4, 6}, {18, 19, 24}, {18, 19, 36}, {2, 18, 19}, {0, 18, 19}, {1, 18, 19}, {21, 22, 23}, {11, 12, 18}, {4, 5, 18}, {0, 11, 13}, {0, 1, 18}, {24, 26, 28}, {1, 2, 11}, {1, 2, 4}, {36, 38, 40}, {15, 16, 17}, {8, 9, 10}, {5, 18, 19}, {18, 24, 26}, {18, 36, 38}, {38, 40, 42}, {12, 18, 19}, {0, 12, 13}, {40, 42, 44}, {28, 30, 32}, {18, 20, 21}, {0, 2, 18}, {26, 28, 30}, {4, 5, 7}, {19, 20, 24}, {11, 12, 14}, {19, 20, 36}, {18, 19, 21}, {6, 18, 19}, {13, 18, 19}, {2, 11, 18}, {0, 5, 6}, {0, 2, 11}, {12, 13, 18}, {19, 36, 38}, {4, 6, 7}, {19, 24, 26}, {2, 4, 18}, {0, 2, 4}, {5, 6, 18}, {3, 19, 20}, {11, 13, 14}, {11, 12, 36}, {2, 3, 19}, {1, 3, 18}, {19, 22, 23}, {2, 18, 24}, {4, 5, 24}, {7, 18, 19}, {20, 22, 23}, {24, 26, 29}, {11, 18, 36}, {19, 20, 22}, {20, 21, 37}, {15, 18, 19}, {2, 18, 36}, {14, 18, 19}, {4, 18, 24}, {8, 18, 19}, {3, 18, 20}, {36, 38, 41}, {20, 21, 25}, {9, 18, 19}, {6, 7, 26}, {13, 14, 36}, {42, 44, 46}, {30, 32, 34}};
  std::set<triple> start_list_set;
  for (int i = 0; i < start_list_len; i++) start_list_set.insert(start_list[i]);
  //const std::vector<FW::triple>* triples = nullptr;
  // if (ctx.eventStore.get("triples", triples) == FW::ProcessCode::ABORT){
  //   ACTS_WARNING("Not able to get triples from eventStore in prunetriples");
  // }
  std::sort(triples.begin(), triples.end());
  for (int i = 0; i < triples.size(); i++) {
    triple&t = triples.at(i);
    int ok = 1;
    if (i && t.x == abs(triples.at(i-1).x) && t.y == triples.at(i-1).y && t.z == triples.at(i-1).z) ok = 0;
    else {
      point p = scoreTripleHitDir(t);
      double f = 4;
      if (metai->at(t.x) < 18) p.x *= f;
      if (metai->at(t.y) < 18) p.y *= f;
      if (metai->at(t.z) < 18) p.z *= f;
      double score = p.x+p.y+p.z;
      if (score >= 3.1) ok = 0;
    }
    if (0) {//Not used, only keep triples that are in the beginning of paths
      triple ii(metai->at(t.x), metai->at(t.y), metai->at(t.z));
      if (!start_list_set.count(ii)) ok = 0;
    }
    if (!ok) {
      t.x *= -1;
    }
  }
  for (int i = 0; i < triples.size(); i++) {
    triple&t = triples.at(i);
    if (t.x < 0) {
      std::swap(t, triples.at(triples.size()-1));
      triples.pop_back();
      i--;
    }
  }
  return triples;
}

FW::point
FW::Topquarks::scoreTripleHitDir(triple t) const
{
  point a = hits->at(t.x), b = hits->at(t.y), c = hits->at(t.z);
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

  point ca = hits->at(t.z)-hits->at(t.x);
  double ang_ca = asin(dist(ca.x, ca.y)*.5*ir)*2;

  double ret[3];
  for (int i = 0; i < 3; i++) {
    int di = i ? i==2 ? t.z : t.y : t.x;

    double rx = hits->at(di).x-p.x, ry = hits->at(di).y-p.y;
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
    double angle = 1e9;
    for (int k = 0; k < 2; k++) {
      point&hd = hit_dir->at(di).at(k);
      double angle_ = acos(fabs(dir*hd)/dist0(dir));
      angle = std::min(angle, angle_);
    }
    ret[i] = angle;
  }
  return point(ret[0], ret[1], ret[2]);
}

void
FW::Topquarks::storeInEventStore(AlgorithmContext ctx) const
{
  int failed = 0;
  if (ctx.eventStore.add("hits", std::move(hits)) == FW::ProcessCode::ABORT) failed++;
  if (ctx.eventStore.add("layer", std::move(layer)) == FW::ProcessCode::ABORT) failed++;
  if (ctx.eventStore.add("metai", std::move(metai)) == FW::ProcessCode::ABORT) failed++;
  if (ctx.eventStore.add("metaz", std::move(metaz)) == FW::ProcessCode::ABORT) failed++;
  if (ctx.eventStore.add("polar", std::move(hits)) == FW::ProcessCode::ABORT) failed++;
  if (ctx.eventStore.add("disc_z", std::move(disc_z)) == FW::ProcessCode::ABORT) failed++;
  //if (ctx.eventStore.add("assignment", std::move(assignment)) == FW::ProcessCode::ABORT) failed++;
}

std::vector<std::vector<int> >
FW::Topquarks::findPaths(FW::AlgorithmContext ctx, std::vector<FW::triple>&triples, std::array<PolarModule, 48> *mod, std::array<int, 200000> &match) const
{
  std::vector<std::vector<int> > paths;
  const std::vector<int>* metai = nullptr;

  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting metai from eventStore in findPaths");
  }

  paths.reserve(triples.size());
  for (auto p : triples) {
    // load(triples.size(), "Find paths", 0);
    std::vector<int> v;
    int ai = p.x, bi = p.y, ci = p.z, misses = 0;

    for (int li = metai->at(ai)-1; li >= 0; li--) {
      if (next_layer[li][metai->at(ai)] < adj_thres) continue;
      int di = extend3(ctx, mod, match, ci, bi, ai, li, 0.5);
      if (di > 0) {
	v.push_back(di);
	ci = bi;
	bi = ai;
	ai = di;
	misses = 0;
      } else if (di == -1) misses++;
      if (misses == 1) break;
    }
    std::reverse(v.begin(), v.end());

    ai = p.x, bi = p.y, ci = p.z;

    v.push_back(ai), v.push_back(bi), v.push_back(ci);

    misses = 0;
    for (int li = metai->at(ci)+1; li < 48; li++) {
      if (next_layer[metai->at(ci)][li] < adj_thres) continue;
      int di = extend3(ctx, mod, match, ai, bi, ci, li, 0.5);
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
    paths.push_back(v);
  }
  return paths;
}

int
FW::Topquarks::extend3(FW::AlgorithmContext ctx, std::array<PolarModule, 48> *mod, std::array<int, 200000> &match, int ai, int bi, int ci, int li, double target = 0.5) const {
  point d, dp, xp, bap;
  if (prepareQuadrupleScoreDefault(ctx, ai,bi,ci,li,d,dp,xp,bap)) return -2;

  //double mins = 0.5;//d*d*1e-4;
  double tt = 400;//findDensity(dp, xp, target, li);
  double mins = 1e9;//target;//1e9;
  // const std::array<int, 200000>* match = nullptr;
  // if (ctx.eventStore.get("match", match) == FW::ProcessCode::ABORT){
  //   ACTS_WARNING("Error in getting match from eventStore in extend3");
  // }
  double fac = getDensity3(ctx, dp, xp, tt, li)/tt;
  tt = target/fac;
  int matches = mod->at(li).getNear(ctx, dp, xp, bap, tt, match, li);
  int mini = -1;
  for (int i = 0; i < matches; i++) {
    int ti = match.at(i);
    // if (assignment[ti]) continue;
    double s = evaluateScore(ctx, ti,dp,xp,bap)*fac;
    if (s < mins) {
      mins = s;
      mini = ti;
    }
  }
  return mini;
}

std::vector<std::vector<int> >
FW::Topquarks::addDuplicates(FW::AlgorithmContext ctx, std::vector<std::vector<int> >&paths, std::array<PolarModule, 48> &mod, std::array<int, 200000> &match) const
{
  std::vector<std::vector<int> > extended;
  const std::vector<int>* metai = nullptr;
  const std::vector<point>* meta = nullptr;
  const std::vector<int>* metaz = nullptr;
  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting metai from eventStore");
  }
  if (ctx.eventStore.get("meta", meta) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting metai from eventStore");
  }
  if (ctx.eventStore.get("metaz", metaz) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting metai from eventStore");
  }
  for (auto&path : paths) {
    // load(paths.size(), "Add duplicates");
    if (path.size() < 3) continue;

    std::vector<int> ext;
    for (int i = 0; i < path.size(); i++) {
      ext.push_back(path[i]);
      int ai, bi, ci;
      if (i < path.size()-2) ai = path[i], bi = path[i+1], ci = path[i+2];
      else if (i == path.size()-1) ai = path[i-2], bi = path[i-1], ci = path[i];
      else ai = path[i-1], bi = path[i], ci = path[i+1];

      int li = metai->at(path[i]);
      FW::point d, dp, xp, bap;

      //Add on average about "target" outliers to each path
      const double target = 0.1;

      if (prepareDuplicateScore(ctx, ai, bi, ci, li, d, dp, xp, bap, match)) continue;

      double tt = findDensity(ctx, dp, xp, target, li);

      int pi = path[i];
      int matches = mod[li].getNear(ctx, dp, xp, bap, tt, match, li);//target/fac

      std::map<int, std::pair<double, int> > mins;

      for (int i = 0; i < matches; i++) {
      	int di = match[i];
      	double s = evaluateScore(ctx, di, dp, xp, bap);//*fac;
      	if (meta->at(di).z != meta->at(pi).z) {
      	  int zi = metaz->at(di);
      	  if (!mins.count(zi) || s < mins[zi].first) {
      	    mins[zi] = std::make_pair(s, di);
      	  }
      	}
      }
      //vector<pair<double, int> > v;
      for (auto&p : mins)
      	ext.push_back(p.second.second);
      /*sort(v.begin(), v.end());
      for (auto&p : v)
      ext.push_back(p.second);*/
    }
    ext.shrink_to_fit();
    path.clear();
    path.shrink_to_fit();
    extended.push_back(ext);
  }
  return extended;
}

int
FW::Topquarks::prepareDuplicateScore(FW::AlgorithmContext ctx, int ai, int bi, int ci, int li, FW::point&d, FW::point&dp, FW::point&xp, FW::point&dirp, std::array<int, 200000> &match) const
{
  const std::array<Layer, 48>* layer = nullptr;
  const std::vector<FW::point>* hits = nullptr;
  const std::vector<int>* metai = nullptr;
  const std::vector<FW::point>* polar = nullptr;
  // const std::vector<FW::point>* meta = nullptr;
  // const std::vector<<int>* metaz = nullptr;
  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting layer from eventStore in prepareDuplicateScore");
  }
  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting metai from eventStore in prepareDuplicateScore");
  }
  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting hits from eventStore in prepareDuplicateScore");
  }
  if (ctx.eventStore.get("polar", polar) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting polar from eventStore in prepareDuplicateScore");
  }
  // if (ctx.eventStore.get("meta", meta) == FW::ProcessCode::ABORT){
  //   ACTS_WARNING("Error in getting meta from eventStore in prepareDuplicateScore");
  // }
  // if (ctx.eventStore.get("metaz", metaz) == FW::ProcessCode::ABORT){
  //   ACTS_WARNING("Error in getting metaz from eventStore in prepareDuplicateScore");
  // }
  const Layer&l = layer->at(li);
  FW::point a = hits->at(ai), b = hits->at(bi), c = hits->at(ci);
  //Find circle with center p, radius r, going through a, b, and c (in xy plane)
  double ax = a.x-c.x, ay = a.y-c.y, bx = b.x-c.x, by = b.y-c.y;
  double aa = ax*ax+ay*ay, bb = bx*bx+by*by;
  double idet = .5/(ax*by-ay*bx);
  FW::point p;
  p.x = (aa*by-bb*ay)*idet;
  p.y = (ax*bb-bx*aa)*idet;
  p.z = 0;
  double r = dist(p.x, p.y), ir = 1./r;
  p.x += c.x;
  p.y += c.y;

  int di = -1;
  if (metai->at(ai) == li) di = ai;
  else if (metai->at(bi) == li) di = bi;
  else if (metai->at(ci) == li) di = ci;
  else {
    // cout << "prepareDuplocateScore given layeri not corresponding to ai, bi or ci" << endl;
    return -1;
  }
  d = hits->at(di);
  dp = polar->at(di);

  double rx = hits->at(di).x-p.x, ry = hits->at(di).y-p.y;

  //TODO: do with respect to nearest circle arc, not ca
  point ca = hits->at(ci)-hits->at(ai);
  double ang_ca = asin(dist(ca.x, ca.y)*.5*ir)*2;
  double cross = rx*ca.y-ry*ca.x;

  FW::point dir;
  if (ir) {
    dir.x =-ry*ang_ca;
    dir.y = rx*ang_ca;
    dir.z = ca.z;
    if (cross < 0) dir.z *= -1;
  } else {
    dir = ca;
  }

  /*
  //dir = truth_mom[di];
  {
    point dir2 = dir*(-1./sqrt(dir*dir));
    point mom2 = truth_mom[di]*(1./dist(truth_mom[di]));
    if (dir2.z * mom2.z < 0) dir2 = dir2*-1;
    cout << dir2 << endl;
    cout << mom2 << endl << endl;
  }
  */

  xp = point(0,0,0);
  dirp = point(d.x*dir.x+d.y*dir.y, d.x*dir.y-d.y*dir.x, dir.z*dp.x);
  //cout << dirp << endl; //dirp.x = l.avgr
  if (l.type == Tube)
    dirp = dirp.x ? dirp*(1./dirp.x) : point(1,0,0);
  else
    dirp = dirp.z ? dirp*(1./dirp.z) : point(0,0,1);

  return 0;
}

std::vector<std::vector<int> >
FW::Topquarks::prunePaths(FW::AlgorithmContext ctx, std::vector<std::vector<int> >&paths, std::array<int, 200000> &match) const
{
  std::vector<std::vector<int> > pruned;

  std::vector<std::vector<int> > hash_list;
  hash_list.resize(1<<17);
  int c = 0;
  double thres = 2.5;
  for (auto&path : paths) {
    // load(paths.size(), "Pruning paths", 0);
    while (path.size() >= 4 && scoreQuadrupleDensity(ctx, path[3], path[2], path[1], path[0]) > thres) {
      path.erase(path.begin());
    }
    /*while (path.size() >= 4 && scoreQuadrupleDensity(path[path.size()-4], path[path.size()-3], path[path.size()-2], path[path.size()-1]) > 1e2) {
      path.pop_back();
    }*/
    double s = -scorepathDensity(ctx, path, match);
    if (s < 10 && path.size() >= 3) {
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
      	hash_list[h].push_back(pruned.size());
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

double
FW::Topquarks::scoreQuadrupleDensity(FW::AlgorithmContext ctx, int ai, int bi, int ci, int di) const
{
  const std::vector<int>* metai = nullptr;
  const std::vector<FW::point>* polar = nullptr;
  const double density_eps = 1e-6;

  if (ctx.eventStore.get("polar", polar) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get polar from eventStore in scoreQuadrupleDensity");
  }
  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get metai from eventStore in scoreQuadrupleDensity");
  }
  point d, dp, xp, bap;
  if (prepareQuadrupleScore(ctx, ai, bi, ci, metai->at(di), d, dp, xp, bap, polar->at(di))) return 1e9;
  double s = evaluateScore(ctx, di, dp, xp, bap);
  s = getDensity3(ctx, dp, xp, s, metai->at(di));
  return s+density_eps;

}

// Score a full path by looking at the probability that it could happen by outliers alone
// Multiply average number of outliers gotten from
// - The first triple
// - All consecutive quadruples
// - Duplicates
// This function is very important for score, and full of tuning opportunities
double
FW::Topquarks::scorepathDensity(FW::AlgorithmContext ctx, std::vector<int>&path, std::array<int, 200000> &match) const
{
  static std::vector<int> u(48);
  u.resize(0);
  //cout << u.capacity() << endl;

  int last_metai = -1;
  for (int i : path) {
    if (i <= 0) continue;
    int mi = metai->at(i);
    if (mi != last_metai) {
      //if (last_metai != -1 && next_layer[last_metai][mi] < adj_thres) return -1e4;
      u.push_back(i);
      last_metai = mi;
    }
  }
  if (u.size() < 3) return -0.1;
  double prod = 1;
  double s = scoreTripleDensity(ctx, u[0], u[1], u[2]);//+1e-4;
  s = std::min(s, 1e4);
  prod *= s;

  //Hugely important tuning parameters
  double quad_off = 4;
  double dup_off = 10;
  double quad_max = 50;
  double quad_max1 = 10;
  double dup_max = 1;
  //if (u.size() >= 4)
  //  prod *= pow(scoreQuadrupleDensity(u[3], u[2], u[1], u[0]), 1);
  for (int i = 3; i < u.size(); i++) {
    double s = std::min(scoreQuadrupleDensity(ctx, u[i-3], u[i-2], u[i-1], u[i]), quad_max1);
    s *= std::min(scoreQuadrupleDensity(ctx, u[i], u[i-1], u[i-2], u[i-3]), quad_max1)*quad_off;
    //s = min(s, quad_max);
    prod *= s;
  }

  int j = 0;
  for (int i = 0; i < path.size(); i++) {
    if (path[i] <= 0) continue;
    while (metai->at(path[i]) != metai->at(u[j])) j++;
    if (path[i] == u[j]) continue;
    int k = std::max(j, 2);
    double s = pow(scoreDuplicateDensity(ctx, u[k-2], u[k-1], u[k], path[i], match), 0.92)*dup_off; // Really important magic constant of 0.92, I have no clue why
    s = std::min(s, dup_max);
    prod *= s;
  }

  return -prod;
}

double
FW::Topquarks::scoreDuplicateDensity(FW::AlgorithmContext ctx, int ai, int bi, int ci, int di, std::array<int, 200000> &match) const
{
  point d, dp, xp, bap;
  const std::vector<int>* metai = nullptr;
  const double density_eps = 1e-6;

  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to load metai from eventStore in scoreDuplicateDensity");
  }
  if (prepareDuplicateScore(ctx, ai, bi, ci, metai->at(di), d, dp, xp, bap, match)) return 1e9;
  double s = evaluateScore(ctx, di, dp, xp, bap);
  //cout << "S0: " << s << endl;
  s = getDensity3(ctx, dp, xp, s, metai->at(di));
  //if (!s) cout << "What?" << endl;
  return s+density_eps;
}

std::vector<std::vector<int> >
FW::Topquarks::findAssignment(FW::AlgorithmContext ctx, std::vector<std::vector<int> >&paths, std::array<PolarModule, 48> *mod, std::array<int, 200000> &match, int use_trash) const
{
  const std::vector<FW::point>* hits;

  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting hits from eventStore in findAssignment");
  }
  paths.insert(paths.begin(), std::vector<int>());
  std::map<int, int> map_assignment;
  myMap2 path_score(paths.size()+hits->size());
  std::vector<std::vector<std::pair<int, int> > > used_by;
  used_by.resize(hits->size());

  for (int i = 1; i < paths.size(); i++) {
    // load(paths.size()-1, "Finding assignment (1/2)", 0);
    double score = scorepathDensity(ctx, paths[i], match);
    path_score.add(i, score);
    for (int j = 0; j < paths[i].size(); j++)
      used_by[paths[i][j]].push_back(std::make_pair(i, j));
  }

  int total = hits->size()-1;
  for (int i = 1; i < hits->size(); i++) {
    used_by[i].shrink_to_fit();
    // if (!assignment[i] && used_by[i].empty()) total--;
  }
  /*
  int added = 0;
  for (int i = 1; i < hits.size(); i++)
    if (!used_by.count(i)) {
      lost += truth_weight[i];
      added++;
      vector<int> v = {i};
      //used_by[i].push_back(make_pair(paths.size(), 0));
      paths.push_back(v);
      double score = scorepath(paths[paths.size()-1]);
      path_score.add(paths.size()-1, score);
    }
    cout << "Added: " << added << endl;*/
  double lost = 0;
  // cout << "Unused: " << lost << endl;
  lost = 0;

  /*for (int i = 1; i < hits.size(); i++)
    if (assignment[i])
    assignment[i] += 1e8;*/

  std::vector<std::vector<int> > solution_paths;
  solution_paths.push_back(std::vector<int>());

  int last = 0, merged = 0, thrown = 0;
  while (path_score.notempty()) {
    std::pair<int, double> pop = path_score.pop();
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
      	// if (hit_id > 0)
      	  // cout << "ERROR: " << hit_id << endl;
      	continue;
      }
      // load(total, "Finding assignment (2/2)", 0);
      map_assignment[hit_id] = trash ? 0 : solution_paths.size();
      // if (trash) lost += truth_weight[hit_id];
    }
    if (!trash)
      solution_paths.push_back(paths[i]);

    for (int hit_id : paths[i]) {
      if (hit_id < 0) continue;

      std::vector<std::pair<int, int> >&used = used_by[hit_id];
      for (int k = 0; k < used.size(); k++) {
        	if (used[k].first == i) continue;
        	std::vector<int>&pi = paths[used[k].first];
        	if (pi.empty()) continue;
        	//cout << assignment[pi[used[k].second]] << ' '<< c << endl;
        	pi[used[k].second] *= -1;
        	double score = scorepathDensity(ctx, pi, match);
        	path_score.update(used[k].first, score);
      }
    }
    paths[i].clear();
    paths[i].shrink_to_fit();
  }
  return solution_paths;
}

std::vector<std::vector<int> >
FW::Topquarks::extendPaths(FW::AlgorithmContext ctx, std::vector<std::vector<int> >&paths, std::array<PolarModule, 48> *mod, std::array<int, 200000>& match) const
{

  // initDensity3(ctx);

  // PolarModule mod[48];
  // for (int i = 0; i < 48; i++)
  //   mod[i] = PolarModule(i);
  //
  // std::vector<int>*tube = readTubes();
  const std::vector<int>* metai = nullptr;
  int adj_thres = 1000;

  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in getting metai from eventStore");
  }

  std::array<int, 200000> assignment;
  for (int k = 1; k < paths.size(); k++) {
    std::vector<int> path;
    int last = -1;
    for (int i : paths[k])
      if (i > 0 && metai->at(i) != last) {
      	last = metai->at(i);
      	path.push_back(i);
        }
    if (path.size() < 3) continue;

    int misses = 0;
    int s = path.size()-3;
    int ai = path[s], bi = path[s+1], ci = path[s+2];
    for (int li = metai->at(ci)+1; li < 48; li++) {
      if (next_layer[metai->at(ci)][li] < adj_thres) continue;
      int di = extend3(ctx, mod, match, ai, bi, ci, li, 0.07+!(paths.size()%2)*0.2);
      if (di > 0) {
      	path.push_back(di);
      	assignment[di] = k;
      	ai = bi;
      	bi = ci;
      	ci = di;
      	misses = 0;
      } else if (di == -1) misses++;
    }

    misses = 0;
    ai = path[0], bi = path[1], ci = path[2];
    for (int li = metai->at(ai)-1; li >= 0; li--) {
      if (next_layer[li][metai->at(ai)] < adj_thres) continue;
      int di = extend3(ctx, mod, match, ci, bi, ai, li, 0.07+!(paths.size()%2)*0.2);
      if (di > 0) {
      	path.insert(path.begin(), di);
      	assignment[di] = k;
      	ci = bi;
      	bi = ai;
      	ai = di;
      	misses = 0;
      } else if (di == -1) misses++;
      //if (misses == 2) break;
    }
    paths[k] = path;
  }
  {
    // initDensity3(ctx);

    // PolarModule mod[48];
    // for (int i = 0; i < 48; i++)
    //   mod[i] = PolarModule(i);
    //
    // vector<int>*tube = readTubes();

    double earn = 0;
    std::vector<std::pair<int, int> > pairs = findDuplicates(ctx, *mod, match);
    //cout << pairs.size() << endl;
    for (std::pair<int, int> p : pairs) {
      if (assignment[p.first] || assignment[p.second]) continue;
      int found = 0;
      {
	int ai = p.first, bi = p.second;
	std::vector<int> s[48];
	for (int li = 0; li < 48; li++) {
	  if (next_layer[li][metai->at(ai)] < adj_thres/2 &&
	      next_layer[metai->at(ai)][li] < adj_thres/2) continue;
	  if (li == metai->at(ai)) continue;

	  FW::point d, dp, xp, bap;
	  if (prepareTripleScoreDefault(ctx, ai, bi, li, d, dp, xp, bap) &&
	      prepareTripleScoreDefault(ctx, bi, ai, li, d, dp, xp, bap)) continue;

	  xp = normalize(xp)*0.99;
	  //xp = point(0,0,0);
	  const double target = 0.5;
	  double mid = findDensity(ctx, dp, xp, target, li);
	  int matches = mod->at(li).getNear(ctx, dp, xp, bap, mid, match, li);

	  double best = target;
	  std::vector<std::pair<double, int> > v;
	  for (int i = 0; i < matches; i++) {
	    int ci = match[i];
	    double s = evaluateScore(ctx, ci, dp, xp, bap);
	    v.push_back(std::make_pair(s, ci));
	  }
	  sort(v.begin(), v.end());
	  for (int i = 0; i < v.size(); i++) {
	    if (i >= 1) break;
	    found = v[i].second;
	    //s[li].push_back(v[i].second);
	    //triples.push_back(triple(ai, bi, v[i].second));
	  }
	  if (v.size()) break;
	}
      }

      std::vector<int> v;
      v.push_back(p.first);
      v.push_back(p.second);
      v.push_back(found);
      for (int i : v)
	       assignment[i] = paths.size();
      paths.push_back(v);
    }
    //cout << "Max earnings: " << earn << endl;
  }
  return paths;
}

  //Init everything needed for fast density calculations, includes most global variables above
void
FW::Topquarks::initDensity3(FW::AlgorithmContext ctx) const
{
  int Tube = 0, Disc = 1;
  int done = 1;
  const std::array<std::vector<int>, 48>* tube = nullptr;// = new std::array<std::vector<int>, 48>();
  const std::array<FW::Topquarks::Layer, 48>* layer = nullptr;
  const std::vector<FW::point>* hits = nullptr;
  const std::vector<FW::point>* polar = nullptr;
  const std::vector<int>* metai = nullptr;

  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get hits from eventStore in initDensity3");
  }
  if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get metai from eventStore in initDensity3");
  }
  if (ctx.eventStore.get("tube", tube) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get tube from eventStore in initDensity3");
  }
  if (ctx.eventStore.get("polar", polar) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get polar from eventStore in initDensity3");
  }
  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get layer from eventStore in initDensity3");
  }
  ACTS_VERBOSE("Got hits, metai, layer, tube and polar from eventStore in initDensity3");

  // const std::array<int, 200000>* assignment = nullptr;
  //
  // if (ctx.eventStore.get("assignment", assignment) == FW::ProcessCode::ABORT){
  //   ACTS_WARNING("Not able to get assignment from eventStore in initDensity3");
  // }

  // for (int i = 1; i < (int) hits->size(); i++) {
  //   // if (!assignment->at(i))
  //   tube->at(metai->at(i)).push_back(i);
  // }
  // ACTS_INFO("Created tubes for initDensity3");
  const int crude_steps = 1<<10;
  std::array<std::vector<double>, 48>* sorted_hits = new std::array<std::vector<double>, 48>();
  std::array<std::pair<double, double>, 48>* crudeIndex_a = new std::array<std::pair<double, double>, 48>();
  std::array<std::array<int, 48>, crude_steps>* crudeIndex = new std::array<std::array<int, 48>, crude_steps>();
  std::array<std::array<std::array<double, 3>, 20000>, 48>* poly = new std::array<std::array<std::array<double, 3>, 20000>, 48>();

  for (int li = 0; li < 48; li++) {
    sorted_hits->at(li).clear();
    if (layer->at(li).type == Tube) {
      for (int i : tube->at(li)){
	      sorted_hits->at(li).push_back(hits->at(i).z);
      }
    }
    else {
      for (int i : tube->at(li)){
	      sorted_hits->at(li).push_back(polar->at(i).x);
      }
    }
    ACTS_VERBOSE("Pushed " << sorted_hits->at(li).size() << " at index " << li << " of sorted_hits");
    sorted_hits->at(li).push_back(-1e50);
    sorted_hits->at(li).push_back(1e50);
    sort(sorted_hits->at(li).begin(), sorted_hits->at(li).end());

    double minx = *(std::next(sorted_hits->at(li).begin()))-1e-8;
    double maxx = *(std::next(sorted_hits->at(li).rbegin()))+1e-8;
    double f = crude_steps/(maxx-minx);
    crudeIndex_a->at(li) = std::make_pair(f, -minx*f);

    for (int i = 0; i < crude_steps; i++) {
      double x = (i+.5)/f+minx;
      crudeIndex->at(li)[i] = (int)(std::upper_bound(sorted_hits->at(li).begin(), sorted_hits->at(li).end(), x)-sorted_hits->at(li).begin());
    }
    double acc[3] = {};
    for (int i = 1; i < sorted_hits->at(li).size(); i++) {
      for (int j = 0; j < 3; j++) poly->at(li)[i][j] = acc[j];
      double x = sorted_hits->at(li)[i];
      for (int j = 0; j < 3; j++)
	      acc[j] += std::pow(x, j);
    }
  }

  if (ctx.eventStore.add("sorted_hits", std::move(*sorted_hits)) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in adding sorted_hits to eventStore in initDensity3");
  }
  if (ctx.eventStore.add("crudeIndex", std::move(*crudeIndex)) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in adding crudeIndex to eventStore in initDensity3");
  }
  if (ctx.eventStore.add("crudeIndex_a", std::move(*crudeIndex_a)) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in adding crudeIndex_a to eventStore in initDensity3");
  }
  if (ctx.eventStore.add("poly", std::move(*poly)) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Error in adding poly to eventStore in initDensity3");
  }
  ACTS_VERBOSE("Stored sorted_hits, crudeIndex, crudeIndex_a and poly in eventStore");
  // delete[]tube;

  for (double j:sorted_hits->at(18)) ACTS_VERBOSE(j << "\n");
}

void
FW::Topquarks::writeSubmission(FW::AlgorithmContext ctx, std::map<int, int>&assignment) const
{
  char filename[1000];
  sprintf(filename, "tmp/submission0.csv");
  FILE*fp = fopen(filename, "w");
  if (!fp) {
    ACTS_INFO("Could not open " << filename << " for writing submission" << "\n");
    return;
  }
  fprintf(fp, "event_id,hit_id,track_id\n");
  const std::vector<FW::point>* hits = nullptr;

  if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT){
    ACTS_WARNING("Not able to get hits frmo eventStore in writeSubmission");
  }
  for (int i = 1; i < hits->size(); i++) {
    int a = 0;
    if (assignment.count(i)) a = assignment[i];
    fprintf(fp, "%d,%d,%d\n", 0, i, a);
  }
  fclose(fp);
}
