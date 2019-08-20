#include "main.hpp"
#include "PolarModule.hpp"
#include "Point.hpp"

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/GeometryID.hpp>
#include "Acts/Utilities/Logger.hpp"
#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include <stdexcept>

#include <cmath>
#include <string>
#include <vector>

FW::PolarModuleInternal::PolarModuleInternal(int li, FW::AlgorithmContext ctx, int zi = 0) : layeri(li), zi(zi) {
    // // std::cout << "\n\n\n\n\n\n";
    layeri = li;
    int failed = 0;

    const std::vector<FW::point>* hits = nullptr;
    const std::vector<FW::point>* polar = nullptr;
    const std::vector<int>* metai = nullptr;
    const std::vector<int>* metaz = nullptr;
    const std::array<FW::Topquarks::Layer, 48>* layer = nullptr;
    const std::array<std::array<double, 4>, 48> *disc_z = nullptr;
    // ACTS_INFO("Not sure what the error is!");
    if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT);// // std::cout << "Couldn't obtain layer from eventStore in PolarModuleInternal instantiation ";

    // if (ctx.eventStore.get("assignment", assignment) == FW::ProcessCode::ABORT) failed++;

    if (ctx.eventStore.get("metai", metai) == FW::ProcessCode::ABORT);// // std::cout << "Couldn't obtain metai from eventStore in PolarModuleInternal instantiation ";

    if (ctx.eventStore.get("polar", polar) == FW::ProcessCode::ABORT);// // std::cout << "Couldn't obtain polar from eventStore in PolarModuleInternal instantiation ";

    if (ctx.eventStore.get("hits", hits) == FW::ProcessCode::ABORT);// // std::cout << "Couldn't obtain hits from eventStore in PolarModuleInternal instantiation ";

    if (ctx.eventStore.get("metaz", metaz) == FW::ProcessCode::ABORT);// // std::cout << "Couldn't obtain metaz from eventStore in PolarModuleInternal instantiation ";

    if (ctx.eventStore.get("disc_z", disc_z) == FW::ProcessCode::ABORT);// // std::cout << "Couldn't obtain disc_z from eventStore in PolarModuleInternal instantiation ";

    // std::cout << "Error in loading" << failed << " objects from eventStore\n" ;
    const FW::Topquarks::Layer &l = layer->at(layeri);
    // std::cout << l.type;
    if (l.type == Disc) {
        aspectx = 1;
        aspecty = 2*M_PI/log(l.maxr/l.minr)+.5;
    } else {
        double a = 2*M_PI*l.avgr/(l.maxz-l.minz);
        aspectx = 1;
        aspecty = 1;
        if (a > 1)
            aspecty = ceil(a);
        else
            aspectx = ceil(1./a);
    }
    // std::cout << "Gave values to aspectx and aspecty. Initializing num and ind\n";
    std::array<std::vector<int>, lods> num = std::array<std::vector<int>, lods>();
    std::array<std::vector<int>, lods> ind = std::array<std::vector<int>, lods>();
    // std::cout << "Assigning num and ind\n";
    for (int lod = 0; lod < lods; lod++) {
        ind[lod] = std::vector<int>(aspectx*aspecty<<lod*2, 0);
        num[lod] = std::vector<int>(aspectx*aspecty<<lod*2, 0);
    }
    // std::cout << "Assigned outer vectors\n";
    for (int i = 1; i < (int)hits->size(); i++) {
        // if (assignment[i]) continue;
        if (metai->at(i) == layeri && (l.type == Tube || metaz->at(i) == zi)) {
            point temp_point = polar->at(i);
            double x = calcX(ctx, temp_point, li, zi), y = polar->at(i).y*.5/M_PI+.5;
            if (!(x >= 0 && y >= 0 && x < 1 && y < 1)) {
                // // std::cout << "ERROR: hit " << i << " outside of bounds of detector " << layeri << std::endl;
                // // std::cout << x << ' ' << y << std::endl;
                continue;
            }
            for (int lod = 0; lod < lods; lod++) {
                int ix = x*(aspectx<<lod);
                int iy = y*(aspecty<<lod);
                num[lod].at(ix+iy*(aspectx<<lod))++;
            }
        }
    }
    // std::cout << "Assigning mem\n";
    int tot = 0;
    for (int j = 0; j < aspecty; j++)
        for (int i = 0; i < aspectx; i++)
            tot += num[0].at(i+aspectx*j);
    int k = 0;
    // std::cout << "Check\n";
    for (int i = 0; i < aspectx; i++)
        for (int j = 0; j < aspecty; j++)
            recIndex(ctx, k, ind, num, i, j, 0);
    // std::cout << "Performed recIndex loop";
    // std::cout << "Called recIndex\n";
    // mem = new std::array<int, tot>();
    // mem = new int[tot];
    std::map<int, int> mem;
    for (int i = 1; i < hits->size(); i++) {
      // std::cout << i << " out of " << hits->size() << "\n";
      // std::cout << metai->at(i) << "  " << layeri << "  " << l.type << "  "<< metaz->at(i) << "  " << zi << "\n";
      // if (assignment[i]) continue;
      if (metai->at(i) == layeri && (l.type == Tube || metaz->at(i) == zi)) {
          // std::cout << "1\n";
          FW::point t_point = polar->at(i);
          // std::cout << "2\n";
          double x = calcX(ctx, t_point, li, zi);
          double y = polar->at(i).y*.5/M_PI+.5;
          // std::cout << "3\n";
          int ix = x*(aspectx<<(lods-1));
          // std::cout << "4\n";
          int iy = y*(aspecty<<(lods-1));
          // std::cout << "5\n";
          mem[--ind[lods-1][ix+iy*(aspectx<<(lods-1))]] = i;
          // std::cout << "6\n";
      }
    }
    std::vector<int> data;
    data.push_back(layeri);
    data.push_back(zi);
    data.push_back(aspectx);
    data.push_back(aspecty);

    std::string str = "data_vector_";

    str += std::to_string(layeri);
    str += "_";
    str += std::to_string(zi);

    if (ctx.eventStore.add(str, std::move(data)) == FW::ProcessCode::ABORT){
      std::cout << "Error in adding data vector of PolarModules to eventStore for layer " << layeri;
    }
    std::string ind_str = "ind_";
    ind_str += std::to_string(layeri);
    ind_str += "_";
    ind_str += std::to_string(zi);

    std::string num_str = "num_";
    num_str += std::to_string(layeri);
    num_str += "_";
    num_str += std::to_string(zi);

    std::string mem_str = "mem_";
    mem_str += std::to_string(layeri);
    mem_str += "_";
    mem_str += std::to_string(zi);

    if (ctx.eventStore.add(ind_str, std::move(ind)) == FW::ProcessCode::ABORT) std::cout << "Error in storing " << ind_str << " in eventStore";
    if (ctx.eventStore.add(num_str, std::move(num)) == FW::ProcessCode::ABORT) std::cout << "Error in storing " << num_str << " in eventStore";
    if (ctx.eventStore.add(mem_str, std::move(mem)) == FW::ProcessCode::ABORT) std::cout << "Error in storing " << mem_str << " in eventStore";

    std::cout << "Successful instantiation of PolarModuleInternal\n";



    // if (ctx.eventStore.add("mem", std::move(*mem)) == FW::ProcessCode::ABORT){
    //   failed++;
    // }
    //
    // if (ctx.eventStore.add("num", std::move(*num)) == FW::ProcessCode::ABORT){
    //   failed++;
    // }
    //
    // if (ctx.eventStore.add("ind", std::move(*ind)) == FW::ProcessCode::ABORT){
    //   failed++;
    // }

    /*if (l.type == Tube)
     cout << "Tube: " << l.avgr*2*M_PI << ' ' << l.maxz-l.minz << endl;
     else
     cout << "Disc: " << l.minr << ' ' << l.maxr << ' ' << M_PI*2*l.maxr << endl;*/
}

void
FW::PolarModuleInternal::recIndex(FW::AlgorithmContext ctx, int&i, std::array<std::vector<int>, lods>&ind, std::array<std::vector<int>, lods>&num, int ix = 0, int iy = 0, int lod = 0) {

  int failed = 0;
  if (lod == lods-1) {
      i += num[lod].at(ix+iy*(aspectx<<lod));
      ind[lod].at(ix+iy*(aspectx<<lod)) = i;
  } else {
      ind[lod].at(ix+iy*(aspectx<<lod)) = i;
      for (int y = 0; y < 2; y++)
          for (int x = 0; x < 2; x++)
              recIndex(ctx, i, ind, num, ix*2+x, iy*2+y, lod+1);
  }
}

FW::PolarModuleInternal::~PolarModuleInternal() {
  std::cout<< "Destructor of PolarModuleInternal called\n";
  // if (mem) {
  //     // for (int lod = 0; lod < lods; lod++) {
  //     //     delete[]ind[lod];
  //     //     delete[]num[lod];
  //     // }
  //     delete[]mem;
  //     mem = NULL;
  // }
}

inline double
FW::PolarModuleInternal::calcX(FW::AlgorithmContext ctx, point&p, int li, int zi) const
{
  const std::array<FW::Topquarks::Layer, 48>* layer;
  const int Tube = 0, Disc = 1;

  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    // std::cout << "Not able to get layer from eventStore in double calcX";
  }
  // std::string str = "data_vector_";
  // str += std::to_string(layeri);
  // str += "_";
  // str += std::to_string(zi);
  // const std::vector<int>* temp_data_vec = nullptr;
  // if (ctx.eventStore.get(str, temp_data_vec) == FW::ProcessCode::ABORT){
  //   std::cout << "Error in getting " << str << " from eventStore in getNear";
  // }
  int layeri = li;
  // std::cout << "Got layer from eventStore in calcX";
  const FW::Topquarks::Layer&l = layer->at(layeri);
  if (l.type == Disc) //Logarithmic scaling for polar coordinates
      return p.x > l.minr ? log(p.x/l.minr)/log(l.maxr/l.minr) : 0.;//(p.x-l.minr)/(l.maxr-l.minr);
  else
      return (p.z-l.minz)/(l.maxz-l.minz);
}

inline double
FW::PolarModuleInternal::calcX(FW::AlgorithmContext ctx, double r, int li, int zi) const
{
  const std::array<FW::Topquarks::Layer, 48>* layer;
  const int Tube = 0, Disc = 1;
  // const std::vector<int>* temp_data_vec = nullptr;


  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    // std::cout << "Not able to get layer from eventStore in double calcX";
  }
  int layeri = li;
  // std::cout << "Got layer from eventStore in calcX";
  const FW::Topquarks::Layer&l = layer->at(layeri);
  if (l.type == Disc)
      return r > l.minr ? log(r/l.minr)/log(l.maxr/l.minr) : 0.;//(r-l.minr)/(l.maxr-l.minr);
  else
      return (r-l.minz)/(l.maxz-l.minz);
}

//Find all hits di with "evaluateScore(di, dp, xp, bap) < tt" in array match, return number of matches
int FW::PolarModuleInternal::getNear(FW::AlgorithmContext ctx, point&dp0, point&xp, point&bap, double tt, std::array<int, 200000> &match, int li, int zi, int pointer) const{
  const std::array<FW::Topquarks::Layer, 48>* layer;
  const int Tube = 0, Disc = 1;
  int layeri = li;
  const std::vector<int>* temp_data_vec = nullptr;
  std::string str = "data_vector_";
  str += std::to_string(layeri);
  str += "_";
  str += std::to_string(zi);
  if (ctx.eventStore.get(str, temp_data_vec) == FW::ProcessCode::ABORT){
    std::cout << "Error in getting " << str << " from eventStore in getNear";
  }
  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    // std::cout << "Not able to get layer from eventStore in double calcX";
  }
  int aspectx = temp_data_vec->at(2);
  int aspecty = temp_data_vec->at(3);

  // std::cout << "Checkin\n";
  // std::cout << layeri << zi << std::endl;
  // std::cout << "Got layer from eventStore\n";
  const FW::Topquarks::Layer&l = layer->at(layeri);
  // if (l.type == Disc) {
  //     aspectx = 1;
  //     aspecty = 2*M_PI/log(l.maxr/l.minr)+.5;
  // } else {
  //     double a = 2*M_PI*l.avgr/(l.maxz-l.minz);
  //     aspectx = 1;
  //     aspecty = 1;
  //     if (a > 1)
  //         aspecty = ceil(a);
  //     else
  //         aspectx = ceil(1./a);
  // }
  const std::vector<FW::point>* polar = nullptr;
  const std::array<std::array<double, 4>, 48>* disc_z = nullptr;
  const std::map<int, int>* temp_mem = nullptr;
  const std::array<std::vector<int>, lods>* temp_num = nullptr;
  const std::array<std::vector<int>, lods>* temp_ind = nullptr;
  std::string ind_str = "ind_";
  ind_str += std::to_string(layeri);
  ind_str += "_";
  ind_str += std::to_string(zi);

  std::string num_str = "num_";
  num_str += std::to_string(layeri);
  num_str += "_";
  num_str += std::to_string(zi);

  std::string mem_str = "mem_";
  mem_str += std::to_string(layeri);
  mem_str += "_";
  mem_str += std::to_string(zi);
  int failed = 0;

  if (ctx.eventStore.get("polar", polar) == FW::ProcessCode::ABORT) std::cout << "Error in getting polar from eventStore in getNear";
  if (ctx.eventStore.get("disc_z", disc_z) == FW::ProcessCode::ABORT) std::cout << "Error in getting disc_z from eventStore in getNear";
  if (ctx.eventStore.get(ind_str, temp_ind) == FW::ProcessCode::ABORT) std::cout << "Error in getting " << ind_str << " from eventStore in getNear";
  if (ctx.eventStore.get(num_str, temp_num) == FW::ProcessCode::ABORT) std::cout << "Error in getting " << num_str << " from eventStore in getNear";
  if (ctx.eventStore.get(mem_str, temp_mem) == FW::ProcessCode::ABORT) std::cout << "Error in getting " << mem_str << " from eventStore in getNear";

  // std::cout << "Got objects from eventStore in getNear of PolarModuleInternal\n";


  //cout << layeri << ' ' << l.minr << ' ' << l.maxr << endl;

  double yscale = 1./(dp0.x*M_PI*2);

  double ext_xm = 0, ext_xp = 0, ext_ym = 0, ext_yp = 0;
  FW::point dp = dp0;
  // const int Tube = 0, Disc = 1;
  if (l.type == Disc) {
      point off = bap*(disc_z->at(layeri)[zi]-dp0.z);
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
  if (l.type == Disc) {
      x = dp.x;
      dx = xp.x;
  } else {
      x = dp.z;
      dx = xp.z;
  }
  // std::cout << "No errors0\n";
  int lod = -log2(tt*pow(yscale*aspecty, 2))*.5;
  if (lod > lods-1) lod = lods-1;//{cout << "Too large lod: " << lod << endl; lod = lods-1;}
  if (lod < 0) lod = 0;
  int stride = aspectx<<lod;
  // std::cout << "You found it\n";
  double inv = 1./(1-dx*dx-dy*dy);
  double rx = sqrt(tt*(1-dy*dy)*inv);
  double ry = yscale*sqrt(tt*(1-dx*dx)*inv);
  if (rx != rx || ry != ry || x != x || y != y) {
      /*            if (DEBUG) {
       // std::cout << "NaN in PolarModule" << std::endl;
       // std::cout << layeri << ' ' << dp.x << ' ' << x << ' ' << y << ' ' << rx << ' ' << ry << std::endl;
       exit(0);
       } else
       */
      return 0;
  }
  // std::cout << "No errors till now";

  int sx =  calcX(ctx, x-rx-ext_xm, li, zi)*(aspectx<<lod);
  // std::cout << "First call works\n";
  int sy = floor((y-ry-ext_ym)*(aspecty<<lod));
  // std::cout << "Floored"
  int ex =  calcX(ctx, x+rx+ext_xp, li, zi)*(aspectx<<lod)+1;
  // std::cout << "Second call works\n";
  int ey = floor((y+ry+ext_yp)*(aspecty<<lod))+1;
  // std::cout << "No errors till now2";
  //cout << sx << ' ' << ex << ' ' << sy << ' ' << ey << ' ' << (1<<lod) << endl;
  sx = std::max(sx, 0);
  ex = std::min(ex, stride);

  ey = std::min(ey, sy+(aspecty<<lod));
  // std::cout << "No errors till now3";
  //cout << (ex-sx)*(ey-sy) << endl;
  int considered = 0;
  int matches = 0;

  int ii = sy;
  int i = ii%(aspecty<<lod);
  if (i < 0) i += aspecty<<lod;
  // std::cout << "No errors till now4";
  // std::cout << "Loop length " << ex - sx << "\n";
  do {
    for (int j = sx; j < ex; j++) {
      // std::cout << "Running " << j << " out of " << ex << "\n";
      // const int&n = num[lod].at(i*stride+j);
      const int&n = temp_num->at(lod).at(i*stride+j);
      if (!n) continue;
      //if (rx*rx+ry*ry-dot*dot < t2) //TODO
      considered += n;
      // int k0 = ind[lod].at(i*stride+j);
      int k0 = temp_ind->at(lod).at(i*stride+j);
      for (int k = 0; k < n; k++) {
        // std::cout << "Running inner " << k << " out of " << n << "\n";
        // int hit_id = mem->at(k0 + k);
        // int hit_id = mem[k0+k];
        int hit_id = temp_mem->at(k0 + k);
        if ((hit_id == 40350) || (hit_id==73133) || (hit_id==73138)) std::cout << "In a correct hit_id for hardcoded ai and bi";
        // std::cout << "Yes that is the error";
        point r = polar->at(hit_id);
        //if (r.x < l.minr || r.x > l.maxr) cout << "What?" << endl;
        point err = r-dp0;
        if (err.y > M_PI) err.y -= M_PI*2;
        if (err.y <-M_PI) err.y += M_PI*2;
        err.y *= dp0.x; //Why doesn't dp.x work better than r.x?

        err = err-bap*(l.type == Disc ? err.z : err.x);
        double dot = err*xp;
        double r2 = err*err-dot*dot;
        if (r2 < tt) match.at(pointer + matches++) = hit_id;
      }
    }
    ii++;
    if (++i == (aspecty<<lod)) i = 0;
  } while (ii < ey);
  //cout << considered << endl;
  //cout << endl;
  // std::cout << "Found " << matches << " matches";
  return matches;
}


FW::PolarModule::PolarModule(int li, FW::AlgorithmContext ctx) {
  const std::array<FW::Topquarks::Layer, 48>* layer = nullptr;

  if (ctx.eventStore.get("layer", layer) == FW::ProcessCode::ABORT){
    std::cout << "Error in getting layer from eventStore in PolarModule instantiation";
  }

  std::cout << "Got layer from eventStore in PolarModule instantiation";
  layeri = li;
  const FW::Topquarks::Layer&l = layer->at(li);
  if (l.type == PolarModuleInternal::Disc) {
      internals = 4;
      // std::array<PolarModuleInternal, 4>* internal = {new PolarModuleInternal{{li, ctx, 0}, {li, ctx, 1}, {li, ctx, 2}, {li, ctx, 3}}};
      std::vector<PolarModuleInternal> internal_data = {PolarModuleInternal(li, ctx, 0), PolarModuleInternal(li, ctx, 1), PolarModuleInternal(li, ctx, 2), PolarModuleInternal(li, ctx, 3)};
      std::vector<PolarModuleInternal>* internal = &internal_data;
  } else {
      internals = 1;
      // std::array<PolarModuleInternal, 1>* internal = PolarModuleInternal[1]{{li,ctx,0}};
      std::vector<PolarModuleInternal> internal_data = {PolarModuleInternal(li, ctx, 0)};
      std::vector<PolarModuleInternal>* internal = &internal_data;
      // std::array<PolarModuleInternal, 1>* internal = {PolarModuleInternal(li, ctx, 0)};
  }
  // if (ctx.eventStore.add("internal0", std::move(internal_data)) == FW::ProcessCode::ABORT){
  //   std::cout << "Error in moving internal to eventStore in PolarModule instantiation";
  // }
}

int
FW::PolarModule::getNear(FW::AlgorithmContext ctx, point&dp, point&xp, point&bap, double tt, std::array<int, 200000> &match, int li) const{
    int matches = 0;
    int pointer = 0;
    for (int i = 0; i < internals; i++) {
      // std::cout << "Running loop " <<  i << " out of " << internals << " in wrapper getNear\n";
        pointer += matches;
        matches += internal->at(i).getNear(ctx, dp, xp, bap, tt, match, li, i, pointer);
    }
    return matches;
}

FW::PolarModule::~PolarModule(){
  std::cout << "Destructor of PolarModule called\n";
}
