#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#include <iostream>
#include <vector>


using namespace std;

struct FitHit{
#ifdef SGSafeRun
FitHit() : x(0.f), y(0.f), z(0.f), BzFwd(0.f), BzMid(0.f), BzBck(0.f)  {}
#endif
  float x;
  float y;
  float z;
  float BzFwd;
  float BzMid;
  float BzBck;  
};

 
struct Hit // container for hit information
{
  FitHit fitHit;
  float x() const { return fitHit.x; }
  float y() const { return fitHit.y; }
  float z() const { return fitHit.z; }
  float BzFwd() const { return fitHit.BzFwd;  }
  float BzMid() const { return fitHit.BzMid;  }
  float BzBck() const { return fitHit.BzBck;  }  
  float r;
  float phi;
  float t;
  int layerID;
  int isUsed;
  int trackID;
  void Print() const;
};


constexpr int gMaxTracktHits = 24;

struct Track
{
  int nHits;
  int nLayers; 
  double chi2;
  int hitIDs[gMaxTracktHits];  
};

struct HitMC // hit truth data, stored in an array parallel to hits
{  
  int hitID; // hit id
  double x;
  double y;
  double z;
  int partID;
  double w; // weight 
  double px;
  double py;
  double pz;
  double pt;
  double dOrigin2; // squared distance to particle origin 
  double p; // momentum for sorting
  double q;
  double r;
  double phi;
  double t;
  double random;
  void Print() const;

  bool operator< ( const HitMC &h ) const
  {
    return ( p > h.p ) || (p==h.p && dOrigin2 < h.dOrigin2);
  }
};


struct Particle // structure for truth particle info
{
  Particle(int nhits=0)
    : hits(nhits)
  {
    hits.clear();
  }
  ~Particle() = default;
  void Print() const;

  long unsigned int origID;
  int countID;
  int pid;
  int isSelected;
  double x;
  double y;
  double z;
  double r;  
  double px;
  double py;
  double pz;
  double pt;
  double p;
  double q;
  double w;
  double dcaR;
  double dcaZ;
  double xl;
  double yl;
  double zl;
  double rl;

  int nLayers;
  std::vector<int> hits;
  std::vector<int> hitClusterIds; // along the trajectory
  int prim; // is it coming from the origin
  bool baseV;  // some combinations of layers it crosses
  bool baseV0; 
  bool baseV1; 
  bool baseV2; 
  bool baseA0;
  bool baseA1;
  bool baseA2;
  bool baseA3;
  double ptV;
  double ptV0;
  double ptV1;
  double ptV2;
  double ptA0;
  double ptA1;
  double ptA2;
  double ptA3;
  double phiV;
  double phiV0;
  double phiV0L12;
  double phiV0L23;
  double phiV0V1;
  double phiV1;
  double phiV2;
  double phiV1L23;
  double phiV2L23;
  double phiV1L45;
  double phiV2L45;
  double phiV1L56;
  double phiV2L56;
  double phiV0V3;
  double phiV3L01;
 double phiV3L12;

  double phiA0;
  double phiA1;
  double phiA2;
  double phiA3;
};


#endif
