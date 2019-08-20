#ifndef Engine_H
#define Engine_H

#include "util.h"
#include <iostream>
#include "TrackModelPhysical.h"
#include "SearchLayer.h"
#include "DataStructures.h"

struct Tracker;
struct SearchHit;

#include "Cuts.h"
#include <list>

using namespace std;



struct TrackCandidate
{
  int nLayers;
  int nHits;
  double chi2;
  TrackModelPhysical track3;
  FitHit hit3;
  int hitIDs[gMaxTracktHits];
};


struct CombinatorialKnot
{
  CombinatorialKnot() = default;

  int iLayer=0;
  // parameters of the tracklet, extrapolated to the layer
  double extrPhi=0;
  double extrT=0;
  int localGap=0;  
  int totalGap=0;  
  int nMaterialsCrossed=0;  
  // hits on the layer
  SearchLayerArea hitArea;
  // hit(s) for the current combinatorial branch
  FitHit hit;
  //int iHit1=0;
  int hitID;
  int nClones=0;
  TrackModelPhysical updatedTrack; // updated with current hit
  double totalChi2=0.; // updated with current hit
};


class Engine
{
 public:

  Engine() = default;
  ~Engine() = default;


  enum EndPointFlag {Off, On };

  std::vector<Track> mTracks;
  Tracker *mTracker=0;  
  const RecoPassParameters *mRecoPassParameters=0;

  CombinatorialKnot mTree[Geo::NLayers*2];  
  TrackCandidate mCurrentTrack;
  TrackCandidate mBestTrack;

  int mLayerStart=0;
  int mLayerEnd=0 ;
  int mLayerIncr=1 ;
  EndPointFlag  mIsEndpointSet=EndPointFlag::Off;
  FitHit mEndPoint;
  int mNCombLayers=0;
  int mMinNLayers=0;  
  int mLastExtrLayer = 0;
  int mStatFitErrors=0;  
  int mNBranches=0;

  bool mDoPrint=0;


  void operator()( int iThread, int nThreads,		 
		 int zSide);


  bool findNextLayer( bool &isLayerInnerCrossed );
  bool findNextHit();
  void pickUpDuplicates( int iCombLayer );
 
  void Exec( int iThread, int nThreads,
	     int zSide );

  void prolongateTrack( int layerStart, int layerEnd, EndPointFlag isEndPointSet, int minNLayers );

};




#endif
