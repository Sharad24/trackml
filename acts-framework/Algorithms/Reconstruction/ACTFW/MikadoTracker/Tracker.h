#ifndef Tracker_H
#define Tracker_H

#include "util.h"
#include <iostream>
#include "TrackModelPhysical.h"
#include "DataStructures.h"
#include "Geo.h"
#include "SearchLayer.h"
#include "AccuracyEvaluator.h"
#include "Learning.h"
#include "Engine.h"
#include "TrackSelector.h"

class Tracker
{
 public:
  Tracker() = default;

  ~Tracker(){
    delete[] mRecThreads;
  }
  /*
  int readEvent( float *x, float *y, float *z, int *id, int *vol, int *layer, int nHits );
  int readEvent( const char *directory, int event, bool loadMC );
  int readMCEvent( const char *directory, int event );
  */

  void reconstruct( int LearningMode );

  void doReconstructionPass( int iPass, bool force = 0 );

  void doReconstructionPassForce( int iPass ){
    doReconstructionPass( iPass, 1 );
  }


  void analyzeGeometry( bool endOfData );
  void analyzeField( bool endOfData );
  bool checkLayers( const int *baseLayers, int nBaseLayers, Particle &p, double &fitpt, double &phi );

  void TrackFitTest( const RecoPassParameters &cuts );
  void TrackFitTest( int ilayer1, int ilayer2, int ilayer3 );

  void Print( const Track &t ); 


  std::vector<Hit> mHits;

  std::vector<HitMC> mHitsMC;
  std::vector<Particle> mParticles;
  
  std::vector<SearchHit> mLayerHits[Geo::NLayers];

  bool mcFlag = 0;
  int mNThreads = 1;
  AccuracyEvaluator mQA = AccuracyEvaluator(this);

  const RecoPassParameters *mRecoPassParameters=0;

  int mMaxNBranches=-1;
  int mFirstNonStopVolume = 6;
  
  bool doPrint=1;  
  bool mUseMC=0;
  double wrongHitMix = 0.;


  std::vector<SearchLayer> mSearchLayers = std::vector<SearchLayer>(Geo::NLayers);
  std::vector<SearchLayer> mPickUpLayers = std::vector<SearchLayer>(Geo::NLayers);

  std::vector<Track> mTrackCandidates;
  std::vector<Track> mTracks;


  bool isInMix( const HitMC &h ) {
    return ( h.random < wrongHitMix );
  }

  Engine *mRecThreads=0;
  int mNRecThreads = 0;  
  TrackSelector selector = TrackSelector( this );
  int mRecoPassForLearning=-1;
  bool  mIsLearningPointReached = 0;
};



#endif
