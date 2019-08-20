#ifndef Learning_H
#define Learning_H

#include <iostream>
#include "DataStructures.h"
#include "Geo.h"
#include "Tracker.h"
#include "Cuts.h"

struct Tracker;

using namespace std;

class Learning
{
 public:

  struct Event{
    vector<Hit> mHits;
    vector<HitMC> mHitsMC;
    vector<Particle> mParticles;
    vector<Track> mTrackCandidates;
    vector<Track> mTracks;
    vector<SearchHit> mLayerHits[Geo::NLayers];
  };

  Learning();  
  ~Learning() = default;

  void init();
  void StoreEvent( const Tracker &tracker );
  void RestoreEvent( int iev, Tracker &tracker );

  void SetCalibrationObject( RecoPassParameters &cuts, bool sskipDead );

  void selectAll( bool v );
  void selectParameters( bool v );
  void selectParameter( int iPar, bool v );
  void selectVolume( int iVolume, bool v );
  void selectSearchLayer( int il, bool v );
  void selectPickUpLayer( int il, bool v );

  void selectSearchLayers( bool v );
  void selectPickUpLayers( bool v );
  void selectPhiT( bool v );
  void selectUV( bool v );
  void selectUVms( bool v );

  void selectDead( RecoPassParameters &p, double multSearch, double multFit, double multPick );

  void extendSearchLayers( RecoPassParameters &cuts, double v );
  void extendPickUpLayers( RecoPassParameters &cuts, double v );

  void setDefaultCutsOpen( const RecoPassParameters &par, double multSearch, double multFit, double multPick );
  void setDefaultCutsClosed();
  void setDefaultCuts( const RecoPassParameters &par, double multParam, double multSearch, double multFit, double multPick );
  void setDefaultCuts( const RecoPassParameters &par, double mult ){
    setDefaultCuts( par, mult, mult, mult, mult );
  } 
  
  void extendToMin( RecoPassParameters &p ) const;
  void resetPhi( RecoPassParameters &p, double mult ) const ;
  void resetT( RecoPassParameters &p, double mult ) const ;
  void resetUV( RecoPassParameters &p, double mult ) const;


  void writeDefaultCuts( const char *file );
  void readDefaultCuts( const char *file );


  void StartLearning();

  bool startNewTry(); 
  void setTryResult( double eff, double info[], int nnInfo );

  void PrintStatus( double effTry );

  vector<Event> mEvents;  

  RecoPassParameters *mCuts=0;
  RecoPassParameters &mCutsOrig = * new RecoPassParameters;
  RecoPassParameters &mCutsBest = * new RecoPassParameters;

  double infoOrig[10];
  double infoBest[10];
  double infoSample[10][10];
  int nInfo=0;

  double mEffBest = -1.;
  double mEffOrig = -1;
  double mEffCurrent = -1;

  int mCurrentParameter = -1;

  double *mCurrentParameterP = nullptr;
  double *mBestParameterP    = nullptr;

  int mSample;
  int mSampleElement = 0;
  double mSampleValues[4];
  double mSampleEff[4];
  bool mSkipDead;  
  bool parameterSelection[ RecoPassParameters::nPar() ];

  RecoPassParameters mDefaultCuts;

  int startSample=0;

  vector<int> vSelectedParameterIDs;

 private:
  
  void setupNewSample( int iPar, int iSample, bool tryStrickCut );

};


#endif
