
/*
  TNtuple *nt = (TNtuple*) gDirectory->FindObjectAny("edgeHits");
  "part:nhits:pt:p:w:side:x:y:z:r"
  TNtuple *nt = (TNtuple*) gDirectory->FindObjectAny("fitHits");
"part:nhits:pt:p:w:fittype:type:x:y:z:r:du:dv:dz:dr:fitpt");
TH1F *h = (TH1F*) gDirectory->FindObjectAny("weights");
 */

#include "Tracker.h"
#include "util.h"
#include "TrackModelPhysical.h"
#include "Engine.h"
#include "AccuracyEvaluator.h"

#include "Cuts.h"
#include <thread>
#include <sys/time.h>
using namespace std;



void Tracker::reconstruct( int LearningMode )
{
  if( doPrint ) cout<<"Reconstruction.."<<endl;

  clock_t startClock = clock();
  struct timeval startTime;
  gettimeofday( &startTime, NULL );

  if( LearningMode ==2 ) goto LearningEntryPoint;

  //ReconstructMC();
  mTrackCandidates.clear();
  mTracks.clear();

  mMaxNBranches=-1;

  mFirstNonStopVolume = 6;

  mIsLearningPointReached = 0;

  doReconstructionPass( 0 );
  doReconstructionPass( 1 );
  doReconstructionPass( 5 );
  doReconstructionPass( 7 );

  doReconstructionPass( 2 );
  doReconstructionPass( 3 );
  doReconstructionPass( 6 );
  doReconstructionPass( 8 );

  doReconstructionPass( 4 );
  doReconstructionPass( 9 );
  doReconstructionPass( 10 );

  // V0 Pass 2: gap 1, stoped
  doReconstructionPass( 11 );
  doReconstructionPass( 61 );
  doReconstructionPass( 12 );
  doReconstructionPass( 13 );
  doReconstructionPass( 14 );
  doReconstructionPass( 15 );
  doReconstructionPass( 16 );
  doReconstructionPass( 17 );

  // V0 Pass 3
  doReconstructionPass( 18 );
  doReconstructionPass( 19 );

  // V1 Pass 2

  doReconstructionPass( 20 );
  doReconstructionPass( 21 );

  // V1 Pass 3
  doReconstructionPass( 25 );
  doReconstructionPass( 27 );
  doReconstructionPass( 26 );
  doReconstructionPass( 22 );

  // V1 L2&3
  doReconstructionPass( 30 );
  doReconstructionPass( 62 );

  // V1 L4&5
  doReconstructionPass( 31 );

  // V1 Pass 4
  doReconstructionPass( 29 );

  // V2 Pass 2
  doReconstructionPass( 32 );
  doReconstructionPass( 33 );

  // V2 Pass 3
  doReconstructionPass( 34 );
  doReconstructionPass( 35 );
  doReconstructionPass( 36 );
  doReconstructionPass( 37 );

  // V2 L2&3
  doReconstructionPass( 38 );
  doReconstructionPass( 39 );

  // V2 L4&5
  doReconstructionPass( 40 );

  // V2 Pass 4
  doReconstructionPass( 41 );

  // V0 Pass 3 continuation..
  doReconstructionPass( 42 );
  doReconstructionPass( 43 );

  // V0 L1+L2
  doReconstructionPass( 44 );
  doReconstructionPass( 45 );
  doReconstructionPass( 46 );

  // V0 Pass 4: non-stop before 3 (does not really needed )

  doReconstructionPass( 47 );
  doReconstructionPass( 48 );
  doReconstructionPass( 49 );

  // V0 Pass 5: non-stop before 3 (does not really needed )
  doReconstructionPass( 50 );

  if(1){ // not tuned

    // V0 next pass
    // cut n hits is wrong (has been tuned to 2)

    doReconstructionPass( 51 );
    doReconstructionPass( 52 );

    // V0 last pass
    doReconstructionPass( 53 );
    doReconstructionPass( 54 );

    // V0 L12 last pass
    doReconstructionPass( 55 );

    // V0 L23 last pass
    doReconstructionPass( 56 );

    // V0 L4 V1 L0 last pass
    doReconstructionPass( 57 );

    // V3 L01 last pass
    doReconstructionPass( 58 );

    // V3 L12 last pass
    doReconstructionPass( 59 );

  }


  if( mRecoPassForLearning>=0 ){
    mTracks.clear();
  }

 LearningEntryPoint: if( LearningMode==1 ) return;

  if( mRecoPassForLearning >= 0 ){
    mUseMC = 0;
    wrongHitMix = 0.;
    mTracks.clear();
    doReconstructionPassForce( mRecoPassForLearning );
    const RecoPassParameters &cuts = Cuts::sliceCuts[mRecoPassForLearning];
    mQA.SelectParticles( cuts, cuts.l2PhiMin, cuts.l2Phi );
  }

  clock_t stopClock = clock();
  struct timeval stopTime;
  gettimeofday( &stopTime, NULL );

  double runClock = ((double) (stopClock-startClock))/CLOCKS_PER_SEC;
  struct timeval diffTime;
  timersub( &stopTime, &startTime, &diffTime );

  double runTime = diffTime.tv_sec * 1. + diffTime.tv_usec * 1.e-6;

  static double totalClock=0;
  static double totalTime=0;
  static double totalCalls=0;

  totalTime+=runTime;
  totalClock+=runClock;
  totalCalls++;

  if( doPrint ){
    static int nEvents = -1;
    nEvents++;
    //std::cout << " number of cores: " << std::thread::hardware_concurrency() << endl;
    cout<<"\nEvent "<<nEvents<<", reconstructed "<<mTracks.size()<< " tracks, processing time: "<<runTime<<" total: "<<totalTime/totalCalls <<" cpu: "<<totalClock/totalCalls<<endl;
  }
  //doPrint = 0;

  mQA.Evaluate( doPrint );

  //doPrint=1;

  //analyzeGeometry();
}

inline void createLayers( Tracker &tracker, int iThread, int nThreads )
{
  const RecoPassParameters &par = *tracker.mRecoPassParameters;

  for( int id = iThread; id<Geo::NLayers; id+= nThreads){

    SearchLayer &layer = tracker.mSearchLayers[ id ];
    SearchLayer &pickUpLayer = tracker.mPickUpLayers[ id ];
    layer.mIsOn = 0;
    pickUpLayer.mIsOn = 0;

    const LayerCuts &cuts = par.searchCuts[id];
    double phi = cuts.cutPhi;
    double t = cuts.cutT;
    if( id==par.baseLayer1 || id==par.baseLayer2 ){ // base layers
      phi = par.l2Phi;
      t = par.l2T;
    } else if( phi < 1.e-4 ) { // phi == 0 => layer not in use for the current reconstruction pass
#ifndef SgSafeRun
      continue;
#endif
    }
    layer.mIsOn = 1;
    bool makeClones = (id!=par.baseLayer1);

    SearchLayer::CleanUp( &tracker, id );

    layer.Create( &tracker, id, 2*phi, 2*t, makeClones );
    // corresponding pick up layer
    const LayerCuts &pickUpCuts = par.pickUpCuts[id];
    if( pickUpCuts.cutPhi <  1.e-6 ) {
#ifndef SgSafeRun
      continue;
#endif
    }
    pickUpLayer.mIsOn = 1;
    //pickUpLayer.CreatePickUpLayer( layer, 2*pickUpCuts.cutPhi, 2*pickUpCuts.cutT );
    pickUpLayer.Create( &tracker, id, 2*pickUpCuts.cutPhi, 2*pickUpCuts.cutT, 1 );
  }
}



void Tracker::doReconstructionPass( int iPass, bool force )
{
  // if( doPrint ) cout<<"----------------- track constructor .."<<endl;

  if( !force ){
    if( iPass == mRecoPassForLearning ){
      mIsLearningPointReached=1;
    }
    if( mIsLearningPointReached ) return;
  }

  const RecoPassParameters &par = Cuts::sliceCuts[ iPass ];

  //struct timeval startTime;
  //gettimeofday( &startTime, NULL );

  mRecoPassParameters = &par;
  mFirstNonStopVolume = mRecoPassParameters->firstNonStopVolume;

  int nThreads = mNThreads;
  //nThreads = 2;

  // two selection layers

  SearchLayer &baseLayer1 = mSearchLayers[par.baseLayer1];
  SearchLayer &baseLayer2 = mSearchLayers[par.baseLayer2];

  int zSide = 0;
  if( baseLayer1.mType==0 && baseLayer2.mType!=0 ){
    zSide = baseLayer2.mType;
  }


  // combinatorial layers

  int nThreads0 = nThreads;
  //nThreads0 = 1;
  std::thread *tThreads0[nThreads0];
  for( int i=0; i<nThreads0; i++ ){
    tThreads0[i] = new std::thread( createLayers, std::ref(*this), i, nThreads0 );
    if( !tThreads0[i] ){
      std::cout<<"Can not create a thread!!! "<<std::endl;
      exit(-1);
    }
  }
  for( int i=0; i<nThreads0; i++ ){
    tThreads0[i]->join();
    delete tThreads0[i];
    tThreads0[i] = nullptr;
  }

  //nThreads = 1;

  if( mNRecThreads < nThreads ){
    mNRecThreads = nThreads;
    delete[] mRecThreads;
    mRecThreads = new Engine[mNRecThreads];// = new Engine[nThreads];
    if( !mRecThreads ){
      std::cout<<"Can not create thread objects!!! "<<std::endl;
      exit(-1);
    }
    uint capacity = 10000;
    if( capacity<mTrackCandidates.capacity() ) capacity = mTrackCandidates.capacity();
    for( int i=0; i<nThreads; i++ ){
      mRecThreads[i].mTracker = this;
      mRecThreads[i].mRecoPassParameters = mRecoPassParameters;
      mRecThreads[i].mDoPrint = ( doPrint && (nThreads==1) );
      mRecThreads[i].mTracks.reserve( capacity );
    }
  }

  std::thread *tThreads[nThreads];

  for( int i=0; i<nThreads; i++ ){
    mRecThreads[i].mTracker = this;
    mRecThreads[i].mRecoPassParameters = mRecoPassParameters;
    mRecThreads[i].mDoPrint = ( doPrint && (nThreads==1) );

    tThreads[i] = new std::thread( std::ref(mRecThreads[i]), i, nThreads, zSide  );

    if( !tThreads[i] ){
      std::cout<<"Can not create a thread!!! "<<std::endl;
      exit(-1);
    }
  }

  int ntr=0;
  for( int i=0; i<nThreads; i++ ){
    tThreads[i]->join();
    ntr+=mRecThreads[i].mTracks.size();
    delete tThreads[i];
  }
  mTrackCandidates.reserve(mTrackCandidates.size()+ntr);
  for( int i=0; i<nThreads; i++ ){
    mTrackCandidates.insert( mTrackCandidates.end(), mRecThreads[i].mTracks.begin(), mRecThreads[i].mTracks.end() );
  }

  selector.SelectTracks();
  mTrackCandidates.clear();

  //struct timeval stopTime;
  //gettimeofday( &stopTime, NULL );

  //struct timeval diffTime;
  //timersub( &stopTime, &startTime, &diffTime );
  //double runTime = diffTime.tv_sec * 1. + diffTime.tv_usec * 1.e-6;

  if( doPrint ){
    // cout<<"slice reconstruction time "<< runTime<<" sec "<<endl;
  }
}


void Tracker::Print( const Track &t )
{
  int nMC[Geo::NLayers];
  int nRec[Geo::NLayers];
  int nRecF[Geo::NLayers];

  HitMC &mc0 = mHitsMC[ t.hitIDs[0] ];
  if( mc0.partID<0 ) return;
  Particle &p = mParticles[mc0.partID];

  for( int i=0; i<Geo::NLayers; i++ ){
    nMC[i] = 0;
    nRec[i] = 0;
    nRecF[i] = 0;
  }

  int nFakeHits=0;
  int nTrueHits=0;
  for( int i=0; i<t.nHits; i++){
    Hit &h = mHits[ t.hitIDs[i] ];
    HitMC &mc = mHitsMC[ t.hitIDs[i] ];
    if( mc.partID == mc0.partID ){
      nRec[h.layerID]++;
      nTrueHits++;
    } else {
      nRecF[h.layerID]++;
      nFakeHits++;
    }
  }

  for( uint i=0; i<p.hits.size(); i++){
    Hit &h = mHits[ p.hits[i] ];
    //HitMC &mc = mHitsMC[ t.vHits[i] ];
    //mc.Print();
    nMC[h.layerID]++;
  }

  cout<<"Track: part "<<mc0.partID<<" nhits mc "<<p.hits.size()<<" rec "<<nTrueHits<<" + fake "<< nFakeHits<<" : "<<endl;
  for( int i=0; i<Geo::NLayers; i++ ){
    if( !nMC[i] && !nRec[i] && !nRecF[i]) continue;
    int iv = Geo::getVolume(i);
    int il = Geo::getVolumeLayer(i);
    cout<<" v "<<iv<<" l "<<il<<": mc "<<nMC[i]<<" rec "<<nRec[i]<<" fake "<<nRecF[i]<<endl;
  }
}
