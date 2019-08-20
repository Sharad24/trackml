#include "TrackSelector.h"
#include "Tracker.h"
#include "util.h"
#include <iostream>
#include "Cuts.h"

using namespace std;


struct SortObject
{
  int iTrack;
  int nActiveHits;
  double chi2;   
  bool operator<( const SortObject o ) const {
    if( nActiveHits < o.nActiveHits ) return 1;
    if( nActiveHits > o.nActiveHits ) return 0;       
    if( chi2 > o.chi2 ) return 1;
    if( chi2 < o.chi2 ) return 0;
    return ( iTrack > o.iTrack );
  }
};


struct SortHit
{
  int isUsed; 
  int layer;  
  SortObject s;
};

struct SortTrack
{
  bool skip;
  int nHits;
  SortObject s;
};


void TrackSelector::SelectTracks( )
{  
  // if( doPrint) cout<<"Sorting "<<mTrackCandidates.size()<<" tracks... "<<endl;

  int minTrackNLayers =  mTracker->mRecoPassParameters->minTrackNLayers;
  if( mTracker->mTrackCandidates.size()==0 ) return;
  
  int nHits = (int) mTracker->mHits.size();
  int nTracks = (int)  mTracker->mTrackCandidates.size();

  SortHit *sHits = new SortHit[nHits];
  SortTrack *sTracks = new SortTrack[nTracks];
 
  // init hits

#ifdef SGSafeRun    
  for( int i=0; i<nHits; i++ ){
    Hit &h = mTracker->mHits[i];
    SortHit &sh = sHits[i];
    sh.isUsed = h.isUsed;    
    sh.layer = h.layerID;
    sh.s.iTrack = -1;
    sh.s.nActiveHits = -1;
    sh.s.chi2 = 1.e10;
  }
#endif

  // init tracks

  for( int itr=0; itr<nTracks; itr++ ){
    Track &t =  mTracker->mTrackCandidates[itr];
    SortTrack &st = sTracks[itr];    
    st.nHits = t.nHits;    
    st.s.iTrack = itr;    
    st.s.chi2 = t.chi2/t.nHits;    
    st.s.nActiveHits = 0;
    for( int i=0; i<t.nHits; i++ ){
      int id = t.hitIDs[i];
      Hit &h = mTracker->mHits[id];
      SortHit &sh = sHits[id];
      sh.isUsed = h.isUsed;    
      sh.layer = h.layerID;
      sh.s.iTrack = -1;
      sh.s.nActiveHits = -1;
      sh.s.chi2 = 1.e10;
      if( !sh.isUsed ) st.s.nActiveHits++;
    }
    st.skip = 0;
    if( st.s.nActiveHits < 3 || st.s.nActiveHits < (int) t.nHits - 2 ) st.skip = 1;
    st.skip = st.skip || ( st.s.nActiveHits < minTrackNLayers );   
#ifdef SGSafeRun    
    if(1){ // check
      for( int i=0; i<t.nHits; i++ ) {      
	for( int j=i+1; j<t.nHits; j++ ) {      
	  if( t.hitIDs[i] == t.hitIDs[j] ){
	    cout<<"Sorting: Hit is assigned twice!!"<<endl;	  
	    for( int k=0; k<t.nHits; k++ ) {
	      mTracker->mHits[t.hitIDs[k]].Print();
	    }
	  }
	}      
      }
    }
#endif
  }  

  // competition between tracks

  for( int itr=0; itr<nTracks; itr++ ){    
    Track &t =  mTracker->mTrackCandidates[itr];
    SortTrack &st = sTracks[itr];
    if( st.skip ) continue;    
    for( int i=0; i<t.nHits; i++ ){
      SortHit &sh = sHits[t.hitIDs[i]];
      if( sh.isUsed ) continue;
      if( st.s < sh.s ) continue;      
      sh.s = st.s; // mark the hit
    }
  }
  
  // kill some tracks

  constexpr int nLayers = Geo::NLayers;
  int nTrHitsPerLayer[nLayers];

  for( int itr=0; itr<nTracks; itr++ ){    
    Track &t =  mTracker->mTrackCandidates[itr];
    SortTrack &st = sTracks[itr];
    if( st.skip ) continue;        
   
    // number of layers 
    for( int il=0; il<nLayers; il++ ){
      nTrHitsPerLayer[il]=0;
    }
    
    int nActive=0;
    for( int i=0; i<t.nHits; i++ ){
      SortHit &sh = sHits[t.hitIDs[i]];
      if( sh.s.iTrack!=itr ) continue;
      nActive++;
      nTrHitsPerLayer[sh.layer]++;
    }

    int tnLayers=0;
    for( int il=0; il<nLayers; il++ ){
      if( nTrHitsPerLayer[il]>0 ) tnLayers++;
    }   
       
    if( nActive < 3 || nActive < st.nHits - 2 || tnLayers < minTrackNLayers ){ // kill 
      st.skip=1; 
      continue;
    }
    
    // mark the hits as used
    for( int i=0; i<t.nHits; i++ ){
      SortHit &sh = sHits[t.hitIDs[i]];
      if( sh.s.iTrack ==itr ) sh.isUsed = 1;
  
    }
  }

  // store tracks

  for( int itr=0; itr<nTracks; itr++ ){
    Track &t =  mTracker->mTrackCandidates[itr];
    SortTrack &st = sTracks[itr];
    if( st.skip )     continue;
    
    Track tt = t;
    tt.nHits=0;    
    st.skip = 1;
    for( int ih=0; ih<t.nHits; ih++ ){
      SortHit &sh = sHits[t.hitIDs[ih]];
      Hit &h = mTracker->mHits[t.hitIDs[ih]];
      if( h.isUsed ) continue;
      if( sh.s.iTrack == itr || !sh.isUsed ){
	sh.isUsed = 1;
	h.isUsed = 1;
	tt.hitIDs[tt.nHits++] = t.hitIDs[ih];
      }
    }
     mTracker->mTracks.push_back(tt);
  }
  mTracker->mTrackCandidates.clear();

  delete[] sHits;
  delete[] sTracks;

  //if( doPrint) cout<<"...sorting"<<endl;
}




