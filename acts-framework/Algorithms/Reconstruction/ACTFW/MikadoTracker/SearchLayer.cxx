
#include "SearchLayer.h"

#include "util.h"
#include <iostream>
#include "TrackModelPhysical.h"
#include "Tracker.h"

using namespace std;

constexpr double PI = 3.14159265359;

void  SearchHit::Print() const
{
  //cout<<"fit hit x "<<x<<" y "<<y<<" z "<<z<<" r "<<sqrt(x*x+y*y)<<" phi "<<phi<<" t "<<t<<" search phi "<<phi<<endl;  
}

void SearchLayer::Print() const
{
  cout<<"SearchLayer v "<<mVolume<<" l "<< mLayer<<" r "<<mR<<" z "<<mZ<<endl;
}


void SearchLayer::CleanUp( Tracker *tracker, int layerID )
{
  std::vector<SearchHit> &layerHits = tracker->mLayerHits[layerID];
  
  int nHits0 = layerHits.size();
  int ih=0;
  for( ; ih<nHits0; ih++ ){
    const SearchHit &fh = layerHits[ih];
    const Hit &h = tracker->mHits[fh.id];
    if( h.isUsed ) break;
  }

  int nHits = ih;

  for( ; ih<nHits0; ih++ ){
    const SearchHit &fh = layerHits[ih];
    const Hit &h = tracker->mHits[fh.id];
    if( !h.isUsed ) layerHits[nHits++] = fh;
  }

  layerHits.resize(nHits);
}



void SearchLayer::Create( Tracker *tracker, int layerID, double searchWindowPhi, double searchWindowT, bool makeClones )
{
  //cout<<"Layer create : layer "<<layerID<<endl;
  mTracker = tracker;  
  mLayerID = layerID;
  
  constexpr bool prn = 0;

  //prn = (layerID<4); //SG!!!

  std::vector<SearchHit> &layerHits = mTracker->mLayerHits[mLayerID];

  Layer &geoLayer = Geo::layers[layerID];   
  
  // find detector geometry from hits
  
  mType = Geo::layers[layerID].type;
  mVolume = Geo::layers[layerID].volume;
  mLayer = Geo::layers[layerID].layer;
  mTimeMark = 0;

  mTmin = geoLayer.tMin;
  if( mTmin<0 ){
    mTmin*=1.01;
  } else {
    mTmin*=0.99;
  }
  mTmax = geoLayer.tMax*1.01;
  mR = geoLayer.r;
  mZ = geoLayer.z;
  mPitchU=geoLayer.pitchU;
  mPitchV=geoLayer.pitchV;
 
  int nHits = layerHits.size();

  // add some margin for T and PI, that exact max values stays inside the intervals
  mTmax+=.001;
  mPhiMin = -PI;
  double phiMax = PI+ 0.00001;
  int binsMax = (int) sqrt(nHits/4);
  if( binsMax < 1 ) binsMax=1;
  mNbinsPhi = (searchWindowPhi>1.e-6)  ?int( (phiMax - mPhiMin)/searchWindowPhi ) :binsMax;
  mNbinsT   = (searchWindowT>1.e-6)  ?int( (mTmax - mTmin)/searchWindowT ) :binsMax;    
  if( mNbinsPhi < 1 ){
    //cout<<"mNbinsPhi = "<<mNbinsPhi<<endl;
    //exit(-1);
    mNbinsPhi = 1;
  }
  if( mNbinsT < 1 ){
    //cout<<"mNbinsT = "<<mNbinsT<<endl;
    //cout<<"layerID = "<<layerID<<" TMin = "<<mTmin<<" tmax = "<<mTmax<<" searchWindowT = "<<searchWindowT <<" nHits= "<<nHits<<endl;
    //exit(-1);
    mNbinsT = 1;
  }
  if( mNbinsPhi > binsMax ) mNbinsPhi = binsMax;
  if( mNbinsT   > binsMax ) mNbinsT   = binsMax;

  mBinPhi = (phiMax - mPhiMin)/mNbinsPhi;
  mBinT = (mTmax - mTmin)/mNbinsT;

  mBinPhiInv = 1./mBinPhi;
  mBinTInv = 1./mBinT;

  mNbinsPhi+=2;
  mNbinsT+=2;

  if(prn) cout<<"fill layer hits.. layer "<<layerID<<":"<<endl;
  if(prn) cout<<" mNbinsPhi "<<mNbinsPhi<<" mBinPhi "<<mBinPhi<<endl;
  // 

  mNbinsTotal = mNbinsPhi*mNbinsT;

  mBins.resize(mNbinsTotal);  

  for( int ib=0; ib<mNbinsTotal; ib++ ){
    mBins[ib].firstHit=0;
    mBins[ib].nHits=0;
  }

#ifdef SGSafeRun
  for( uint iih=0; iih<layerHits.size(); iih++ ){
    const SearchHit &fh = layerHits[iih];
    const Hit &h = mTracker->mHits[fh.id];
    
    if( h.layerID != layerID ||  h.isUsed!=0 ) {
      cout<<"Create layer: wrong indexing!!! mark 1"<<endl;
      exit(-1);
    }
    int it = 1 + int ( (fh.t - mTmin)*mBinTInv ) ;
    int iphi = 1 + int ( (fh.phi - mPhiMin )*mBinPhiInv );
    
    if( it<1 || it>=mNbinsT-1 || iphi<1 || iphi>= mNbinsPhi-1 ){
      cout<<"fit layer: wrong indexing: mark 2"<<endl;
      cout<<" it "<<it<<" of "<<mNbinsT<<" iphi "<<iphi<<" of "<<" "<<mNbinsPhi<<endl;
      cout<<"phi "<<fh.phi<<" +pi="<< fh.phi + PI <<endl;
      cout<<"t "<<fh.t<<" tmin "<<mTmin<<" tmax "<<mTmax<<" tbin "<<mBinT<<endl;
      cout<<"it float = "<<(fh.t - mTmin)*mBinTInv<<endl;
      exit(0);
    }
  }
#endif
  
  for( uint iih=0; iih<layerHits.size(); iih++ ){
    const SearchHit &h = layerHits[iih];
    int it = 1 + int ( (h.t - mTmin)*mBinTInv ) ;
    int iphi = 1 + int ( (h.phi - mPhiMin )*mBinPhiInv );    
    int binInd = it*mNbinsPhi + iphi;
    mBins[binInd].nHits++;
  }
    
  if( makeClones ){ // duplicate hits at the phi edges (phi=+-PI)    
    for( int it=0; it<mNbinsT; it++ ){
      int binTf = it*mNbinsPhi;
      int binTl = binTf + mNbinsPhi - 1;
      mBins[binTf].nHits = mBins[binTl-1].nHits;
      mBins[binTl].nHits = mBins[binTf+1].nHits;
    }
  }
  
  nHits = 0;
  for( int ib=0; ib<mNbinsTotal; ib++ ){
    SearchLayerBin &bin = mBins[ib];
    bin.firstHit = nHits;
    nHits += bin.nHits;
    bin.nHits = 0;
  }
   
  // store hits

  mFitHits.resize( nHits );

  for( uint ih=0; ih<layerHits.size(); ih++ ){
    const SearchHit &h = layerHits[ih];

    int it = 1 + int ( (h.t - mTmin)*mBinTInv ) ;
    int iphi = 1 + int ( (h.phi - mPhiMin )*mBinPhiInv );
    
    int binT0 = it*mNbinsPhi;
    int binInd = binT0 + iphi;

    SearchLayerBin &bin = mBins[binInd];
    int ind = bin.firstHit + bin.nHits;
    bin.nHits++;
    mFitHits[ind] = h;
  }

  if( makeClones ){ // duplicate hits at the phi edges (phi=+-PI)    
    for( int it=0; it<mNbinsT; it++ ){
      int iBinTf = it*mNbinsPhi;
      int iBinTl = iBinTf + mNbinsPhi - 1;

      SearchLayerBin &binf  = mBins[iBinTf];
      SearchLayerBin &binff = mBins[iBinTl-1];

      for( int i=0; i<binff.nHits;  i++ ) {
	SearchHit &h = mFitHits[binf.firstHit + i];
	h = mFitHits[binff.firstHit + i];
	if( h.phi > PI - mBinPhi  ){
	  if( h.phi > 0 ) h.phi = -2*PI + h.phi;
	}	
      }
      binf.nHits = binff.nHits;

      SearchLayerBin &binl  = mBins[iBinTl];
      SearchLayerBin &binll = mBins[iBinTf+1];
      for( int i=0; i<binll.nHits;  i++ ) {
	SearchHit &h = mFitHits[binl.firstHit + i];
	h = mFitHits[binll.firstHit + i];
	if( h.phi < -PI + mBinPhi ){
	  if( h.phi < 0) h.phi = 2*PI+h.phi;	
	}
      }
      binl.nHits = binll.nHits;
    }
  }
  
#ifdef SGSafeRun
  int nHits2=0;
  for( int ib=0; ib<mNbinsTotal; ib++ ){
    SearchLayerBin &bin = mBins[ib];
    if( bin.firstHit!= nHits2){
      cout<<" wrong fit layer indexing: mark 3 "<<endl;
      exit(0);
    }
    nHits2+=bin.nHits;
  }
  
  if( nHits2!=nHits ){
    cout<<" wrong fit layer indexing: mark 4 "<<endl;
    exit(0);
  }
#endif
  
}



