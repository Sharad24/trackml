#ifndef SearchLayer_H
#define SearchLayer_H

#include "util.h"
#include <iostream>
#include "TrackModelPhysical.h"
#include "Cuts.h"
#include "DataStructures.h"

struct Tracker;

using namespace std;


struct SearchHit{
#ifdef SGSafeRun
  SearchHit() :phi(0.f), t(0.f), id(0) {}
#endif  
  float phi;
  float t;  
  int id;
  void Print() const;  
};

struct SearchLayerBin{
  int firstHit;
  int nHits;
};

struct SearchLayer
{
  // phi + t (== z or r) coordinates

  static void CleanUp( Tracker *tracker, int layerID );

  void Create( Tracker *tracker, int layerID, double searchWindowPhi, double searchWindowT, bool makeClones );    

  int getBinPhi( double phi ) const { // get bin for phi with searchWindowPhi/2 marging
    int bin = int ( 0.5 + (phi - mPhiMin) * mBinPhiInv );
#ifdef SGSafeRun
    if( bin < 0 || bin >= mNbinsPhi-1 ){
      cout<<"Wrong bin Phi : phi "<<phi<<" bin "<< bin <<" out of "<< mNbinsPhi <<endl;
      return -1;//exit(0);
    }
#endif  
    return bin;
  }
  
  int getBinT( double t ) const {// get bin for t with searchWindowT/2 marging
    int bin = int ( 0.5 + (t-mTmin) * mBinTInv );
    // TODO: can be removed later

#ifdef SGSafeRun
    if( bin < 0 ) bin=0;
    if( bin > mNbinsT-2 ) bin = mNbinsT-2;

    if( bin < 0 || bin > mNbinsT-2 ){
      cout<<"Wrong bin T "<< bin <<" out of "<< mNbinsT <<endl;
      cout<<" t "<<t<<" tmin "<<mTmin<<" mBinTInv "<<mBinTInv<<" mNbinsT "<<mNbinsT<<endl;
      Print();
      return -1;//exit(0);
    }
#endif
    return bin;
  }

  int getBin( double phi, double t ) const {
    int it = getBinT(t);
    int iphi = getBinPhi( phi );
#ifdef SGSafeRun
   if( it<0 || iphi<0 ) return -1;
#endif
    return it*mNbinsPhi + iphi;
  }

  void Print() const;

  Tracker *mTracker = 0;
  bool mIsOn;
  int mType;
  int mVolume;
  int mLayer;
  int mLayerID;
  double mBinPhi;
  double mBinT;
  double mBinPhiInv;
  double mBinTInv;
  int mNbinsPhi;
  int mNbinsT;
  int mNbinsTotal;
  double mTmin;
  double mTmax;
  double mPhiMin;

  double mR;
  double mZ;
  
  unsigned int mTimeMark;
  
  double mPitchU=0;
  double mPitchV=0;

  vector<SearchLayerBin> mBins;
  vector<SearchHit> mFitHits;
};


class SearchLayerArea
{
 public:
  SearchLayerArea() = default;
  void init( const SearchLayer &l, double phi, double t );
  int currentHitIndex() const;
  const SearchHit &currentHit() const;
  bool next();

 private:
  const SearchLayer *mLayer = 0;
  int mHit = 0;
  int mHitEnd = 0;
  int mHitStart2 = 0;
  int mHitEnd2 = 0;
};


inline void SearchLayerArea::init( const SearchLayer &layer, double phi, double t )
{
  // TODO: store nHits in 2 bins in the fit layer
  mLayer = &layer;
  int binT = layer.getBinT( t );
  int binPhi = layer.getBinPhi( phi );
  if( binT<0 || binPhi<0 ){
    mHit = 0;
    mHitEnd = 0;
    mHitStart2 = 0;
    mHitEnd2 = 0;
    return;
  }
  int bin1 = binT*layer.mNbinsPhi + binPhi;
  int bin2 = bin1 + mLayer->mNbinsPhi;

  mHit = layer.mBins[bin1].firstHit -1;
  mHitEnd  = layer.mBins[bin1+1].firstHit + layer.mBins[bin1+1].nHits;

  mHitStart2 = layer.mBins[bin2].firstHit;
  mHitEnd2  = layer.mBins[bin2+1].firstHit + layer.mBins[bin2+1].nHits; 
}

inline int SearchLayerArea::currentHitIndex() const
{
  return mHit;
}

inline const SearchHit &SearchLayerArea::currentHit() const
{
  return mLayer->mFitHits[mHit];
}

inline bool SearchLayerArea::next()
{ 
  mHit++;
  if( mHit < mHitEnd ) return 1;
  if( mHitEnd < mHitEnd2 ){
    mHit = mHitStart2;
    mHitEnd = mHitEnd2;
  }
  return ( mHit < mHitEnd );
}
  
  
  

#endif
