
#include "Learning.h"
#include "Tracker.h"

Learning::Learning()
  :mEvents(), mCuts(nullptr)
{    
  init();
}

void Learning::init()
{    
  mEvents.clear();
  mCuts = 0;
  selectAll(1);
  mSkipDead=0;  
  setDefaultCutsClosed();
}

void Learning::selectAll( bool v )
{
  for( int i=0; i<RecoPassParameters::nPar(); i++ ) parameterSelection[i]=v;
}


void Learning::selectVolume( int iVolume, bool v )
{
  constexpr int nLayerParam = LayerCuts::nPar();
  constexpr int nLayerParamTot = Geo::NLayers*nLayerParam;
  const Volume &vol = Geo::volumes[iVolume];
  for( int il=0; il<vol.nLayers; il++){
    int iLayer = vol.layerIDs[il];
    int first  = RecoPassParameters::nGlobalPar() + iLayer*nLayerParam;
    for( int i=0; i<nLayerParam; i++ ){
      parameterSelection[first + i] = v;
      parameterSelection[nLayerParamTot + first + i] = v;
    }
  }
}

void Learning::selectParameter( int iPar, bool v )
{
  parameterSelection[iPar] = v;
}

void Learning::selectParameters( bool v )
{
  for( int i=0; i< RecoPassParameters::nGlobalPar(); i++ ){
    parameterSelection[i] = v;
  }
}

void Learning::selectSearchLayer( int il, bool v )
{
  constexpr int nLayerParam = LayerCuts::nPar();
  for( int i=0; i<nLayerParam; i++ ){
    parameterSelection[RecoPassParameters::nGlobalPar() + il*nLayerParam + i ] = v;
  }
}

void Learning::selectPickUpLayer( int il, bool v )
{
  constexpr int nLayerParam = LayerCuts::nPar();
  constexpr int nLayerParamTot = Geo::NLayers*nLayerParam;
  for( int i=0; i<nLayerParam; i++ ){
    parameterSelection[RecoPassParameters::nGlobalPar() + nLayerParamTot + il*nLayerParam + i ] = v;
  }
}

void Learning::selectSearchLayers( bool v )
{
  constexpr int nLayerParam = LayerCuts::nPar();
  constexpr int nLayerParamTot = Geo::NLayers*nLayerParam;
  for( int i=RecoPassParameters::nGlobalPar(); i< RecoPassParameters::nGlobalPar() + nLayerParamTot; i++ ){
    parameterSelection[i] = v;
  }
}

void Learning::selectPickUpLayers( bool v )
{
  constexpr int nLayerParam = LayerCuts::nPar();
  constexpr int nLayerParamTot = Geo::NLayers*nLayerParam;
  for( int i=RecoPassParameters::nGlobalPar() + nLayerParamTot; i<RecoPassParameters::nPar(); i++ ){
    parameterSelection[i] = v;
  }
}

void Learning::selectPhiT( bool v )
{
  constexpr int nLayerParam = LayerCuts::nPar();
  for( int i=RecoPassParameters::nGlobalPar(); i< RecoPassParameters::nPar(); i+=nLayerParam  ){
    parameterSelection[i+0] = v;
    parameterSelection[i+1] = v;
  }
}

void Learning::selectUV( bool v )
{
  constexpr int nLayerParam = LayerCuts::nPar();
  for( int i=RecoPassParameters::nGlobalPar(); i< RecoPassParameters::nPar(); i+=nLayerParam  ){
    parameterSelection[i+2] = v;
    parameterSelection[i+3] = v;
  }
}

void Learning::selectUVms( bool v )
{
  constexpr int nLayerParam = LayerCuts::nPar();
  for( int i=RecoPassParameters::nGlobalPar(); i< RecoPassParameters::nPar(); i+=nLayerParam  ){
    parameterSelection[i+4] = v;
    parameterSelection[i+5] = v;
  }
}

void Learning::extendSearchLayers( RecoPassParameters &cuts, double v )
{
  for( int il=0; il<Geo::NLayers; il++){
    LayerCuts &c = cuts.searchCuts[il];
    c.cutPhi*= v;
    c.cutT*=v;
    c.cutV*=v;
    c.cutUZ*=v;
    c.cutVms*=v;
    c.cutUZms*=v;
  }
}

void Learning::extendPickUpLayers( RecoPassParameters &cuts, double v )
{
  for( int il=0; il<Geo::NLayers; il++){
    LayerCuts &c = cuts.pickUpCuts[il];
    c.cutPhi*= v;
    c.cutT*=v;
    c.cutV*=v;
    c.cutUZ*=v;
    c.cutVms*=v;
    c.cutUZms*=v;
  }
}

void Learning::StoreEvent( const Tracker &tracker )
{
  Event tmp;
  mEvents.push_back(tmp);
  Event &ev = mEvents.back();
  ev.mHits = tracker.mHits;
  ev.mHitsMC =  tracker.mHitsMC;
  ev.mParticles = tracker.mParticles;
  ev.mTrackCandidates = tracker.mTrackCandidates;
  ev.mTracks = tracker.mTracks;
  for( int i=0; i<Geo::NLayers; i++ ){
    ev.mLayerHits[i] = tracker.mLayerHits[i];
  } 
}

void Learning::RestoreEvent( int iev, Tracker &tracker )
{
  Event &ev = mEvents[iev];
  tracker.mHits = ev.mHits;
  tracker.mHitsMC = ev.mHitsMC;
  tracker.mParticles = ev.mParticles;
  tracker.mTrackCandidates = ev.mTrackCandidates;
  tracker.mTracks = ev.mTracks;
  for( int i=0; i<Geo::NLayers; i++ ){
    tracker.mLayerHits[i] = ev.mLayerHits[i];
  } 
}

void Learning::SetCalibrationObject( RecoPassParameters &cuts, bool skipDead)
{
  mCuts = &cuts;
  mCutsOrig = cuts;
  mCutsBest = cuts; 
  mSkipDead = skipDead;
}

void Learning::StartLearning()
{
  vSelectedParameterIDs.clear();
  for( int i=0; i<RecoPassParameters::nPar(); i++){
    if( !parameterSelection[i] ) continue;

    if( mSkipDead ){    
      if( i >= RecoPassParameters::nGlobalPar() ){
	double cut = 0.;//mDefaultCuts.par( i );
	if( fabs(mDefaultCuts.par( i )) <= 1.e-4 && fabs( mCutsBest.par( i ) - cut) <= 1.e-4 ) continue;      
      }
    }  
    
    vSelectedParameterIDs.push_back( i );
  }
  
  mEffBest = -1.;
  mEffOrig = -1.;
  mEffCurrent = -1.;
  mCurrentParameter = 0;
  
  if( vSelectedParameterIDs.size()>0 ){
    setupNewSample(0,startSample,0);
    if( mCurrentParameter < (int) vSelectedParameterIDs.size() ){
      *mCurrentParameterP = *mBestParameterP; // set back to the original value
      mCurrentParameter = -1;
    }
  }
}

bool Learning::startNewTry()
{ 
  if( mCurrentParameter < (int) vSelectedParameterIDs.size() ) return 1;
  return 0; 
} 

void Learning::setupNewSample( int iPar, int iSample, bool tryStrickCut )
{
  *mCuts = mCutsBest;
  
  for( mCurrentParameter = iPar; mCurrentParameter < (int) vSelectedParameterIDs.size(); mCurrentParameter++ ){

    int parID = vSelectedParameterIDs[mCurrentParameter];
    mCurrentParameterP = &(    mCuts->par( parID ) );
    mBestParameterP    = &(mCutsBest.par( parID ) );

    double cutDef = mDefaultCuts.par( parID );  
    {
      int id = parID - RecoPassParameters::nGlobalPar();
      if( id > 0 ){
	int ip = id % LayerCuts::nPar();
	if( ip == 4 ||  ip == 5 ){ // UVms parameters : set default to UV
	  cutDef = mCutsBest.par( parID-2 );
	  int layerID = id / LayerCuts::nPar();
	  if( layerID > Geo::NLayers ) layerID -= Geo::NLayers;
	  int vol = Geo::getVolume( layerID );
	  if( 0 && mCuts->gapTotal==0 && vol < mCuts->firstNonStopVolume ){ // does not work this way
	    *mCurrentParameterP = cutDef;	    
	    *mBestParameterP = cutDef;	    
	    continue;
	  }
	}
      }
    }
      
    mSample = iSample;
    double val = *mBestParameterP; 
    double parMin = RecoPassParameters::parRound( parID );
 
    double dv = 0.05;
    if( mSample==0 ) dv = 0.2;  
    if( mSample==1 ) dv = 0.05;
    if( mSample==2 ) dv = 0.01;
    if( mSample==3 ) dv = 0.0025;
    dv = fabs( dv*val );
  
    // new 
    dv = 0.1;
    if( mSample==0 ) dv = 1000*parMin;//0.5 *RecoPassParameters::parMaxStep( parID ) ;
    if( mSample==1 ) dv = 100*parMin;
    if( mSample==2 ) dv = 10*parMin;
    if( mSample==3 ) dv = parMin;
    dv*=100; //SG!!!
    dv = fabs( dv );
    if( dv < parMin ){
      dv = parMin;    
    }

    mSampleValues[0] = val - 3*dv;
    mSampleValues[1] = val -   dv;
    mSampleValues[2] = val +   dv;
    mSampleValues[3] = val + 3*dv;  
   
    for( int i=0; i<4; i++ ){
      mSampleValues[i] = round(mSampleValues[i]/parMin)*parMin;
    }

  
    if( tryStrickCut ){
      if( cutDef > 0. ) mSampleValues[0] = cutDef;
      else mSampleValues[3] = cutDef;
    }

    mSampleElement = 0;
    for( int i=0; i<4; i++) mSampleEff[i] = -1;
    *mCurrentParameterP = mSampleValues[ mSampleElement ];
    //cout<<" set parameter to "<<*mCurrentParameterP<<endl;
    break;
  }
}


void Learning::setTryResult( double eff, double info[], int nnInfo )
{
  if( mCurrentParameter >= (int) vSelectedParameterIDs.size() ) return;

  nInfo = nnInfo;
  mEffCurrent = eff;  
  
  if( mCurrentParameter < 0 ){ // first call, store original efficiency & start
    mEffOrig = eff;
    mEffBest = eff;
    for( int i=0; i<nInfo; i++ ){
      infoOrig[i] = info[i];
      infoBest[i] = info[i];
    }
    setupNewSample( 0, startSample, 0 );
    return;
  }
  
  mSampleEff[ mSampleElement ] = eff;
  for( int j=0; j<nInfo; j++ ){
    infoSample[mSampleElement][j] = info[j];
  }

  if( mSampleElement < 3 ){ // try next element  
    mSampleElement++;
    *mCurrentParameterP = mSampleValues[ mSampleElement ];
    //cout<<" set parameter to "<<*mCurrentParameterP<<endl;
   return;
  }

  // sample completed, pick up the best value
  //bool isDifference=0;
  bool isChanged =0;
  bool tryStrickCut=0;
  int parID = vSelectedParameterIDs[mCurrentParameter];  
  double cutDef = mDefaultCuts.par( parID );
  {
    int id = parID - RecoPassParameters::nGlobalPar();
    if( id > 0 ){
      int ip = id % LayerCuts::nPar();
      if( ip == 4 ||  ip == 5 ){ // UVms parameters : set default to UV
	cutDef = mCutsBest.par( parID-2 );
      }
    }
  }

  for( int i=0; i<4; i++ ){
    double eff = mSampleEff[i];
    // when same eff, prefer tighter cut
    if( mEffBest != eff ){
      //isDifference=1;
    }
    
    if( mEffBest > eff )   continue;
    
    if( mEffBest==eff ){
      double d0 = fabs( *mBestParameterP - cutDef );
      double d1 = fabs( mSampleValues[i] - cutDef );
      if( d0 <= d1 ) continue;
      if(  mSample==startSample ){
	tryStrickCut = 1; // try to cut everything away next time      
      }
    }

    mEffBest = eff;
    for( int j=0; j<nInfo; j++ ){
      infoBest[j] = infoSample[i][j];
    }
    *mBestParameterP = mSampleValues[i];
    //if( *mBestParameterP <= 0. ) *mBestParameterP = -1.;
    isChanged=1;  
  }

  // store the best parameter set
  *mCuts = mCutsBest;

  //if( ( fabs( *mBestParameterP - cutDef ) > 1.e-4) && (isChanged || mSample<3 ) ){ // SG!! try this parameter again  
  
  //if(  ( isDifference && (isChanged || mSample<3) ) 
  //   ||(!isDifference && tryStrickCut )              ){ // SG!! try this parameter again  
  if(  isChanged || (mSample<3) ){ // SG!! try this parameter again  
    if( !isChanged ) mSample++;
    setupNewSample( mCurrentParameter, mSample,  tryStrickCut );    
    return;
  }

  // try next parameter
 
  mCurrentParameter++;
  /*
    if( mSkipDead ){    
    for( ; mCurrentParameter < (int) vSelectedParameterIDs.size();  mCurrentParameter++ ){
    int id = vSelectedParameterIDs[mCurrentParameter];
    double cut =  0;//mDefaultCuts.par( id );
    if( fabs( mCuts->par( id ) - cut) > 1.e-4 ) break;      
    }
    }
  */
  if( mCurrentParameter < (int) vSelectedParameterIDs.size() ){
    setupNewSample( mCurrentParameter, startSample, 0 );
  }
}

void Learning::PrintStatus( double eff )
{
  if( mCurrentParameter < 0 || mCurrentParameter >= (int) vSelectedParameterIDs.size() ) return;

  int id = vSelectedParameterIDs[mCurrentParameter];
  cout<<"Learning p "<<mCurrentParameter
      <<"("<<vSelectedParameterIDs[mCurrentParameter]<<")"
      <<" s "<<mSample<<" e "<<mSampleElement;
  cout
    <<" val: "<< mCutsOrig.par( id )
    <<" -> "<< *mBestParameterP
    <<" ->? "<< *mCurrentParameterP
    <<",  QA: "<<mEffOrig<<"(";
  for( int i=0; i<nInfo-1; i++ ) cout << infoOrig[i]<<" ";
  cout << infoOrig[nInfo-1]<<") -> "
       <<mEffBest<<"(";
  for( int i=0; i<nInfo-1; i++ ) cout << infoBest[i]<<" ";
  cout << infoBest[nInfo-1]<<") ->? "<<eff;

  cout<<endl;
  
}


void Learning::extendToMin( RecoPassParameters &p ) const
{
  const double minPhi[Geo::NVolumes] = {0.015000, 0.012000, 0.012000, 
					0.050000, 0.090000, 0.090000,  
					0.050000, 0.300000, 0.300000  };
  const double minT[Geo::NVolumes] = {  1.0000,   0.1200,   0.1200, 
					2.0000,   1.5000,   1.5000,
					1.6000,   4.0000,   4.0000  };
 
  const double minU[Geo::NVolumes] = { 0.3500, 0.1200, 0.1200,
				       1.3000, 1.2000, 1.2000,
				       1.6000, 4.0000, 4.0000 };

  const double minV[Geo::NVolumes] = { 0.2100, 0.1600, 0.1600,
				       1.3000, 2.1000, 2.1000, 
				       2.6000, 2.5000, 2.5000 };

  const double minPu[Geo::NVolumes] = { 0.0500, 0.0220, 0.020, 
					0.1200, 0.1200, 0.1200,
					1.2000, 0.9000, 0.9000 };
  
  const double minPv[Geo::NVolumes] = { 0.006, 0.0060, 0.0060,
					0.0600, 0.0800, 0.0800,
					0.0600, 1.0000, 1.0000 };

  for( int i=0; i<Geo::NLayers; i++ ){
    int v = Geo::getVolume(i);
    //int vl = Geo::getVolumeLayer(i);
    LayerCuts &ls  = p.searchCuts[i];
    LayerCuts &lp  = p.pickUpCuts[i];    

    if( ls.cutPhi < minPhi[v] ) ls.cutPhi = minPhi[v];         
    if( ls.cutT   < minT[v]   ) ls.cutT   = minT[v];    
    if( ls.cutUZ  < minU[v]   ) ls.cutUZ  = minU[v];
    if( ls.cutV   < minV[v]   ) ls.cutV   = minV[v];
    if( ls.cutUZms  < minU[v]   ) ls.cutUZms  = minU[v];
    if( ls.cutVms   < minV[v]   ) ls.cutVms   = minV[v];
    if( lp.cutPhi < minPhi[v] ) lp.cutPhi = minPhi[v];
    if( lp.cutT   < minT[v]   ) lp.cutT   = minT[v];      
    if( lp.cutUZ  < minPu[v]  ) lp.cutUZ  = minPu[v];
    if( lp.cutV   < minPv[v]  ) lp.cutV   = minPv[v];    
    if( lp.cutUZms  < minPu[v]  ) lp.cutUZms  = minPu[v];
    if( lp.cutVms   < minPv[v]  ) lp.cutVms   = minPv[v];    
  }
}


void Learning::selectDead( RecoPassParameters &p, double multSearch, double multFit, double multPick )
{
  const double minPhi[Geo::NVolumes] = {0.015000, 0.012000, 0.012000, 
					0.050000, 0.090000, 0.090000,  
					0.050000, 0.300000, 0.300000  };
  const double minT[Geo::NVolumes] = {  1.0000,   0.1200,   0.1200, 
					2.0000,   1.5000,   1.5000,
					1.6000,   4.0000,   4.0000  };
  
  const double minU[2][Geo::NVolumes] = { { 0.3500, 0.1200, 0.1200,
					    1.3000, 1.2000, 1.2000,
					    1.6000, 4.0000, 4.0000 },
					  { 0.0500, 0.0220, 0.020, 
					    0.1200, 0.1200, 0.1200,
					    1.2000, 0.9000, 0.9000 } };

  const double minV[2][Geo::NVolumes] = { { 0.2100, 0.1600, 0.1600,
					    1.3000, 2.1000, 2.1000, 
					    2.6000, 2.5000, 2.5000 },
					  { 0.006, 0.0060, 0.0060,
					    0.0600, 0.0800, 0.0800,
					    0.0600, 1.0000, 1.0000 } };

 
  constexpr int nLayerParam = LayerCuts::nPar();
  bool *s = parameterSelection + RecoPassParameters::nGlobalPar();
  
  for( int itype=0; itype<2; itype++){
    for( int iL=0; iL<Geo::NLayers; iL++,  s+= nLayerParam ){
      if( itype==0 ){
	if( iL == p.baseLayer1 ) continue;
	if( iL == p.baseLayer2 ) continue;
      }
      int v = Geo::getVolume(iL);
      LayerCuts &l  = (itype==0) ? p.searchCuts[iL] : p.pickUpCuts[iL];      
      double mult = (itype==0) ?multFit :multPick;
      if( l.cutPhi < 1.e-4 ){ s[0] = 1; l.cutPhi = minPhi[v]*multSearch; }
      if( l.cutT   < 1.e-4 ){ s[1] = 1; l.cutT   = minT  [v]*multSearch; }
      if( l.cutV   < 1.e-4 ){ s[2] = 1; l.cutV   = minV[itype][v]*mult;  }
      if( l.cutUZ  < 1.e-4 ){ s[3] = 1; l.cutUZ  = minU[itype][v]*mult;  }    
      if( l.cutVms < 1.e-4 ){ s[4] = 1; l.cutVms = minV[itype][v]*mult;  }      
      if( l.cutUZms< 1.e-4 ){ s[5] = 1; l.cutUZms= minU[itype][v]*mult;  }
    }
  }
}


void Learning::resetPhi( RecoPassParameters &p, double mult ) const
{

  for( int i=0; i<Geo::NLayers; i++ ){
    //int v = Geo::getVolume(i);
    //int vl = Geo::getVolumeLayer(i);
    //if( v!=4 && v!=5 ) continue;
    //if( vl==4 ) continue;
    LayerCuts &ls  = p.searchCuts[i];
    LayerCuts &lp  = p.pickUpCuts[i];
    
    ls.cutPhi *= mult;
    lp.cutPhi *= mult;

  }
}

void Learning::resetT( RecoPassParameters &p, double mult ) const
{

  for( int i=0; i<Geo::NLayers; i++ ){
    //int v = Geo::getVolume(i);
    //int vl = Geo::getVolumeLayer(i);
    //if( v!=4 && v!=5 ) continue;
    //if( vl==4 ) continue;
    LayerCuts &ls  = p.searchCuts[i];
    LayerCuts &lp  = p.pickUpCuts[i];
    
    ls.cutT *= mult;
    lp.cutT *= mult;

  }
}

void Learning::resetUV( RecoPassParameters &p, double mult ) const
{

  for( int i=0; i<Geo::NLayers; i++ ){
    //int v = Geo::getVolume(i);
    //int vl = Geo::getVolumeLayer(i);
    LayerCuts &ls  = p.searchCuts[i];
    LayerCuts &lp  = p.pickUpCuts[i];
    
    ls.cutUZ *= mult;
    ls.cutV *= mult;
    ls.cutUZms *= mult;
    ls.cutVms *= mult;
    lp.cutUZ *= mult;
    lp.cutV *= mult;
    lp.cutUZms *= mult;
    lp.cutVms *= mult;
  }
}




void Learning::setDefaultCuts( const RecoPassParameters &par, double multParam, double multSearch, double multFit, double multPick )
{
  RecoPassParameters &p = mDefaultCuts;
  p = par;

  double multParamInv = 1.-(multParam-1.);
  if( multParamInv<0. ) multParamInv=0.;

  p.l1AbsTMin = par.l1AbsTMin * multParamInv;
  p.l1AbsTMax = par.l1AbsTMax * multParam;
  p.l2PhiMin = par.l2PhiMin * multParamInv;
  p.l2Phi    = par.l2Phi * multParam;
  p.l2T    = par.l2T * multParam;
  p.l2Pt    = par.l2Pt * multParamInv;
  p.l3Pt    = par.l3Pt * multParamInv;


  for( int itype=0; itype<2; itype++){
    for( int i=0; i<Geo::NLayers; i++ ){
      //int v = Geo::getVolume(i);
      const LayerCuts &l0 = (itype==0) ?par.searchCuts[i] :par.pickUpCuts[i];
      LayerCuts &l  = (itype==0) ? p.searchCuts[i] : p.pickUpCuts[i];

      double mult = (itype==0) ?multFit :multPick;
      
      l.cutPhi = l0.cutPhi*multSearch;      
      l.cutT = l0.cutT*multSearch;      
      l.cutV = l0.cutV*mult;      
      l.cutUZ = l0.cutUZ*mult;      
      l.cutVms = l0.cutVms*mult;      
      l.cutUZms = l0.cutUZms*mult; 
    }
  }
}


void Learning::setDefaultCutsClosed()
{
  RecoPassParameters &p = mDefaultCuts;
  p.ID = -1;
  p.name = "closed cuts";
  p.baseLayer1 = -1;
  p.baseLayer2 = -1;
  p.baseLayer3 = -1;
  p.model2 = 0;
  p.model = 0;
  p.doBranching3=0;
  p.doBranching=0;
  p.gap=0;
  p.gapTotal=0;
  p.keepStopped=0;

  p.l1AbsTMin = 100.;
  p.l1AbsTMax = 0.;
  p.l2PhiMin = 3.;
  p.l2Phi = 0.;
  p.l2T = 0.;

  p.l2Pt = 1;
  p.l3Pt = 1;

  for( int i=0; i<Geo::NLayers; i++ ){
    LayerCuts &l  = p.searchCuts[i];
    l.cutPhi = 0.;
    l.cutT = 0.;
    l.cutV = 0.;
    l.cutUZ = 0.;
    l.cutVms = 0.;
    l.cutUZms = 0.;
  }
  for( int i=0; i<Geo::NLayers; i++ ){    
    LayerCuts &l  = p.pickUpCuts[i];
    l.cutPhi = 0.;
    l.cutT = 0.;
    l.cutV = 0.;
    l.cutUZ = 0.; 
    l.cutVms = 0.;
    l.cutUZms = 0.; 
  }
}


void Learning::writeDefaultCuts( const char *file )
{
  ofstream out(file);    
  if( !out.is_open() ){
    cout<<"Cuts:: Can not open output file"<<endl;
    exit(-1);
  }
  out << mDefaultCuts;
  out.close();
}

void Learning::readDefaultCuts( const char *file )
{
  ifstream in(file);    
  if( !in.is_open() ){
    cout<<"Cuts:: Can not open input file"<<endl;
    exit(-1);
  }
  in >> mDefaultCuts;
  in.close();
}

