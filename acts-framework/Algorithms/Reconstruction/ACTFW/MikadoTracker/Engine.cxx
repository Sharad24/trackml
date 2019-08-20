#include "Engine.h"
#include <iostream>
#include "TrackModelPhysical.h"
#include "SearchLayer.h"
#include "Tracker.h"

constexpr bool prn=0;

using namespace std;


void Engine::operator()( int iThread, int nThreads, int zSide )
{  
  mTracks.clear();

  const RecoPassParameters &par = *mTracker->mRecoPassParameters;

  // two selection layers
    
  SearchLayer &baseLayer1 = mTracker->mSearchLayers[par.baseLayer1]; 
  SearchLayer &baseLayer2 = mTracker->mSearchLayers[par.baseLayer2];
  
  // mc quality control
#ifdef SGUseMCInReco
  int nTracks2Good = 0;
  int *mcParticlesNRec2=nullptr;
  if( mTracker->mcFlag ){
    mcParticlesNRec2 = new int[mTracker->mParticles.size()];  
    for( uint i=0; i<mTracker->mParticles.size(); i++ ){
      mcParticlesNRec2[i] = 0;
    }
  }
#endif

  int nTracks2 = 0;

  long int statN1=0;
  long int statN2=0;

  mStatFitErrors = 0;

  // init combinatorial tree 

  mTree[0].iLayer = par.baseLayer1;
  mTree[0].localGap = 0;
  mTree[0].totalGap = 0;
  mTree[0].nMaterialsCrossed = 0;
  mTree[0].totalChi2=0;
  //mTree[0].iHit1=-1;
  mTree[0].hitID=-1;

  mTree[1].iLayer = par.baseLayer2;
  mTree[1].localGap = 0;
  mTree[1].totalGap = 0;
  mTree[1].nMaterialsCrossed = 0;
  mTree[1].totalChi2=0;
  //mTree[1].iHit1=-1;
  mTree[1].hitID=-1;

  // init current tracklet

  mCurrentTrack.nLayers = 2;
  mCurrentTrack.nHits = 2;
  mCurrentTrack.chi2 = 0;
  mCurrentTrack.hitIDs[0] = -1;
  mCurrentTrack.hitIDs[1] = -1;

  bool doLine = ( par.model2 == 1);
  // selection of the first hit


  for( uint ih1=iThread; ih1<baseLayer1.mFitHits.size(); ih1+=nThreads ){    
  
    //cout<<"ih1 = "<<ih1<<endl;
    statN1++;

    const SearchHit h1 = baseLayer1.mFitHits[ih1];
    const Hit &hh1 = mTracker->mHits[h1.id];

    // skip cloned hits in first and last phi search bins
    //if( h1.isClone ) continue;
    if( hh1.z()*zSide < 0 ) continue;

#ifdef SGUseMCInReco
    int mcID1 = -1;
    if( mTracker->mcFlag ) mcID1 = mTracker->mHitsMC[h1.id].partID;  	    
#endif    
    //cout<<"mc track "<<mcID1<<endl;
    // cut on T-area 
    double tabs = fabs(h1.t);
    if( tabs < par.l1AbsTMin || tabs > par.l1AbsTMax ) continue;

    // put hit 1 to the branch
    
    //mTree[0].iHit1 = ih1;
    mTree[0].hitID = h1.id;
    mTree[0].hit = hh1.fitHit;
    mCurrentTrack.hitIDs[0] = h1.id;
    mTree[0].nClones=1;

    
    // propagate a straight line (origin<->hit1) to the layer2
    
    double extr2T;
    if( baseLayer2.mType==0 ){
      extr2T = baseLayer2.mR * hh1.z()/hh1.r;
    } else {
      extr2T = baseLayer2.mZ * hh1.r/hh1.z();
    }
   
    // cut on the layer size in T

    if( extr2T < baseLayer2.mTmin - par.l2T || extr2T > baseLayer2.mTmax + par.l2T ) continue;

    // look for the second hit 


    //cout<<" phiSide = "<<phiSide<<endl;

    double extr2Phi =  hh1.phi;

    extr2Phi = hh1.phi;      
      
    SearchLayerArea area;
    area.init( baseLayer2, extr2Phi, extr2T );
    while( area.next() ){
      statN2++;  
	
      //int ih2 = area.currentHitIndex();
      //cout<<"ih2 = "<<ih2<<endl;
      const SearchHit &h2 = area.currentHit();      
      const Hit &hh2 = mTracker->mHits[h2.id];
      //cout<<"m0"<<endl;
      //cout<<"extr t = "<<extr2T<<" hit t = "<<h2.t <<endl;
      // accurate cut on t and phi
      if( fabs( extr2T   - h2.t   ) > par.l2T   ) continue;

      double dPhi = fabs( extr2Phi - h2.phi );

      if( dPhi >= par.l2Phi ) continue;
      if( dPhi <  par.l2PhiMin ) continue;


      for( int i=0; i<Geo::NLayers*2; i++ ){
	//mTree[i].updatedTrack.init();
      }
      mTree[0].updatedTrack.init();
      mTree[1].updatedTrack.init();
      mTree[2].updatedTrack.init();

	
      // angular cuts passed, create a helix or a line
      int err=0;
      if( doLine ){
	err = mTree[1].updatedTrack.createXYLine( hh1.x(), hh1.y(), hh1.z(),
						  hh2.x(), hh2.y(), hh2.z()  );
      } else {
	err = mTree[1].updatedTrack.createXYLast( 0., 0., 0., 
						  hh1.x(), hh1.y(), hh1.z(),
						  hh2.x(), hh2.y(), hh2.z() , hh1.BzMid()  );
      }

      if( err!=0 ){
	//cout<<"Error by tracklet creation "<<" err= "<<err<<endl;
	//h1.Print();
	//h2.Print();	  
	mStatFitErrors++;
	continue;
      }

      nTracks2++;      
       
#ifdef SGUseMCInReco
      // check mc info
      int mcID2=-1;
      if( mTracker->mcFlag ){
	HitMC &mc2 = mTracker->mHitsMC[h2.id];
	mcID2 = -1;
	if( mcID1 >=0 && mc2.partID == mcID1  ){
	  mcID2 = mcID1;
	  nTracks2Good++;
	  mcParticlesNRec2[mcID2]++;
	}	

	if( mTracker->mUseMC && mcID2<0 ){
	  if( !mTracker->isInMix( mc2 ) )  continue;
	}
      }
#endif
      if( mTree[1].updatedTrack.pt < par.l2Pt ) continue;

      //mTree[0].iHit1 = ih1;
      mTree[0].hitID = h1.id;
      mTree[0].hit = hh1.fitHit;
      //mTree[1].iHit1 = ih2;
      mTree[1].hitID = h2.id;
      mTree[1].hit = hh2.fitHit;
      mTree[1].nClones=1;
      mTree[1].totalChi2 = 0;
      mTree[1].localGap = 0;
      mTree[1].totalGap = 0;
      mTree[1].nMaterialsCrossed = 0;

      mCurrentTrack.hitIDs[0] = h1.id;      
      mCurrentTrack.hitIDs[1] = h2.id;            
      mCurrentTrack.nLayers = 2;
      mCurrentTrack.nHits = 2;
      mCurrentTrack.chi2 = 0;
      
      mBestTrack.nHits = 0;
      mBestTrack.nLayers = 0;

	
      // prolongate forward
      {
	int layer3 = (par.baseLayer3>=0) ?par.baseLayer3 :par.baseLayer2+1;		
	int minLayers =  par.minTrackNLayers;
	//cout<<"Reconstruct slice: call prolongateTrack .."<<endl;
	//prolongateTrack( layer3, Geo::NLayers, EndPointFlag::Off, 3 );	  
	prolongateTrack( layer3, Geo::NLayers, EndPointFlag::Off, minLayers ); //SG!!! min nlayers==2
	//cout<<"Reconstruct slice: end of prolongateTrack, best tracklet n hits = :"<<mBestTrack.nHits<<endl;
	if( mBestTrack.nHits < minLayers ) continue;
      }

      bool isHit3 = ( mBestTrack.nHits > 2 );

      if( mBestTrack.nHits < 2 ) continue;
	
      if( isHit3 && mBestTrack.track3.pt < par.l3Pt ) continue;

      // recalculate track parameters for the 2-nd base layer      
 
      // search for hits between 2 and 3
      if( par.baseLayer3>par.baseLayer2+1 ){
	mCurrentTrack = mBestTrack;
	mEndPoint = mBestTrack.hit3;
	mTree[1].updatedTrack = mBestTrack.track3;	  
	prolongateTrack( par.baseLayer2+1, par.baseLayer3, EndPointFlag::On, 0);
      }
      
      if( isHit3 && baseLayer1.mType==0 ){ // recalculate track parameters for the first 2 layers w/o vertex constraint
	
	if( doLine ){
	  err = mTree[1].updatedTrack.createXYLine( hh1.x(), hh1.y(), hh1.z(),
						    hh2.x(), hh2.y(), hh2.z() );
	} else {
	  err = mTree[1].updatedTrack.createXYFirst( hh1.x(), hh1.y(), hh1.z(),
						     hh2.x(), hh2.y(), hh2.z(),
						     mBestTrack.hit3.x, mBestTrack.hit3.y, mBestTrack.hit3.z,
						     hh2.BzMid()  );
	}
	if( err!=0 ){
	  //cout<<"Error by tracklet creation mark 1 "<<" err= "<<err<<endl;
	  //h1.Print();
	  //h2.Print();
	  //mBestTrack.hit3.Print();
	  mStatFitErrors++;
	  continue;
	}	
      }

      // pick up duplicates at 1 & 2      
      {
	// cout<<"pick up duplicates at base layers.."<<endl;
	mCurrentTrack = mBestTrack;
	mTree[0].updatedTrack = mTree[1].updatedTrack;
	pickUpDuplicates(0);
	pickUpDuplicates(1);
	mBestTrack = mCurrentTrack;
      }

      // search for hits between 1 and 2
      if( par.baseLayer2>par.baseLayer1+1 ){
	mCurrentTrack = mBestTrack;
	mEndPoint = hh1.fitHit;
	//prolongateTrack(par.baseLayer1+1, par.baseLayer2,  EndPointFlag::On, 0);	  
	prolongateTrack( par.baseLayer2-1, par.baseLayer1,  EndPointFlag::On, 0);
      }
      //if( mBestTrack.nHits < 4 ) continue;
      
      // prolongate backwards
      if( par.baseLayer1>0 ){
	mCurrentTrack = mBestTrack;
	mTree[0].hit = hh2.fitHit;
	mTree[1].hit = hh1.fitHit;	  
	if( baseLayer1.mType==0 ){
	  prolongateTrack( par.baseLayer1-1, -1, EndPointFlag::Off, 0);
	} else {	    
	  mEndPoint.x = 0;
	  mEndPoint.y = 0;
	  mEndPoint.z = 0;
	  prolongateTrack( par.baseLayer1-1, -1, EndPointFlag::On, 0); //SG!!!
	}
	mTree[0].hit = hh1.fitHit;
      }       

      if( mBestTrack.nHits < 4 ) continue;
      if( mBestTrack.nHits < par.minTrackNLayers ) continue;
      // store the track
      {
	Track t;
	t.nHits = mBestTrack.nHits;
	for( int i=0; i<mBestTrack.nHits; i++) t.hitIDs[i] = mBestTrack.hitIDs[i];
	t.nLayers = mBestTrack.nHits;
	t.chi2 = mBestTrack.chi2;	  
	//mymutex.lock();
	mTracks.push_back(t);
	//mymutex.unlock();
	//Print(t);
      }      
    } // hit 2
  } // hit 1
 
#ifdef SGUseMCInReco
  int nParticlesRec2=0;

  if( mTracker->mcFlag ){
    for( uint i=0; i<mTracker->mParticles.size(); i++ ){
      Particle &p = mTracker->mParticles[i];
      if( p.w>0 && mcParticlesNRec2[i]>0 ) nParticlesRec2++;
    }
  }
  delete[] mcParticlesNRec2;
  
  if( mDoPrint ){
    cout<<"Created "<<nTracks2<<" duplets  ("<<  nTracks2Good <<" good )"<<endl;
    cout<<"Reconstructed "<<nParticlesRec2<<" particles at layer 2"<<endl;
    cout<<"\nfit errors: "<<mStatFitErrors<<endl;    
    cout<<"entries in 1 loop: "<<statN1<<endl;
    cout<<"entries in 2 loop: "<<statN2<<endl;
  }
#endif
}



bool Engine::findNextLayer( bool &isLayerInnerCrossed )
{
  // try to prolongate the track to the next layer

  isLayerInnerCrossed = 0;

  CombinatorialKnot &knot = mTree[mNCombLayers-1];
  TrackModelPhysical &t = knot.updatedTrack;
  
  double Bz = (mLayerIncr>0) ?knot.hit.BzFwd :knot.hit.BzBck;
  //SG!!!
  double knotR2 = t.x0*t.x0+t.y0*t.y0;
  CombinatorialKnot &extrKnot = mTree[mNCombLayers];
  double &extrT   = extrKnot.extrT;
  double &extrPhi = extrKnot.extrPhi;
  
  
  if( mLayerIncr>0 ){
    // new scheme. Code does not work with 3 base layers!!
  
  
    const Layer &geoLayer = Geo::layers[ mLastExtrLayer ];
    int zSide = 0 ;
    if( geoLayer.type == -1 ) zSide = 0;
    else if( geoLayer.type == 1 ) zSide = 1;
    else zSide = (t.pz>=0 );

    //const LayerNeighbour *layerTree = (mLayerIncr>0 ) ?geoLayer.nextLayers[zSide] :geoLayer.prevLayers[zSide];
    const LayerNeighbour *layerTree = geoLayer.nextLayers[zSide];
 
    for( int inext=0; inext<5; inext++ ){

      const LayerNeighbour &next = layerTree[inext];
      int newILayer = next.id;
      if( newILayer<0 || newILayer==mLayerEnd) break;
    
      const SearchLayer &layer = mTracker->mSearchLayers[newILayer ];
      const LayerCuts &cuts = mRecoPassParameters->searchCuts[newILayer];
    
      isLayerInnerCrossed = 0;
    
      // if the layer should be skipped
      if( cuts.cutPhi <= 0. || cuts.cutT <= 0. ) continue;    
      //if( !layer.mIsOn ) continue;
    
      bool isRadialLayer = ( layer.mType == 0 );
      int err = t.getPhiT( isRadialLayer,  layer.mZ, layer.mR, extrPhi, extrT, Bz );    
    
      if( err!=0 ){	
	mStatFitErrors++;
	continue;
      }
    
      // check layer crossing
      if( extrT < layer.mTmin - cuts.cutT || extrT > layer.mTmax + cuts.cutT ) continue;	
    
      // store the layer
      double margin = 0.05*(layer.mTmax-layer.mTmin);
      isLayerInnerCrossed = ( (extrT >= layer.mTmin + cuts.cutT + margin ) && (extrT <= layer.mTmax - cuts.cutT - margin) );
    
      mLastExtrLayer = newILayer;
      extrKnot.hitArea.init( layer, extrPhi, extrT );
      extrKnot.iLayer = mLastExtrLayer;
      //extrKnot.iHit1 = -1;
      extrKnot.hitID = -1;
      extrKnot.nClones = 0;
      mNCombLayers++;
      return 1;
    }    

    return 0;
  }

  // old scheme: inward extrapolation
  
  int newILayer = mLastExtrLayer+=mLayerIncr;
  
  // obligatory layer3
  
  if( mLastExtrLayer == mRecoPassParameters->baseLayer2 && mRecoPassParameters->baseLayer3>=0  ){
    newILayer = mRecoPassParameters->baseLayer3;
    mLayerEnd = newILayer + mLayerIncr;
  }

  for( ; newILayer!=mLayerEnd; newILayer+=mLayerIncr ){    
    
    const LayerCuts &cuts = mRecoPassParameters->searchCuts[newILayer];

    isLayerInnerCrossed = 0;

    // if the layer should be skipped
    if( cuts.cutPhi <= 0. || cuts.cutT <= 0. ) continue;	    
    
    const SearchLayer &layer = mTracker->mSearchLayers[ newILayer ];

    if( !layer.mIsOn ) continue;
    /*
    if(prn){
      int iv = Geo::getVolume(newILayer);
      int il = Geo::getVolumeLayer(newILayer);
      cout<<"try to extr to layer v "<<iv<<" l "<<il<<".."<<endl;
    }
    */
    // skip vertical layer when it is on another z-side 
    if( ((layer.mType==-1)&&(t.pz>=0)) || ((layer.mType== 1)&&(t.pz<=0))  ) continue;   
  
    bool isRadialLayer = ( layer.mType == 0 );
    int err = t.getPhiT( isRadialLayer,  layer.mZ, layer.mR, extrPhi, extrT, Bz );    
    
    //if( err==-1234 && iv>6 ) err = 0;
    if( err!=0 ){	
      /*
      if(prn){
	cout<<"extrapolation error"<<endl;
	//t.Print();
	//cout<<"  Bz = "<<Bz/Geo::CLight<<endl;      
	}*/
      mStatFitErrors++;
      continue;
    }
  
    // check layer crossing
    if( extrT < layer.mTmin - cuts.cutT || extrT > layer.mTmax + cuts.cutT ) continue;	
    
    //if(prn) cout<<"phiT cut passed"<<endl;
 
    // SG
     if(0){
       int iv = Geo::getVolume(newILayer);
       int il = Geo::getVolumeLayer(newILayer);
       double r2 = layer.mR*layer.mR;
       if( !isRadialLayer ) r2 = extrT*extrT;
       if( ( (mLayerIncr>0 && r2<=knotR2)  || (mLayerIncr<0 && r2>=knotR2) ) ){
	 HitMC &mc = mTracker->mHitsMC[ knot.hitID];      
	 cout<<"radius!!! part = "<<mc.partID<<" v "<<iv<<" l "<<il<<endl;
	 t.Print();
	 cout<<"knot r = "<<sqrt(knotR2)<<" extr r = "<<sqrt(r2)<<endl;
	 continue;
       }
     }
    
    // store the layer
    double margin = 0.05*(layer.mTmax-layer.mTmin);
    isLayerInnerCrossed = ( (extrT >= layer.mTmin + cuts.cutT + margin ) && (extrT <= layer.mTmax - cuts.cutT - margin) );

    mLastExtrLayer = newILayer;
    extrKnot.hitArea.init( layer, extrPhi, extrT );
    extrKnot.iLayer = mLastExtrLayer;
    //extrKnot.iHit1 = -1;
    extrKnot.hitID= -1;
    extrKnot.nClones = 0;
    mNCombLayers++;
    return 1;
  }
  return 0;
}


bool Engine::findNextHit()
{  
  CombinatorialKnot &knot0 = mTree[mNCombLayers-3];
  CombinatorialKnot &knot1 = mTree[mNCombLayers-2];
  CombinatorialKnot &knot  = mTree[mNCombLayers-1];
 
  FitHit &hh1 = knot1.hit;
 
  // remove the old hit
  bool firstBranch=1;
  if( knot.hitID >= 0 ){
    mCurrentTrack.nHits -= knot.nClones;
    //knot.iHit1 = -1;
    knot.hitID = -1;
    firstBranch=0;
  }
  knot.nClones = 0;

  bool doBranching = mRecoPassParameters->doBranching;
  bool doLine = (mRecoPassParameters->model == 1);
  if( knot.iLayer==mRecoPassParameters->baseLayer3 ){
    doBranching = mRecoPassParameters->doBranching3;
  }

  const LayerCuts &cuts = mRecoPassParameters->searchCuts[knot.iLayer];

  //  open cuts when there was extra material crossed SG!!: adjust the value

  double cutUZ = ( knot.nMaterialsCrossed == 0 ) ?cuts.cutUZ : cuts.cutUZms;
  double cutV = ( knot.nMaterialsCrossed == 0 ) ?cuts.cutV : cuts.cutVms;

  
  const SearchLayer &layer = mTracker->mSearchLayers[ knot.iLayer ];
  bool isRadialLayer = ( layer.mType == 0 );

  int bestHit = -1;
  double bestChi2 = 1.e10;

#ifdef SGUseMCInReco
  HitMC &mc0 = mTracker->mHitsMC[ mTree[0].hitID];
  if(prn){
    cout<<" extr phi, T = "<< knot.extrPhi<<" "<< knot.extrT<<" n material crossed "<<knot.nMaterialsCrossed<<endl;
  }
#endif

  while( knot.hitArea.next() ){
    const SearchHit &h = knot.hitArea.currentHit();
    //SG!!! if( h.timeMark == layer.mTimeMark ) continue; // already identified as duplicate

#ifdef SGUseMCInReco
    if(mTracker->mUseMC){ 
      HitMC &mc = mTracker->mHitsMC[ h.id];
      if( mc.partID != mc0.partID ){
	if( !mTracker->isInMix(mc) )  continue;
      }
    }
    if(prn){ 
      cout<<"hit area hit: "<<endl;      
      h.Print();
    }
#endif
    // accurate cut on z and phi
    if( fabs( h.t - knot.extrT   ) > cuts.cutT   ) continue;
    if( fabs( h.phi - knot.extrPhi ) > cuts.cutPhi ) continue;
  
    FitHit &hh = mTracker->mHits[ h.id ].fitHit;
   
    // get hit deviations from the trajectory
    double Bz = (mLayerIncr>0)  ?(hh1.BzFwd + hh.BzBck)/2 :(hh1.BzBck + hh.BzFwd)/2;
    double duz=0, dv=0;
    int err = knot1.updatedTrack.getDistanceAt( isRadialLayer, hh.x, hh.y, hh.z, duz, dv, Bz );
    if( err!=0 ){
      mStatFitErrors++;
      continue;
    }
    
    //if(prn) cout<<"hit duz dv "<< duz<<" "<<dv<<" cuts: "<<cutUZ<<" "<<cutV<<endl;
    // strick cuts on deviation from the extrapolated trajectory    
    if( fabs( duz  ) > cutUZ ) continue;
    if( fabs( dv   ) > cutV ) continue;
    
    // hit found
    
    //if( h.timeMark != layer.mTimeMark && !firstBranch ) mNBranches++;
    if( !firstBranch ) mNBranches++;
 
    double chi2 = duz*duz + dv*dv;
    
    if( doBranching || bestChi2 > chi2 ){
      bestChi2 = chi2;
      bestHit = knot.hitArea.currentHitIndex();
    }
    if( doBranching ) break;    
  }// area

  
  if( bestHit>=0 ){ // try to refit the trajectory with the new hit

    //knot.updatedTrack.init();

    knot.hit = mTracker->mHits[layer.mFitHits[bestHit].id].fitHit;

    const FitHit &h0 = knot0.hit;
    const FitHit &h1 = knot1.hit;
    const FitHit &h  = knot.hit;
    const FitHit &e = mEndPoint;

    // TODO: unify these cases 

  
    int err = -1;
    if( doLine ){
      if( mIsEndpointSet==EndPointFlag::On ){
	err = knot.updatedTrack.createXYLine( h.x, h.y, h.z,
					      e.x, e.y,  e.z );
      } else {
	if( mLayerIncr>0 ){
	  err = knot.updatedTrack.createXYLine( h1.x, h1.y, h1.z,						
						h.x,  h.y,  h.z  );
	} else {
	  err = knot.updatedTrack.createXYLine( h.x,  h.y,  h.z ,			      
						h1.x, h1.y, h1.z  );
	}
      }
    } else {
      if( mIsEndpointSet==EndPointFlag::On ){
	if( mLayerIncr>0 ){
	  err = knot.updatedTrack.createXYLast( h1.x, h1.y, h1.z,			      
						h.x, h.y, h.z,
						e.x, e.y, e.z  , h.BzMid );	  
	} else {
	  if( e.x==0. && e.y==0.){

	    err = knot.updatedTrack.createXYLast( e.x, e.y, e.z,			      
						  h.x, h.y, h.z,
						  h1.x, h1.y, h1.z, h.BzMid );
	  } else {
	    err = knot.updatedTrack.createXYFirst( e.x, e.y, e.z,			      
						   h.x, h.y, h.z,
						   h1.x, h1.y, h1.z, h.BzMid );
	  }
	}
      } else {
	if( mLayerIncr>0 ){
	  if( fabs( h1.BzMid/Geo::CLight ) > 5. ){	    
	    err = knot.updatedTrack.createXYLast( h0.x, h0.y, h0.z,			      
						  h1.x, h1.y, h1.z,
						  h.x,  h.y,  h.z  , h1.BzMid );	    
	  } else {
	    err = knot.updatedTrack.createXYLastFix( h1.x, h1.y, h1.z,
						  h.x,  h.y,  h.z  , 
						  h1.BzMid,  knot1.updatedTrack.pt, knot1.updatedTrack.q );	    
	  }
	  //if(prn) cout<<"try to create: Bz = "<<h1.BzMid/Geo::CLight<<" err = "<<err<<endl;


	} else {
	  err = knot.updatedTrack.createXYFirst( h.x,  h.y,  h.z ,			      
						 h1.x, h1.y, h1.z,
						 h0.x, h0.y, h0.z, h1.BzMid );	
	}
      }
    }
    if( err!=0 ){
      //if( prn ) cout<<"can not refit with the new hit, err = "<<err<<endl;
      //h0.Print();
      //h1.Print();
      //mEndPoint.Print();
      mStatFitErrors++;
      bestHit = -1; 
    }
  }

  if( bestHit<0 ){ // no new hit found at this layer    
    return 0;
  }

  // store the hit
  
  knot.hitID = layer.mFitHits[bestHit].id;
  knot.nClones = 1;
  knot.totalChi2 = knot1.totalChi2 + bestChi2;  
  if( mCurrentTrack.nHits < gMaxTracktHits){
    mCurrentTrack.hitIDs[mCurrentTrack.nHits++] = knot.hitID;
  }
  // layer.mFitHits[bestHit].timeMark = layer.mTimeMark;
  
  // pick up duplicates now in order to exclude them from the combinatorics
  pickUpDuplicates( mNCombLayers-1 );
  return 1;
}


void Engine::pickUpDuplicates( int iCombLayer )
{  
  CombinatorialKnot &knot  = mTree[iCombLayer];  
  const LayerCuts &cuts = mRecoPassParameters->pickUpCuts[knot.iLayer];

  // if the layer should be skipped
  if( cuts.cutPhi <= 0. || cuts.cutT <= 0. ) return; 

  SearchLayer &fitLayer = mTracker->mSearchLayers[ knot.iLayer ];
  SearchLayer &pickLayer = mTracker->mPickUpLayers[ knot.iLayer ];
  if( !pickLayer.mIsOn ) return;

  bool isRadialLayer = ( fitLayer.mType == 0 );
  double Bz = knot.hit.BzMid;

  // TODO: take the values from the knot hit
  double pickT = 0, pickPhi = 0;      

  /*  
  int err = knot.updatedTrack.getPhiT( isRadialLayer, pickLayer.mZ, pickLayer.mR, pickPhi, pickT, Bz );
  if( err!=0 ){ 
    mStatFitErrors++;
    return;
  }  

  if( iCombLayer>1 ) {
    //pickT = knot.hit.t;
    //pickPhi = knot.hit.phi;
    
     pickT = knot.extrT;
     pickPhi = knot.extrPhi;
  }
  */
  {
    Hit &hh = mTracker->mHits[ knot.hitID ];
    pickT = hh.t;
    pickPhi = hh.phi;
  }
  // check layer crossing
  if( pickT < pickLayer.mTmin - cuts.cutT || pickT > pickLayer.mTmax + cuts.cutT ) return;	 


  SearchLayerArea pickArea;
  pickArea.init( pickLayer, pickPhi, pickT );

#ifdef SGUseMCInReco
  HitMC &mc0 = mTracker->mHitsMC[ mTree[0].hitID];
#endif
  double cutUZ = ( knot.nMaterialsCrossed == 0 ) ?cuts.cutUZ : cuts.cutUZms;
  double cutV = ( knot.nMaterialsCrossed == 0 ) ?cuts.cutV : cuts.cutVms;


  // find one best duplicate
  //int bestID = -1;
  //double bestD = 1.e20;
  while( pickArea.next() ){
    
    const SearchHit &h = pickArea.currentHit();
    if( h.id == knot.hitID ) continue;
 
    FitHit &hh = mTracker->mHits[ h.id ].fitHit;
 
    //const SearchHit &hh = fitLayer.mFitHits[h.id];

 #ifdef SGUseMCInReco
   if(mTracker->mUseMC){ 
      HitMC &mc = mTracker->mHitsMC[ h.id];
      if( mc.partID != mc0.partID ){
	if( !mTracker->isInMix( mc ) )  continue;
      }
    }
#endif

   // accurate cut on z and phi
   if( fabs( h.t - pickT   ) > cuts.cutT   ) continue;

   if( fabs( h.phi - pickPhi ) > cuts.cutPhi ) continue;      

   // get hit deviations from the trajectory
    double duz=0, dv=0;
    int err = knot.updatedTrack.getDistanceAt( isRadialLayer, hh.x, hh.y, hh.z, duz, dv, Bz );
    if( err!=0 ){	
      mStatFitErrors++;
      continue;
    }    

    // strick cuts on deviation from the extrapolated trajectory    
    if( fabs( duz  ) > cutUZ ) continue;
    if( fabs( dv   ) > cutV  ) continue;    
   
    
    //double d = duz*duz + dv*dv;
    
    // store the hit
    //hh.timeMark = fitLayer.mTimeMark;
    if( mCurrentTrack.nHits < gMaxTracktHits ){
      knot.nClones++;
      mCurrentTrack.hitIDs[mCurrentTrack.nHits++] = h.id;
    }
      /*      
    if( d< bestD ){
      bestD = d;
      bestID = h.id;
    }
      */
  }// area  
  /*
  if( bestID>=0 ){ // store the hit
    SearchHit &hh = fitLayer.mFitHits[bestID];
    hh.timeMark = fitLayer.mTimeMark;
    knot.nClones++;
    mCurrentTrack.hitIDs[mCurrentTrack.nHits++] = hh.id;
  }
  */
}



void Engine::prolongateTrack( int layerStart, int layerEnd, EndPointFlag isEndPointSet, int minNLayers )
{  
  // if( prn ) cout<<"PROLONGATE tracklet from l "<<layerStart<<" to "<<layerEnd<<endl;
  mLayerStart = layerStart;
  mLayerEnd = layerEnd;
  mLayerIncr = (mLayerStart <= mLayerEnd) ?1 :-1;
  mMinNLayers = minNLayers;
  mIsEndpointSet = isEndPointSet;
  mNBranches=0;

  mNCombLayers = 2;
  
  mTree[1].totalChi2 = mCurrentTrack.chi2;
  int initialNLayers = mCurrentTrack.nLayers;

  TrackCandidate &tracklet = mCurrentTrack;

  bool doBranching = mRecoPassParameters->doBranching3 || mRecoPassParameters->doBranching;   

  // (optional) loop over combinatorial branches
  while( 1 ){ 
    
    // 1. the main loop
    // the state: tracklet is fitted at layer = tracklet.nLayers-1
    // try to extend it to other layers and fit it there
    CombinatorialKnot &knotStart  = mTree[mNCombLayers-1];

    if( mNCombLayers == 2 ) mLastExtrLayer = mLayerStart - mLayerIncr;
    else mLastExtrLayer = knotStart.iLayer;
  
    int localGap = 0;
    int totalGap = knotStart.totalGap;
    int nMaterialsCrossed = 0;


    while( localGap<=mRecoPassParameters->gap && totalGap<=mRecoPassParameters->gapTotal ){ //SG!!!
      bool isLayerInnerCrossed=0;
      //if(prn) cout<<"find next layer.."<<endl;
      if( !findNextLayer( isLayerInnerCrossed ) ) break;  // no new layer exist
      if( mNCombLayers < 3 ){
	cout<<"bug!! wrong N combinatorial layers!!!"<<endl;
	exit(-1);
      }
      CombinatorialKnot &knot  = mTree[mNCombLayers-1];
      //if(prn) cout<<"new layer: knot "<<mNCombLayers-1<<" layer "<<knot.iLayer<<" inner crossed= "<<isLayerInnerCrossed<<endl;
      //if(prn) cout<<"find next hit.."<<endl;
      int oldErr = mStatFitErrors;
      knot.nMaterialsCrossed = nMaterialsCrossed;
#ifdef  SGSafeRun
      knot.updatedTrack.init( INFINITY ); // SG!!! initialisation must not be needed, but there was a bug 
#endif
      if( !findNextHit() ){ // no hit on this layer, try next layer
	//knot.updatedTrack.init(0); // SG!!! initialisation must not be needed, but there was a bug 	
	nMaterialsCrossed++; // layer material crossed, open the cuts SG!!! adjust the value
	if( isLayerInnerCrossed ){
	  if( mStatFitErrors == oldErr){
#ifdef SGUseMCInReco
	    if( 1||iv==8 ){	      
	      int il = Geo::getVolumeLayer(knot.iLayer);
	      CombinatorialKnot &knot0  = mTree[0];
	      CombinatorialKnot &knot1  = mTree[mNCombLayers-2];
	      HitMC &mc0 = mTracker->mHitsMC[ knot0.hitID];
	      HitMC &mc = mTracker->mHitsMC[ knot1.hitID];
	      if( prn  && mc.partID == mc0.partID ){
		cout<<"hit not found: mc "<<mc.partID<<" v "<<Geo::getVolume(knot.iLayer)<<" l "<<Geo::getVolumeLayer(knot.iLayer)<<endl;
		/*
		cout<<"   extr Phi = "<<knot.extrPhmTree[mNCombLayers-1].i<<" extr T = "<<knot.extrT  <<endl;
		
		cout<<"  layer r = "<<Geo::layers[knot.iLayer].r<<" extr T = "<<knot.extrT<<endl;		
		knot1.updatedTrack.Print();
		cout<<"mc pt "<<mc.pt<<" mc pz "<<mc.pz<<" mc q "<<mc.q<<endl;
		cout<<"  track r = "<<knot1.updatedTrack.r0()<<endl;	      
		cout<<"Bz "<<knot1.hit.BzFwd/Geo::CLight<<endl;
		*/
	      }
	    }
#endif
	    if( Geo::getVolume(knot.iLayer) < mTracker->mFirstNonStopVolume ){
	      localGap++;
	      totalGap++;
	    } 
	    if( 0 && mNCombLayers<=3 ){ // SG!! no gap at the 3-d layer allowed
	      mNCombLayers--; 
	      break;
	    }
	  }
	}
	mNCombLayers--;// continue extention from the old layer
	if( mTree[mNCombLayers].iLayer==mRecoPassParameters->baseLayer3 ) break; // obligatory layer 3
      } else { // hit found, extend further from the new layer

#ifdef SGSafeRun	
	if( !std::isfinite( mTree[mNCombLayers-1].updatedTrack.pz)  ){
	  cout<<"something wrong.. "<<endl;
	  knot.updatedTrack.Print();
	  exit(-1);
	}
#endif
	//cout<<"hit found"<<endl;
#ifdef SGUseMCInReco
	if( prn ){
	  //HitMC &mc0 = mTracker->mHitsMC[ knot0.hit.id];
	  HitMC &mc = mTracker->mHitsMC[ knot.hitID];
	  cout<<"Found hit: "<<"mc "<<mc.partID<<" v "<<iv<<" l "<<il<<" nclones: "<<knot.nClones<<endl;
	  knot.updatedTrack.Print();
	  cout<<" Bz = "<<knot1.hit.BzMid/Geo::CLight<<endl;      
	}
#endif
	mTree[mNCombLayers-1].localGap = localGap;
	mTree[mNCombLayers-1].totalGap = totalGap;
	mTree[mNCombLayers-1].nMaterialsCrossed = nMaterialsCrossed;
	localGap = 0;
	nMaterialsCrossed = 0;
	if( initialNLayers + mNCombLayers - 2 ==3 ){
	  tracklet.track3 = mTree[mNCombLayers-1].updatedTrack;
	  tracklet.hit3 = mTree[mNCombLayers-1].hit;
	}
      }
    }

    bool isStopped =  (localGap > mRecoPassParameters->gap || totalGap > mRecoPassParameters->gapTotal );
     
    // 2. state: can not extend the tracklet anymore. 
    // Store the tracklet
    
    tracklet.chi2 = mTree[mNCombLayers-1].totalChi2;
    tracklet.nLayers = initialNLayers + mNCombLayers-2;
    if( tracklet.nLayers >= mMinNLayers) {
      if( (tracklet.nHits > mBestTrack.nHits) || 
	  ( (tracklet.nHits == mBestTrack.nHits) && (tracklet.chi2*mBestTrack.nLayers < mBestTrack.chi2*tracklet.nLayers) )
	  ){
	if( mRecoPassParameters->keepStopped || !isStopped ){    //SG!!!
	  mBestTrack = tracklet;
	}
      }
    }

    // 3. try next combinatorial branch (optional)
    	
    if( !doBranching ) break;    
    
    // try to update the track with the next hit on its last layer
    
    while(mNCombLayers > 2){
      if( mTree[mNCombLayers-1].iLayer==mRecoPassParameters->baseLayer3 ){
	if( !mRecoPassParameters->doBranching3 ){ mNCombLayers--; continue; }
      } else {
	if( !mRecoPassParameters->doBranching ){ mNCombLayers--; continue; }
      }
      if( !findNextHit() ){ // SG!! add the multiplier
	mNCombLayers--;
      } else break;
    }
    if( mNCombLayers<=2 ) break; // end of the combinatorial tree
  }
  
  if( mTracker->mMaxNBranches>=0 && mNBranches > mTracker->mMaxNBranches ){
    mBestTrack.nHits=0;
  }
}

