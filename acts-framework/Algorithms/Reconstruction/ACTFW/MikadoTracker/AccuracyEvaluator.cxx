/*
 
  root -l statRecEff.root

  TNtuple *eff = (TNtuple*) gDirectory->FindObjectAny("eff");

  eff->Draw("r:z","w>0&&prim&&bV0==1&&pt<.3&&pt>=0.2&&hitRec!=1",""); 

  "part:nhits:pt:p:w" + ":prim:r0:z0:nl:vol" + 
  ":layer:ihit:clust:x:y" + ":z:r:phi:fitpt:partRec" + 
  ":hitRec:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3:partID" +
  ":module:cellsU:cellsV:pitchU:pitchV" + ":px:py:pz:x0:y0" + 
  ":dcaR:dcaZ:track:phiV0:phiV1" + ":phiV2:phiV1L23:phiV2L23:phiV1L45:phiV2L45";

  // partRec: 0 or 1 
  // hitRec: 
  // 0 not reconstructed
  // 1 reconstructed
  // 2 wrongly assigned
  // 3 assigned to ghost
  // 4 clone track

 */

#include "AccuracyEvaluator.h"
#include "Tracker.h"
#include "TrackModelPhysical.h"
#include "util.h"
#include <iostream>
#include <iomanip>
#include "util.h"
#include <iostream>
#include "TrackModelPhysical.h"
#include "DataStructures.h"


#ifdef USEROOT
#include "TFile.h"
#include "TNtuple.h"
#endif

bool prn=0;

using namespace std;




void AccuracyEvaluator::SelectParticles( const RecoPassParameters &par, double phiMin, double phiMax )
{ 
  if( !mTracker->mcFlag ) return;

  int iLayer1 = par.baseLayer1;
  int iLayer2 = par.baseLayer2;
  int iLayer3 = par.baseLayer3;
  
  const int nParticles = (int) mTracker->mParticles.size();  
  int nSelected=0;
  for( int ipart=0; ipart<nParticles; ipart++ ){
    Particle &p = mTracker->mParticles[ipart];         
    p.isSelected = 0;
    if( p.w<=0 ) continue;
    //if( p.hits.size()<4 ) continue;
    //if( p.nLayers<3 ) continue;

    int layerHits[3] = {-1,-1,-1};
    for( uint ih=0; ih<p.hits.size(); ih++ ){
      int id = p.hits[ih];
      int layer = mTracker->mHits[id].layerID;
      if( layer == iLayer1 ) layerHits[0] = id;
      if( layer == iLayer2 ) layerHits[1] = id;
      if( layer == iLayer3 ) layerHits[2] = id;
    }
    
    int ok = ( layerHits[0]>=0 && layerHits[1]>=0 && ( iLayer3 < 0 || layerHits[2]>=0)  ) ;
    if( !ok ) continue;
    Hit h1 = mTracker->mHits[layerHits[0] ];
    Hit h2 = mTracker->mHits[layerHits[1] ];
    
    //if( fabs( h1.t ) < par.l1AbsTMin ) continue;
    //if( fabs( h1.t ) >= par.l1AbsTMax ) continue;
    /*
    double extr2T;
    if( Geo::layers[iLayer2].type==0 ){
      extr2T = Geo::layers[iLayer2].r * h1.z/h1.r;
    } else {
      extr2T = Geo::layers[iLayer2].z * h1.r/h1.z;
    }
    */
    double extr2Phi =  h1.phi;
    //if( fabs( extr2T   - h2.t   ) > par.l2T   ) continue;
    double dPhi = fabs( extr2Phi - h2.phi );

    dPhi = -1;
    int v1 = Geo::getVolume( par.baseLayer1 );
    int v2 = Geo::getVolume( par.baseLayer2 );
    int l1 = Geo::getVolumeLayer( par.baseLayer1 );
    int l2 = Geo::getVolumeLayer( par.baseLayer2 );
    if( v1 == 0 && v2 == 0 ){
      if( l1==0 && l2==1 ) dPhi = p.phiV0;
      if( l1==1 && l2==2 ) dPhi = p.phiV0L12;
    }
    if( v1 == 1 && v2 == 1 ){
      if( l1==0 && l2==1 ) dPhi = p.phiV1;
      if( l1==2 && l2==3 ) dPhi = p.phiV1L23;
      if( l1==4 && l2==5 ) dPhi = p.phiV1L45;
    }
    if( v1 == 2 && v2 == 2 ){
      if( l1==0 && l2==1 ) dPhi = p.phiV2;
      if( l1==2 && l2==3 ) dPhi = p.phiV2L23;
      if( l1==4 && l2==5 ) dPhi = p.phiV2L45;
    }
    // SG!!!
    //dPhi = -1;
   
    if( dPhi > M_PI ) dPhi = fabs( 2*M_PI - dPhi);       
    if( (phiMax>0) && dPhi >= phiMax ) continue;
    if( (phiMin>=0) && dPhi <  phiMin ) continue;    
    //if( dPhi <0 ) continue;
    p.isSelected = 1;        
    nSelected++;
  }
  //  cout<<"N selected: "<<nSelected<<endl;
}




void AccuracyEvaluator::Evaluate( bool doPrint )
{ 
  mEfficiency = 0;
  mPurity = 0;
 
  mHitsMissed = 0;
  mHitsFake = 0;

  mTotalEfficiency = 0.;  
  mTotalRemoved = 0.;
  mTotalFakes = 100.;

  if( !mTracker->mcFlag ) return;
  
  static int nEvents = 0;

#ifdef USEROOT
  static TFile *statFile = 0;
  static TNtuple *ntRecHits = 0;

  if( doPrint && nEvents == 0 ){
    statFile = new TFile("statRecEff.root", "RECREATE");      
    string vars;
    vars = vars + 
      "part:nhits:pt:p:w" + ":prim:r0:z0:nl:ihit" + 
      ":vol:layer:clust:x:y" + ":z:r:phi:fitpt:partRec" + 
      ":hitRec:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3:partID" +
      ":module:cellsU:cellsV:pitchU:pitchV" + ":px:py:pz:x0:y0" + 
      ":dcaR:dcaZ:track:phiV0:phiV0V1:phiV1" + ":phiV2:phiV1L23:phiV2L23:phiV1L45:phiV2L45"+
      ":phiV1L56:phiV2L56:phiV0L12:phiV0L23";
    ntRecHits = new TNtuple( "eff","reconstructed particle hits",vars.data() );
    // 0 not reconstructed
    // 1 reconstructed
    // 2 wrongly assigned
    // 3 assigned to ghost
    // 4 clone track
    ntRecHits->SetMarkerStyle(8);
    ntRecHits->SetMarkerSize(.3);//0.3);
  }
#endif

  nEvents++;


  const int nParticles = (int) mTracker->mParticles.size();  
  int partRecN[nParticles];  
  double partRecW[nParticles];  
  int partRecTr[nParticles];
  bool partSelected[nParticles];
  for( int i=0; i<nParticles; i++ ){
    partRecN[i] = 0;
    partRecW[i] = 0;
    partRecTr[i] = -1;
    partSelected[i] = 0;
  }


  for( int ipart=0; ipart<nParticles; ipart++ ){
    Particle &p = mTracker->mParticles[ipart];
    partSelected[ipart] = p.isSelected;
  }
  

 
  double thisFakesW = 0;
  double thisRecW = 0;
  double thisRecWRemoved = 0;
  double thisRecWMissed = 0;
  double thisRecWFake = 0;
  double thisRecWAllHits=0;
  

  vector<Track> *tracks = &mTracker->mTracks; 

#ifdef USEROOT
  const int nTracks = (int) tracks->size();
  int trackRecPID[nTracks];
  int trackRecStatus[nTracks]; // 0 fake 1 rec 2 clone
#endif

  for( uint itr=0; itr<tracks->size(); itr++ ){
    Track &track = (*tracks)[itr];

#ifdef USEROOT
    trackRecPID[itr] = -1;
    trackRecStatus[itr] = 0;
#endif
    // find the majority of mc indices

    int majorID = -1;
    int majorNhits = 0;
    double majorW=0;
    double totalW=0;
    {
      int trRecNHits[ nParticles ];
      double trRecW[ nParticles ];
      for( int i=0; i<nParticles; i++ ){
	trRecNHits[i]=0;
	trRecW[i] = 0;
      }

      for( int ih=0; ih<track.nHits; ih++ ){
	HitMC &mc = mTracker->mHitsMC[track.hitIDs[ih]];
	if( mc.partID<0 ) majorNhits++;
	else {
	  trRecNHits[mc.partID]++;
	  trRecW[mc.partID]+= mc.w;
	}
	totalW+=mc.w;
      }
      
      for( int i=0; i<nParticles; i++ ){
	if( majorNhits < trRecNHits[i] ){
	  majorNhits = trRecNHits[i] ;
	  majorID = i;
	  majorW = trRecW[i];
	}
      }  
    }    
 
    int fakeType = 3;

    if( track.nLayers>=11 ) fakeType = 11;
    else if( track.nLayers==10 ) fakeType = 10;
    else if( track.nLayers==9 ) fakeType = 9;
    else if( track.nLayers==8 ) fakeType = 8;
    else if( track.nLayers==7 ) fakeType = 7;
    else if( track.nLayers==6 ) fakeType = 6;
    else if( track.nLayers==5 ) fakeType = 5;
    else if( track.nLayers==4 ) fakeType = 4;
    else  fakeType = 3;


    thisRecWRemoved+= totalW;
    statRecWRemoved[fakeType]+= totalW;
    statRecWRemoved[0]+= totalW;

    if( majorID>0 ){
      Particle &p = mTracker->mParticles[majorID];
      thisRecWRemoved+= p.w - majorW;
      statRecWRemoved[fakeType]+= p.w - majorW;
      statRecWRemoved[0]+= p.w - majorW;
      thisRecWMissed += p.w - majorW;
      thisRecWFake += totalW - majorW;
      statRecWMissed += p.w - majorW;
      statRecWFake += totalW - majorW;
    }
    
    if( ( majorID<0 ) || 
	( majorNhits <= 0.5*track.nHits ) ){ // fake track
#ifdef USEROOT
      trackRecPID[itr] = -1; 
      trackRecStatus[itr] = 0;        
#endif
      statFakeW[fakeType]+=totalW;
      statFakeN[fakeType]++;
      statFakeW[0]+=totalW;
      statFakeN[0]++;
      thisFakesW+=totalW;
      thisRecWFake+=totalW;
      continue;  
    }

    Particle &p = mTracker->mParticles[majorID];

    if( majorNhits <= 0.5 * p.hits.size() ){ // too short track
 #ifdef USEROOT
      trackRecPID[itr] = majorID;
      trackRecStatus[itr] = 2;
#endif
       statShortW[fakeType]+=totalW;      
       statShortN[fakeType]++;
       statShortW[0]+=totalW;      
       statShortN[0]++;
       if(partSelected[majorID]){
        statShortSelW[fakeType]+=totalW;      
	statShortSelW[0]+=totalW;      
       }
      continue;
    }
  
    statRecW[fakeType]+=majorW;
    thisRecW+=majorW;

    statRecN[fakeType]++;
    statRecW[0]+=majorW;
    statRecN[0]++;

#ifdef USEROOT
    trackRecPID[itr] = majorID;
    trackRecStatus[itr] = 1;
#endif 
    partRecN[majorID]++; 
    if( partRecW[majorID]==0 || partRecW[majorID] < majorW ){
      partRecW[majorID] = majorW;
      partRecTr[majorID] = itr;
    }
  } // tracks

  
#ifdef USEROOT
  if( doPrint ){
 
    int recHits[mTracker->mHits.size()];
    int hitTracks[mTracker->mHits.size()];
    for( uint i=0; i<mTracker->mHits.size(); i++ ){
      hitTracks[i] = -1;
      recHits[i]=0;
    } 

    // fill hit statuses for particle hits
    
    int nRecHits = 0;
    long double recW = 0;
    for( uint itr=0; itr<tracks->size(); itr++ ){
      // 0 not reconstructed
      // 1 reconstructed
      // 2 wrongly assigned
      // 3 assigned to ghost
      // 4 clone track

      Track &track = (*tracks)[itr];
    
      int pid = trackRecPID[itr];
      int recStatus = trackRecStatus[itr]; // 0 fake 1 rec 2 clone
 
      for( int ih=0; ih<track.nHits; ih++ ){
	int hid = track.hitIDs[ih];
	HitMC &mc = mTracker->mHitsMC[hid];
	hitTracks[hid] = itr;
	if( recStatus == 0 ){// 0 fake track 
	  recHits[hid] = 3;
	  continue;
	}
	if( recStatus == 1 ){// 1 rec track
	  if( mc.partID == pid ){
	    recHits[hid] = 1;
	    nRecHits++;
	    recW+=mc.w;
	  } else {
	    recHits[hid] = 2;
	  }
	  continue;
	}
	if( recStatus == 2 ){// 2 clone track
	  recHits[hid] = 4;
	  continue;
	}
      }
    }
    //cout<<"N Rec hits = "<<nRecHits<<" w = "<<recW<<endl;

   // fill the ntuple

    for( uint ipart=0; ipart<mTracker->mParticles.size(); ipart++ ){          
      Particle &p = mTracker->mParticles[ipart];         
      float f[100];
      for( int j=0; j<100; j++ ) f[j] = 0.;
      /*
	00 "part:nhits:pt:p:w" + ":prim:r0:z0:nl:ihit" + 
	10 ":vol:layer:clust:x:y" + ":z:r:phi:fitpt:partRec" + 
	20 ":hitRec:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3:partID" +      
	30 ":module:cellsU:cellsV:pitchU:pitchV" + ":px:py:pz:x0:y0" +
	40 ":dcaR:dcaZ";
      */ 
      float fitpt = -1;
      if( partRecTr[ipart]>=0 ){
	//Track &track = (*tracks)[partRecTr[ipart]];
	fitpt = -1;//track.pt;
      }
      f[ 0] = ipart;
      f[ 1] = p.hits.size();
      f[ 2] = p.pt;
      f[ 3] = p.p;
      f[ 4] = p.w; 
      f[ 5] = p.prim;
      f[ 6] = p.r;
      f[ 7] = p.z;
      f[ 8] = p.nLayers; 
      f[35] = p.px;
      f[36] = p.py;
      f[37] = p.pz;
      f[38] = p.x;
      f[39] = p.y;
      f[40] = p.dcaR;
      f[41] = p.dcaZ;

      for( uint ih=0; ih<p.hits.size(); ih++ ){
	int hid = p.hits[ih];
	Hit &h = mTracker->mHits[hid];
	Layer &l = Geo::layers[h.layerID];
	f[ 9] = ih;
	f[10] = Geo::getVolume(h.layerID);
	f[11] = Geo::getVolumeLayer(h.layerID);
	f[12] = p.hitClusterIds[ih];
	f[13] = h.x();
	f[14] = h.y();
	f[15] = h.z();
	f[16] = h.r;
	f[17] = h.phi;
	f[18] = fitpt;
	f[19] = (partRecN[ipart]>0) ?1 :0;
	f[20] = recHits[hid];
	f[21] = p.baseV;
	f[22] = p.baseV0;
	f[23] = p.baseV1;
	f[24] = p.baseV2;
	f[25] = p.baseA0;
	f[26] = p.baseA1;
	f[27] = p.baseA2;
	f[28] = p.baseA3;
	f[29] = p.countID;
	f[30] = 0;//h.module;
	f[31] = 0;//h.nCellsU;
	f[32] = 0;//h.nCellsV;
	f[33] = l.pitchU;
	f[34] = l.pitchV;      
	f[42] = hitTracks[hid];
	f[43] = p.phiV0;
	f[44] = p.phiV0V1;
	f[45] = p.phiV1;
	f[46] = p.phiV2;
	f[47] = p.phiV1L23;
	f[48] = p.phiV2L23;
	f[49] = p.phiV1L45;
	f[50] = p.phiV2L45;
	f[51] = p.phiV1L56;
	f[52] = p.phiV2L56;
	f[53] = p.phiV0L12;
	f[54] = p.phiV0L23;
	ntRecHits->Fill(f);
      }
    }

    // fill the ntuple for fake hits
    {
      float f[100];
      for( int j=0; j<100; j++ ) f[j] = -1.;
      /*
	00 "part:nhits:pt:p:w" + ":prim:r0:z0:nl:ihit" + 
	10 ":vol:layer:clust:x:y" + ":z:r:phi:fitpt:partRec" + 
	20 ":hitRec:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3:partID" +      
	30 ":module:cellsU:cellsV:pitchU:pitchV" + ":px:py:pz:x0:y0" +
	40 ":dcaR:dcaZ";
      */ 
      f[0] = -1;
      for( uint hid=0; hid<mTracker->mHits.size(); hid++ ){
	Hit &h = mTracker->mHits[hid];
	HitMC &mc = mTracker->mHitsMC[hid];
	if( mc.partID>=0 ) continue;
	Layer &l = Geo::layers[h.layerID];
	f[10] = Geo::getVolume(h.layerID);
	f[11] = Geo::getVolumeLayer(h.layerID);
	f[13] = h.x();
	f[14] = h.y();
	f[15] = h.z();
	f[16] = h.r;
	f[17] = h.phi;
	f[20] = recHits[hid];
	f[30] = 0;//h.module;
	f[31] = 0;//h.nCellsU;
	f[32] = 0;//h.nCellsV;
	f[33] = l.pitchU;
	f[34] = l.pitchV;
	f[42] = hitTracks[hid];
	ntRecHits->Fill(f);
      }
    }

    statFile->Write();  
  } // doPrint
#endif

  {
    double particlesWCut = 0;
    double particlesWRecCut = 0;   
    for( int ipart=0; ipart<nParticles; ipart++ ){
      Particle &p = mTracker->mParticles[ipart];         
      if( !partSelected[ipart]) continue; 
      particlesWCut+=p.w;
      if( partRecN[ipart] > 0 ){     
	particlesWRecCut+=partRecW[ipart];
      }
      if( prn ){ 
	if( partRecN[ipart] ==0 ){ 
	  cout<<"not found: part "<<ipart<<" nhits "<<p.hits.size()<<endl;	
	} else {
	  Track &track = (*tracks)[partRecTr[ipart]];

	  if( fabs(partRecW[ipart]- p.w)>1.e-10 ){
	    cout<<" hits not found: part "<<ipart<<"  hits "<<track.nHits<<"/"<<p.hits.size()
		<<" w "<<partRecW[ipart]<<"/"<<p.w <<endl;
	  } 
	}
      }
      
    }
    if( particlesWCut>0. ){
      mEfficiency = 100.*particlesWRecCut/particlesWCut;
      mFakes = 100.*thisFakesW/particlesWCut;
    } else {
      mEfficiency = 100.;
      mFakes = 0.; 
    } 

    thisRecWAllHits = thisRecW + thisRecWMissed;
    statRecWAllHits = statRecW[0] + statRecWMissed;

    if( thisRecWAllHits > 0. ){
      mHitsMissed = 100.*thisRecWMissed/thisRecWAllHits;
      mHitsFake = 100.*thisRecWFake/thisRecWAllHits;
    } else {
      mHitsMissed = 0;
      mHitsFake = 0;
    }
 

    mTotalEfficiency = 100.*thisRecW;
    mTotalRemoved = 100.*thisRecWRemoved;
    mTotalFakes = 100.*thisFakesW;

    if( thisRecWRemoved>0. ) mPurity = 100.*thisRecW/thisRecWRemoved;
    else mPurity = 0.;
    
    if( statRecWRemoved[0]>0. ) statPurity = 100.*statRecW[0]/statRecWRemoved[0];
    else  statPurity = 0.;

    long int iEff = (long int ) ( mTotalEfficiency * 1.e12 );
    mTotalEfficiency = iEff * 1.e-12;    
    
  }


  // print efficiency
  
  if( doPrint ){

    //cout<<"total efficiency = "<<std::setprecision(30)<<mTotalEfficiency<<endl;
    static int statParticlesNCut = 0;
    static double statParticlesWCut = 0;
    static int statParticlesNRecCut = 0;
    static double statParticlesWRecCut = 0;
    static double statParticlesWRecFullCut = 0;

    double statHitsMissed = 0;
    double statHitsFake = 0.;
    if( statRecWAllHits > 0. ){
      statHitsMissed = 100.*(1.-statRecWMissed/statRecWAllHits);
      statHitsFake = 100.*(1.-statRecWFake/statRecWAllHits);
    }

    for( int ipart=0; ipart<nParticles; ipart++ ){
      Particle &p = mTracker->mParticles[ipart];         

      if( !partSelected[ipart]) continue;
 
      statParticlesNCut++;
      statParticlesWCut+=p.w;
      if( partRecN[ipart] > 0 ){
	statParticlesNRecCut++;
	statParticlesWRecCut+=partRecW[ipart];
	statParticlesWRecFullCut+=p.w;
      } else if(0) {
	cout<<"particle "<<ipart<<endl;
	p.Print();
	for( uint ih=0; ih<p.hits.size(); ih++ ){
	  Hit &h = mTracker->mHits[p.hits[ih]];
	  h.Print();
	}
      }
    }

    static int statNTracks=0;
  
    statNTracks+=(int) (*tracks).size();


 
    cout<<"------------- reconstruction efficiency in "<<nEvents<<" events ---------"<<endl;
    cout<<"N tracks "<<statNTracks/nEvents<<endl;

    cout<<"Reconstructed w        : "<<100.*statParticlesWRecCut/statParticlesWCut      
	<<"% ("<<statParticlesWRecCut/nEvents<<" out of "<<statParticlesWCut/nEvents<<" weight )"<<endl;
    
    cout<<"HitsMissed w           : "<<statHitsMissed
	<<"% ("<<statRecWMissed/nEvents<<" out of "<< statRecWAllHits/nEvents<<")"<<endl;    
    
    cout<<"HitsFake w             : "<<statHitsFake
	<<"% ("<<statRecWFake/nEvents<<" out of "<<statRecWAllHits/nEvents<<")"<<endl;    
    
    cout<<"Purity                : "<<statPurity<<"% ("<<statRecW[0]/nEvents<<" out of "<<statRecWRemoved[0]/nEvents<<")"<<endl;    
   
 
    cout<<"Reconstructed particles: "<<100.*statParticlesNRecCut/statParticlesNCut      
	<<"% ("<<statParticlesNRecCut/nEvents<<" out of "<<statParticlesNCut/nEvents<<" statParticles )"<<endl;
    cout<<"Reconstructed hits w   : "<<100.*statParticlesWRecCut/statParticlesWRecFullCut
	<<"% ("<<statParticlesWRecCut/nEvents<<" out of "<<statParticlesWRecFullCut/nEvents<<" weight )"<<endl;

 
    cout<<"Recs w                : "<<statRecW[0]/nEvents<<" = "<<100.*statRecW[0]/statParticlesWCut<<"% (";
    for( int i=3; i<=11; i++ ) cout <<i<<":"<<100.*statRecW[i]/statParticlesWCut<<"% ";
    cout<<") out of w "<<statParticlesWCut/nEvents<<endl;
        
    double tmp[12];
    for( int i=0; i<=11; i++ ) tmp[i] = (statRecWRemoved[i]>1.e-8) ?1./statRecWRemoved[i] :0;
    cout<<"Recs purity            : "<<100.*statRecW[0]*tmp[0]<<"% (";
    for( int i=3; i<=11; i++ ) cout <<i<<":"<<100.*statRecW[i]*tmp[i]<<"% ";
    cout<<") out of w "<<statParticlesWCut/nEvents<<endl;

    cout<<"Fakes w                : "<<statFakeW[0]/nEvents<<" = "<<100.*statFakeW[0]/statParticlesWCut<<"% (";
    for( int i=3; i<=11; i++ ) cout <<i<<":"<<100.*statFakeW[i]/statParticlesWCut<<"% ";
    cout<<") out of w "<<statParticlesWCut/nEvents<<endl;

    cout<<"Shorts w               : "<<statShortW[0]/nEvents<<" = "<<100.*statShortW[0]/statParticlesWCut<<"% (";
    for( int i=3; i<=11; i++ ) cout <<i<<":"<<100.*statShortW[i]/statParticlesWCut<<"% ";
    cout<<") out of w "<<statParticlesWCut/nEvents<<endl;

    cout<<"Selected Shorts w     : "<<statShortSelW[0]/nEvents<<" = "<<100.*statShortSelW[0]/statParticlesWCut<<"% (";
    for( int i=3; i<=11; i++ ) cout <<i<<":"<<100.*statShortSelW[i]/statParticlesWCut<<"% ";
    cout<<") out of w "<<statParticlesWCut/nEvents<<endl;

    cout<<"Recs n                : "<<statRecN[0]/nEvents<<" (";
    for( int i=3; i<=11; i++ ) cout <<i<<":"<<statRecN[i]/nEvents<<" ";
    cout<<") out of "<<statNTracks/nEvents<<" tracks:"<<endl;

    cout<<"Fakes n                : "<<statFakeN[0]/nEvents<<" (";
    for( int i=3; i<=11; i++ ) cout <<i<<":"<<statFakeN[i]/nEvents<<" ";
    cout<<") out of "<<statNTracks/nEvents<<" tracks:"<<endl;

    cout<<"Shorts n               : "<<statShortN[0]/nEvents<<" (";
    for( int i=3; i<=11; i++ ) cout <<i<<":"<<statShortN[i]/nEvents<<" ";
    cout<<") out of "<<statNTracks/nEvents<<" tracks:"<<endl;
  }

}


