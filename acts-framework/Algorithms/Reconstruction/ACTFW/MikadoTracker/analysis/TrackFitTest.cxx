/*
  root -l statFit.root

  TNtuple *fit = (TNtuple*) gDirectory->FindObjectAny("fit");

  //fit->Draw("dv/(dS*matThetaV*6.696)","w>0&&fabs(dv/(dS*matThetaV))/6.696<10.");
  fit->Draw("dv/errV","w>0&&fabs(dv/errV)<10.");

  "part:nhits:pt:p:w" + ":partR:partZ:partRL:partZL:prim"+
  ":nl:vol:layer:clust:extrpt"+
  ":dupl:x:y:z:r" + ":phi:t:dPL:dTL:dPH"+
  ":dTH:duz:dv:fitpt:pickDuz" + ":pickDv:dir:bV:bV0:bV1"+
  ":bV2:bA0:bA1:bA2:bA3"+":cellsU:cellsV:pitchU:pitchV:dvDet"+
  ":duzDet:dS:matThick:matL:matThetaV"+":matThetaUZ:sigmaZ:normUZ:normV:picknormUZ"+
  ":picknormV:q";
*/


#include "../Tracker.h"
#include "util.h"

#include <iostream>
#include "../TrackModelPhysical.h"
#include "../SearchLayer.h"
using namespace std;

#ifdef USEROOT
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#endif

void Tracker::TrackFitTest( const RecoPassParameters &cuts )
{
  if( cuts.baseLayer3>=0 ){
    TrackFitTest( cuts.baseLayer1,
		     cuts.baseLayer2,
		     cuts.baseLayer3  );
  } else {
    TrackFitTest( -1,
		     cuts.baseLayer1,
		     cuts.baseLayer2  );
  }

}



void Tracker::TrackFitTest( int ilayer1, int ilayer2, int ilayer3 )
{
#ifdef XXX_USEROOT
  cout<<"\n-----------------\n track fit test  .."<<endl;

  constexpr int NLayers = Geo::NLayers;

  static int nEvents = 0;
  static TFile *statFile = 0;
  static TNtuple *ntFitTracks = 0;

  if( nEvents == 0 ){
    statFile = new TFile("statFit.root", "RECREATE");
    string vars;
    vars = vars +
      "part:nhits:pt:p:w" + ":partR:partZ:partRL:partZL:prim"+
      ":nl:vol:layer:clust:extrpt"+":dupl:x:y:z:r" +
      ":phi:t:dPL:dTL:dPH"+":dTH:duz:dv:fitpt:pickDuz" +
      ":pickDv:dir:bV:bV0:bV1"+":bV2:bA0:bA1:bA2:bA3"+
      ":cellsU:cellsV:pitchU:pitchV:dvDet"+":duzDet:dS:matThick:matL:matThetaV"+
      ":matThetaUZ:sigmaZ:normUZ:normV:picknormUZ"+":picknormV:q:phiV0:phiV1:phiV2";

    ntFitTracks = new TNtuple("fit","fit tracks", vars.data() );
    ntFitTracks->SetMarkerStyle(8);
    ntFitTracks->SetMarkerSize(0.3);
  }
  nEvents++;

  static int statCreateTries=0;
  static int statCreateErrors=0;
  static int statFitErrors=0;

  int baseLayerIDs[3] = {ilayer1, ilayer2, ilayer3};


  for( unsigned int ipart=0; ipart<mTracker->mParticles.size(); ipart++ ){

    Particle &p = mTracker->mParticles[ipart];
    //if( p.r>1 ) continue;
    //if( p.pt<1 ) continue; //SG!!!
    //if( ipart!=301 ) continue;
    bool prn=0;
    vector<int> layerHits[NLayers];

    double pPrevLayer = p.p;

    bool particleOK = 1;
    for( int il=0; il<NLayers && particleOK; il++ ){
      SearchLayer &layer = mSearchLayers[il];
      double pThisLayer = -1;
      for( uint ih = 0; ih<layer.mFitHits.size(); ih++ ){
	SearchLayerHit &h = layer.mFitHits[ih];
	HitMC &mc = mTracker->mHitsMC[h.id];
	if( mc.partID!=(int)ipart ) continue;
	//if( h.isClone ) continue;
	if( prn) cout<<"l "<<il<<" p "<<mc.p<<endl;
	if( prn) mc.Print();
	if( mc.p > pPrevLayer + 0.001 ){ // skip this particle
	  cout<<"wrong hit ordering in fit test, particle "<<ipart<<" dp "<< mc.p - pPrevLayer<<endl;
	  //cout<<"pPrev = "<<pPrevLayer<<endl;
	  if( prn) p.Print();
	  if( prn) mc.Print();
	  //exit(1);
	  //continue;
	  particleOK=0;
	  break;
	}
	if( mc.p > pThisLayer ) pThisLayer = mc.p; // store p for the earliest hit on the layer
	layerHits[il].push_back( ih );
      }
      if( pThisLayer>0 ) pPrevLayer = pThisLayer;
    }

    if( !particleOK ) continue;

    if( baseLayerIDs[0]>=0 && layerHits[baseLayerIDs[0]].size()==0 ) continue;
    if( baseLayerIDs[1]>=0 && layerHits[baseLayerIDs[1]].size()==0 ) continue;
    if( baseLayerIDs[2]>=0 && layerHits[baseLayerIDs[2]].size()==0 ) continue;


    int nPartLayers = 0;
    for( int il=0; il<NLayers; il++ ){
      if( layerHits[il].size()>0 ) nPartLayers++;
    }

    if( nPartLayers < 3 ) continue;

    for( int direction=0; direction<1; direction++ ){

      // fit forward or backward

      SearchLayerHit h1;
      if( ilayer1>=0 ){
	SearchLayer &layer1 = mSearchLayers[ilayer1];
	h1 = layer1.mFitHits[layerHits[ilayer1][0]];
      } else {
	h1.x=p.x; //SG!!!
	h1.y=p.y;
	h1.z=p.z;
	h1.x=0;
	h1.y=0;
	h1.z=0;

	h1.BzFwd = Geo::OriginBzkG;
 	h1.BzMid = Geo::OriginBzkG;
 	h1.BzBck = Geo::OriginBzkG;
      }

      SearchLayer &layer2 = mSearchLayers[ilayer2];
      SearchLayer &layer3 = mSearchLayers[ilayer3];
      SearchLayerHit h2 = layer2.mFitHits[layerHits[ilayer2][0]];
      SearchLayerHit h3 = layer3.mFitHits[layerHits[ilayer3][0]];

      //Module &module3 = Geo::layers[ilayer3].modules[h3.module % Geo::layers[ilayer3].nModules];

      statCreateTries++;
      TrackModelPhysical t;

      int err = 0;
      if( direction==0 ){
	err = t.createXYLast( h1.x, h1.y, h1.z,
			      h2.x, h2.y, h2.z,
			      h3.x, h3.y, h3.z, h2.BzMid );
      } else {
	err = t.createXYMiddle( h1.x, h1.y, h1.z,
				h2.x, h2.y, h2.z,
				h3.x, h3.y, h3.z, h2.BzMid );
      }

      if( err!=0 ){
	statCreateErrors++;
	continue;
      }

      int iLayerStart=ilayer3+1;
      int iLayerEnd=NLayers;
      int iLayerIncr=1;
      if( direction==1 ){
	iLayerStart=ilayer2-1;
	iLayerEnd=-1;
	iLayerIncr=-1;
      }

      int position=-1;
      //int ilayerLast = ilayer3;
      for( int ilayer=iLayerStart; ilayer!=iLayerEnd; ilayer+=iLayerIncr ){

	if( layerHits[ilayer].size()<=0 ) continue;

	SearchLayer &layer = mSearchLayers[ilayer];
 	const bool radialLayer = (layer.mType == 0);

	// extrapolation to the layer nominal Z or R

	double extrBz = h3.BzFwd;
	if( direction==1 ) extrBz = h1.BzBck;

	double extrT = 0;
	double extrPhi = 0;
	err=0;
	if( layer.mType == 0 ){ // radial layer
	  err = t.getPhiZatR( layer.mR, extrPhi, extrT, extrBz );
	} else {
	  if( ( (layer.mType==-1)&&(t.pz>=0)) ||
	      ( (layer.mType== 1)&&(t.pz<=0))    ){
	    continue;
	  }
	  err = t.getPhiRatZ( layer.mZ, extrPhi, extrT, extrBz );
	}
	if( err!=0 ){
	  //cout<<"Error by track extrapolation to layer r: "<<layer.mR<<" z "<<layer.mZ
	  //  <<" part "<<ipart<<" pt "<<p.pt <<" !!!:  err= "<<err<<endl;
	  //p.Print();
	  //t.Print();
	  statFitErrors++;
	  break;
	}
	position++;

	// save extrapolated values
	TrackModelPhysical te = t;

	{ // add new hit and update track & material variables
	  SearchLayerHit &h = layer.mFitHits[layerHits[ilayer][0]];
	  //Module &module = Geo::layers[ilayer].modules[h.module % Geo::layers[ilayer].nModules];
	  TrackModelPhysical tnew = t;
	  if( direction==0 ){
	    err = tnew.createXYLast( h2.x, h2.y, h2.z,
				     h3.x, h3.y, h3.z,
				     h.x, h.y, h.z, h3.BzMid );
	    if( err!=0 ){
	      statCreateErrors++;
	      continue;
	    }
	    h1 = h2;
	    h2 = h3;
	    h3 = h;
	    t = tnew;
	    //ilayerLast = ilayer;
	  } else {
	    err = tnew.createXYFirst( h.x, h.y, h.z,
				      h1.x, h1.y, h1.z,
				      h2.x, h2.y, h2.z, h1.BzMid);
	    if( err!=0 ){
	      statCreateErrors++;
	      continue;
	    }
	    h3 = h2;
	    h2 = h1;
	    h1 = h;
	    t = tnew;
	    //ilayerLast = ilayer;
	  }
	}

	for( uint ihit=0; ihit<layerHits[ilayer].size(); ihit++ ){
	  SearchLayerHit &h = layer.mFitHits[layerHits[ilayer][ihit]];

	  // find detector module
	  //Module &module = Geo::layers[ilayer].modules[h.module % Geo::layers[ilayer].nModules];

	  //  extrapolated phi,t  at the hit position

	  double extrTH = 0;
	  double extrPhiH = 0;
	  err=0;
	  if( layer.mType == 0 ){ // radial layer
	    err = te.getPhiZatR( h.r, extrPhiH, extrTH, extrBz );
	  } else {
	    err = te.getPhiRatZ( h.z, extrPhiH, extrTH, extrBz );
	  }
	  if( err!=0 ){
	    cout<<"Error by track extrapolation to hit  r: "<<h.r<<" z "<<h.z<<" err "<<err<<endl;
	    statFitErrors++;
	    continue;
	  }

	  double dv=0, duz=0; // resulutions hit <-> extrapolated track, XY extrapolation
	  double dvDet=-1, duzDet=-1; // resulutions hit <-> extrapolated track, extrapolation to hit module
	  double pickDuz=0, pickDv=0; // resulutions duplicates <-> track, fitted with the hit

	  int err1=0;
	  int err2=0;
	  double dS=0; // extrapolation distance along track ex,ey
	  if( radialLayer ){ //SG!!
	    err = te.getDistanceAtXY( h.x, h.y, h.z, duz, dv, extrBz,&dS );
	    err1= t. getDistanceAtXY( h.x, h.y, h.z, pickDuz, pickDv, h2.BzMid );
	    //err2 = te.getDistanceAtXY( module.alpha, h.x, h.y, h.z, duzDet, dvDet, extrBz );
	  } else {
	    err = te.getDistanceAtZ( h.x, h.y, h.z, duz, dv, extrBz,&dS );
	    err1= t.getDistanceAtZ( h.x, h.y, h.z, pickDuz, pickDv, h2.BzMid );
	  }

	  if( err!=0 || err1!=0 || err2!=0 ){
	    cout<<"Error by track chi2 calculation!!!: err= "<<err<<endl;
	    statFitErrors++;
	    continue;
	  }

	  /*
	    00 "part:nhits:pt:p:w" + ":partR:partZ:partRL:partZL:prim"+
	    10 ":nl:vol:layer:clust:extrpt"+":dupl:x:y:z:r" +
	    20 ":phi:t:dPL:dTL:dPH"+":dTH:duz:dv:fitpt:pickDuz" +
	    30 ":pickDv:dir:bV:bV0:bV1"+":bV2:bA0:bA1:bA2:bA3"+
	    40 ":cellsU:cellsV:pitchU:pitchV:dvDet"+":duzDet:dS:matThick:matL:matThetaV"+
	    50 ":matThetaUZ:sigmaZ:normUZ:normV:picknormUZ"+":picknormV:q:phiV0:phiV1:phiV2";
	  */
	  float f[100];
	  for( int j=0; j<100; j++ ) f[j] = 0.;
	  f[ 0] = ipart;
	  f[ 1] = p.hits.size();
	  f[ 2] = p.pt;
	  f[ 3] = p.p;
	  f[ 4] = p.w;

	  f[ 5] = p.r;
	  f[ 6] = p.z;
	  f[ 7] = p.rl;
	  f[ 8] = p.zl;
	  f[ 9] = p.prim;

	  f[10] = nPartLayers;
	  f[11] = Geo::layers[ilayer].volume;
	  f[12] = Geo::layers[ilayer].layer;
	  f[13] = position;
	  f[14] = te.getPtkG();

	  f[15] = ihit;
	  f[16] = h.x;
	  f[17] = h.y;
	  f[18] = h.z;
	  f[19] = h.r;

	  f[20] = h.phi;
	  f[21] = h.t;
	  f[22] = h.phi-extrPhi;
	  f[23] = h.t - extrT;
	  f[24] = h.phi-extrPhiH;
	  if( layer.mType == 0 ){ // radial layer
	    f[25] = h.z - extrTH;
	  } else {
	    f[25] = h.r - extrTH;
	  }
	  f[26] = duz;
	  f[27] = dv;
	  f[28] = t.getPtkG();
	  f[29] = pickDuz;
	  f[30] = pickDv;

	  f[31] = direction;
	  f[32] = p.baseV;
	  f[33] = p.baseV0;
	  f[34] = p.baseV1;
	  f[35] = p.baseV2;
	  f[36] = p.baseA0;
	  f[37] = p.baseA1;
	  f[38] = p.baseA2;
	  f[39] = p.baseA3;


	  Layer &l = Geo::layers[ilayer];

	  if( l.type==0 )
	    f[40] = 0;//h.nCellsU;
	  f[41] = 0;//h.nCellsV;
	  f[42] = l.pitchU;
	  f[43] = l.pitchV;
	  f[44] = dvDet;
	  f[45] = duzDet;
	  f[46] = dS;
	  f[47] = 0;
	  f[48] = 0;
	  f[49] = 0;
	  f[50] = 0;
	  f[51] = 0;

	  f[52] = 0;
	  f[53] = 0;
	  f[54] = 0;
	  f[55] = 0;

	  f[ 56] = p.q;
	  f[ 57] = p.phiV0;
	  f[ 58] = p.phiV1;
	  f[ 59] = p.phiV2;
	  ntFitTracks->Fill(f);
	} // layer hits
      } // layers
    } // fit direction

  } // particles

  cout<<"--------\n"<<"track cnstructorfit test: "<<statCreateErrors<<" errors for "<<statCreateTries<<" created particles  ("<< 100.*statCreateErrors / statCreateTries <<" %)"<<endl;

  statFile->Write();
#endif
}
