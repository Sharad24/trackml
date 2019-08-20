/*
  root -l statField.root

  TNtuple *fld = (TNtuple*) gDirectory->FindObjectAny("fld");

  // 0:part :nhits :pt :p :q 5:w :nl :prim :vol :lay
  //10:layID :mx :my :mz :mr 15:mphi :mt :mpt :mp  :hx
  //20:hy :hz :hr :hphi :ht 25:Bz0 :Bz1 :Bz2 :fitBz0 :fitBz1
  //30:fitBz2

*/

#include "../Tracker.h"
#include "../util.h"
#include <iostream>
#include "../TrackModelPhysical.h"
#include "../Geo.h"
#include "../Cuts.h"
#include <iostream>
#include "PolynomFit.h"

#ifdef USEROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include "TMath.h"
#endif

using namespace std;

static constexpr int NFieldPar=9;

class StatLayer
{
public:
  int id=0;
  int vol=0;
  int layer=0;
  PolynomFit fit[3]; // fwd, mid, bck
  double cField[3][NFieldPar]; // Bz field in kGaus (t,phi) for fwd, fit, bckwd
  double fieldSigma[3]={0,0,0};

  double getField( int region, double phi, double t){
    double t2 = t*t;
    double *c=cField[region];
    return
      (c[0] + c[1]*t + c[2]*t2) +
      (c[3] + c[4]*t + c[5]*t2)*sin(phi) +
      (c[6] + c[7]*t + c[8]*t2)*cos(phi);
  }

  double fieldS2[3]={0,0,0};
  long int fieldS2Entries[3]={0,0,0};
};

#ifdef USEROOT

static StatLayer statLayers[Geo::NLayers];

static void resetLayerFit()
{
 for( int il=0; il<Geo::NLayers; il++){
    StatLayer &sl = statLayers[il];
    for( int ir=0; ir<3; ir++ ){
      sl.fit[ir].Reset( NFieldPar );
    }
 }
}

static void initLayers()
{
  for( int il=0; il<Geo::NLayers; il++){
    StatLayer &sl = statLayers[il];
    Layer &l = Geo::layers[il];
    sl.id = il;
    sl.vol = l.volume;
    sl.layer = l.layer;
    for( int ir=0; ir<3; ir++ ){
      for( int i=0; i<NFieldPar; i++){
	sl.cField[ir][i]=0;
      }
      sl.fieldSigma[ir]=0;
      sl.fieldS2[ir]=0;
      sl.fieldS2Entries[ir]=0;
   }
  }
  resetLayerFit();
}


/*
static void writeLayerField()
{
  ofstream out("geoLayerField.txt");
  if( !out.is_open() ){
    cout<<"analyseGeometry:: Can not open output file"<<endl;
    exit(0);
  }

  for( int il=0; il<Geo::NLayers; il++){
    StatLayer &sl = statLayers[il];
    for( int ir=0; ir<3; ir++){
      sl.fieldSigma[ir] = -1;
      if( sl.fieldS2Entries[ir]>0 ) sl.fieldSigma[ir] = sqrt(sl.fieldS2[ir]/sl.fieldS2Entries[ir]);
      out<<il;
      for( int i=0; i<NFieldPar; i++)  out<<" "<<sl.cField[ir][i];
      out<<" "<<sl.fieldSigma[ir]<<endl;
    }
  }
  out.close();
}
*/

static void readLayerField()
{
  ifstream in("geoLayerField.txt");
  if( !in.is_open() ){
    cout<<"analyseGeometry:: Can not open input file"<<endl;
    exit(0);
  }
  for( int il=0; il<Geo::NLayers; il++){
    StatLayer &sl = statLayers[il];
    for( int ir=0; ir<3; ir++){
      int jl;
      in>>jl;
      if( jl!=il ){
	cout<<"geo field file broken"<<endl;
	exit(1);
      }
      for( int i=0; i<NFieldPar; i++ ) in>>sl.cField[ir][i];
      in>>sl.fieldSigma[ir];
      sl.fieldS2[ir]=0;
      sl.fieldS2Entries[ir]=0;
    }
  }
  in.close();
}


static void fitLayerField()
{
  for( int il=0; il<Geo::NLayers; il++){
    StatLayer &sl = statLayers[il];
    for( int ir=0; ir<3; ir++ ){
      int err = sl.fit[ir].Fit(sl.cField[ir]);
      if( err!=0 ){
	cout<<"Can not fit the field!!!! with "<<sl.fieldS2Entries[ir]<<" measurements, error "<<err<<endl;
	sl.fieldSigma[ir] = -1;
      }
    }
  }
}

static void updateLayerField( Tracker *mTracker, vector<double> hitField[3], double cut )
{
  for( uint ih=0; ih<mTracker->mHits.size(); ih++ ){
    const Hit &h = mTracker->mHits[ih];
    for( int ir=0; ir<3; ir++){
      double mesBz = hitField[ir][ih];
      if( mesBz<-100. ) continue;
      //Layer &l=Geo::layers[h.layerID];
      StatLayer &sl = statLayers[h.layerID];
      double polBz = sl.getField( ir, h.phi, h.t );
      double d = mesBz - polBz;
      bool specialRegion = (h.layerID==35) || ((h.layerID==34)&&ir!=2);
      if( !specialRegion && cut>0. && fabs(d)>cut*sl.fieldSigma[ir] ) continue;
      sl.fieldS2[ir]+=d*d;
      sl.fieldS2Entries[ir]++;
      double s = sin(h.phi);
      double c = cos(h.phi);
      double t = h.t;
      double t2 = t*t;
      double f[NFieldPar] = {1, t, t2, s, t*s, t2*s, c, t*c, t2*c};
      sl.fit[ir].AddMeasurement( f, mesBz );
    }
  }
}
#endif


void Tracker::analyzeField( bool endOfData )
{
#ifdef USEROOT

  cout<<"\n-----------------\n analyse geometry  .."<<endl;

  static int nEvents = 0;
  static TFile *statFile = 0;
  static TNtuple *ntField = 0;

  if( nEvents == 0 ){
    statFile = new TFile("statField.root", "RECREATE");

    ntField = new TNtuple
      ("field","field",
       "part:nhits:pt:p:q:w:nl:prim:vol:lay:layID:mx:my:mz:mr:mphi:mt:mpt:mp:hx:hy:hz:hr:hphi:ht:Bz0:Bz1:Bz2:fitBz0:fitBz1:fitBz2");
    ntField->SetMarkerStyle(8);
    ntField->SetMarkerSize(0.3);
    initLayers();
    readLayerField();
  }
  nEvents++;

  if( endOfData ){
    statFile->Write();
    statFile->Close();
    nEvents=0;
    fitLayerField();
    //writeLayerField();
    return;
  }

  // find magnetic field value at each hit

  int nHits = mHits.size();
  vector<double> vBz[3];
  for( int ir=0; ir<3; ir++){
    vBz[ir].resize(nHits);
    for( int ih=0; ih<nHits; ih++) vBz[ir][ih] =-1000;
  }

  for( unsigned int ipart=0; ipart<mParticles.size(); ipart++ ){
    Particle &p = mParticles[ipart];
    if( !p.prim ) continue;

    for( uint ih=0; ih<p.hits.size(); ih++ ){

      int hid = p.hits[ih];
      Hit &h = mHits[hid];

      // search for the next hit
      int iFwd = -1;
      for( uint ih1=ih+1; ih1<p.hits.size(); ih1++ ){
	int hid1 = p.hits[ih1];
	Hit &h1 = mHits[hid1];
	if( (h1.layerID > h.layerID) && (fabs(h1.z()-h.z())<200.) ){
	  iFwd = ih1;
	  break;
	}
      }
      int iBck=-1;
      { // search backwards
	for( int ih1=ih-1; ih1>=0; ih1-- ){
	  int hid1 = p.hits[ih1];
	  Hit &h1 = mHits[hid1];
	  if( (h1.layerID < h.layerID) && (fabs(h1.z()-h.z())<200.)){
	    iBck = ih1;
	    break;
	  }
	}
      }

      HitMC mc = mHitsMC[hid];

      if( mc.pt>1. || mc.pt<.2 ) continue;

      HitMC mcFwd = mc;
      if( iFwd>=0 ) mcFwd = mHitsMC[p.hits[iFwd]];
      HitMC mcBck = mc;
      mcBck.x=0;
      mcBck.y=0;
      mcBck.z=0;
      if( iBck>=0 ) mcBck = mHitsMC[p.hits[iBck]];

      if( mc.partID!=(int)ipart || mc.hitID!=hid || mc.q!=p.q ||
	  (iFwd>=0&&mcFwd.partID!=(int)ipart) || (iBck>=0&&mcBck.partID!=(int) ipart)  ){
	cout<<"geometry: hit indexation is broken!!"<<endl;
	exit(1);
      }

      // fit forward field
      if( iFwd>=0 ){
	if( mc.px*mcFwd.px+mc.py*mcFwd.py < 0 ){
	  cout<<"a kink"<<endl;
	  break;
	}
	double BzkG;
	int err = TrackModelPhysical::estmateBzkG( mc, mcFwd.x, mcFwd.y, BzkG );
	if( err!=0 ){
	  cout<<"geometry: can not estimate forward Bz. vol "
	      <<Geo::getVolume(h.layerID)<<" layer "<<Geo::getVolumeLayer(h.layerID)<<" err = "<<err<<endl;
	  //h.Print();
	  //mc.Print();
	  //mcFwd.Print();
	} else {
	  vBz[0][hid] = BzkG;
	}
      } // fwd

      // fit backward field
      {
	if( iBck>=0 && (mc.px*mcBck.px+mc.py*mcBck.py < 0 ) ){
	  cout<<"a kink"<<endl;
	  break;
	}
	double BzkG;
	int err = TrackModelPhysical::estmateBzkG( mc, mcBck.x, mcBck.y, BzkG );
	if( err!=0 ){
	  cout<<"geometry: can not estimate backward Bz, err="<<err<<endl;
	  //h.Print();
	  //mc.Print();
	  //mcBck.Print();
	} else {
	  vBz[2][hid] = BzkG;
	}
      } // backward

      // fit middle field
      if( iFwd>=0 ){
	TrackModelPhysical tr;
	int err = tr.createXYLast( mcBck.x, mcBck.y, mcBck.z,
				   mc.x, mc.y, mc.z,
				   mcFwd.x, mcFwd.y, mcFwd.z, Geo::OriginBzkG );
	if( err!=0 ){
	  cout<<"geometry: can not estimate  Bz for fit, err="<<err<<endl;
	} else {
	  if( tr.pt>0.1 ){
	    vBz[1][hid] = (Geo::OriginBzkG/Geo::CLight ) * mc.pt / tr.pt*mc.q*tr.q;
	  }
	}
      } // middle
    }
  }


  // add measurements to the field fit
  if(1){
    double cut = 3;
    updateLayerField( this, vBz, cut );
  }

  // fill Bz ntuple
  for( unsigned int ipart=0; ipart<mParticles.size(); ipart++ ){
    Particle &p = mParticles[ipart];
    for( uint ih=0; ih<p.hits.size(); ih++ ){
      int hid = p.hits[ih];
      Hit &h = mHits[hid];
      HitMC &mc = mHitsMC[hid];
      StatLayer &sl = statLayers[h.layerID];

      //if( vBz[ir][hid]<-100 ) continue;
      float f[100];
      for( int j=0; j<100; j++ ) f[j] = -1.;
      // 0:part :nhits :pt :p :q 5:w :nl :prim :vol :lay
      //10:layID :mx :my :mz :mr 15:mphi :mt :mpt :mp  :hx
      //20:hy :hz :hr :hphi :ht 25:Bz0 :Bz1 :Bz2 :fitBz0 :fitBz1
      //30:fitBz2

      f[ 0] = ipart;
      f[ 1] = p.hits.size();
      f[ 2] = p.pt;
      f[ 3] = p.p;
      f[ 4] = p.q;
      f[ 5] = p.w;
      f[ 6] = p.nLayers;
      f[ 7] = p.prim;
      f[ 8] = Geo::getVolume(h.layerID);
      f[ 9] = Geo::getVolumeLayer(h.layerID);
      f[10] = h.layerID;
      f[11] = mc.x;
      f[12] = mc.y;
      f[13] = mc.z;
      f[14] = mc.r;
      f[15] = mc.phi;
      f[16] = mc.t;
      f[17] = mc.pt;
      f[18] = mc.p;
      f[19] = h.x();
      f[20] = h.y();
      f[21] = h.z();
      f[22] = h.r;
      f[23] = h.phi;
      f[24] = h.t;
      f[25] = vBz[0][hid];
      f[26] = vBz[1][hid];
      f[27] = vBz[2][hid];
      f[28] = sl.getField(0,h.phi, h.t);
      f[29] = sl.getField(1,h.phi, h.t);
      f[30] = sl.getField(2,h.phi, h.t);
      ntField->Fill(f);
    }  // particle hits
    //if( ipart%1000==0 ) statFile->Write();
  }// particles
#endif
}
