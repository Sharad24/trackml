/*
  root -l statGeometry.root

  TNtuple *modules = (TNtuple*) gDirectory->FindObjectAny("modules");

  volume,layer,module,
  cx,cy,cz,
  rot_xu,rot_xv,rot_xw,
  rot_yu,rot_yv,rot_yw,
  rot_zu,rot_zv,rot_zw,
  module_t,
  module_minhu,module_maxhu,
  module_hv,
  pitch_u,pitch_v,
  layerID

*/

#include "../Tracker.h"
#include "../util.h"
#include <iostream>
#include <iomanip>
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

#ifdef USEROOT

class StatLayer
{
public:
  Layer geoLayer;
  int id=0;
  int vol=0;
  int layer=0;
  double tMin=1.e10;
  double tMax=1.e-10;
  double r=0;
  double z=0;
  long int nEntriesRZ=0;
};

StatLayer statLayers[Geo::NLayers];

class StatModule
{
public:
  int volume,
    layer,
    module;
  double cx,cy,cz,
    rot_xu,rot_xv,rot_xw,
    rot_yu,rot_yv,rot_yw,
    rot_zu,rot_zv,rot_zw;
  double
  module_t, module_minhu, module_maxhu, module_hv;
  double pitch_u,pitch_v;

  int layerID;
};

vector<StatModule> modules;

void initLayers()
{
  for( int il=0; il<Geo::NLayers; il++){
    StatLayer &sl = statLayers[il];
    Layer &l = sl.geoLayer;
    l = Geo::layers[il];
    sl.id = il;
    sl.vol = l.volume;
    sl.layer = l.layer;
   }
}

void writeLayerSizes()
{
  ofstream out("geoLayerSizes.txt");
  if( !out.is_open() ){
    cout<<"analyseGeometry:: Can not open output file"<<endl;
    exit(0);
  }
  for( int i=0; i<Geo::NLayers; i++){
    StatLayer &sl = statLayers[i];
    Layer &l = sl.geoLayer;
    if( sl.nEntriesRZ>0 ){
      sl.r/=sl.nEntriesRZ;
      sl.z/=sl.nEntriesRZ;
      sl.nEntriesRZ=1;
    }
    l.r = sl.r;
    l.z = sl.z;
    l.tMin = sl.tMin;
    l.tMax = sl.tMax;
    l.WriteGeometry( out );
  }
  out.close();
}

void readLayerSizes()
{
  ifstream in("geoLayerSizes.txt");
  if( !in.is_open() ){
    cout<<"analyseGeometry:: Can not open input file"<<endl;
    exit(0);
  }
  for( int i=0; i<Geo::NLayers; i++){
    StatLayer &sl = statLayers[i];
    Layer &l = sl.geoLayer;
    l.ReadGeometry( in );
    if( l.id!=i ){
      cout<<"geo file broken"<<endl;
      exit(1);
    }
    sl.r = l.r;
    sl.z = l.z;
    sl.tMin = l.tMin;
    sl.tMax = l.tMax;
    sl.nEntriesRZ = 1;
  }
  in.close();
}

void updateLayerSizesRZ( Tracker *mTracker)
{
  for( uint ih=0; ih<mTracker->mHits.size(); ih++ ){
    const Hit &h = mTracker->mHits[ih];
    StatLayer &sl = statLayers[h.layerID];
    sl.r+= h.r;
    sl.z+= h.z();
    sl.nEntriesRZ++;
  }
}

void updateLayerSizesT( Tracker *mTracker)
{
  for( uint ih=0; ih<mTracker->mHits.size(); ih++ ){
    const Hit &h = mTracker->mHits[ih];
    Layer &l=Geo::layers[h.layerID];
    StatLayer &sl = statLayers[h.layerID];

    double t = 0;
    if( l.type==0 ) t = h.z();
    else  t = h.r;
    if( sl.tMin > t ) sl.tMin = t;
    if( sl.tMax < t ) sl.tMax = t;

    t = 0;
    if( l.type==0 ) t = h.z()/h.r*sl.r;
    else  t = h.r/h.z()*sl.z;
    if( sl.tMin > t ) sl.tMin = t;
    if( sl.tMax < t ) sl.tMax = t;

  }
}


static void readModules()
{
  TString fname = "data/detectors.csv";
  ifstream in(fname.Data());
  if( !in.is_open() ){
    cout<<"File "<<fname.Data()<<" does not exist!!"<<endl;
    exit(0);
  }
  char tmpLine[256];
  in.getline(tmpLine,256);
  cout<<tmpLine<<endl;

  string vars;
  vars = vars +
    "volume:layer:module:" +
    "cx:cy:cz:"+
    "rot_xu:rot_xv:rot_xw:"+
    "rot_yu:rot_yv:rot_yw:"+
    "rot_zu:rot_zv:rot_zw:"+
    "module_t:module_minhu:module_maxhu:module_hv:"+
    "pitch_u:pitch_v:"+
    "layerID";

  TNtuple *ntModules = new TNtuple ("modules","modules",vars.data() );

  ntModules->SetMarkerStyle(8);
  ntModules->SetMarkerSize(0.3);

  while (1) {
    float d[22];
    if( !readLine(in,d,21) ) break;
    StatModule m;
    m.volume = (int) d[0];
    m.layer = (int) d[1];
    m.module = (int) d[2];
    m.cx = d[3]*0.1; // convert to [cm]
    m.cy = d[4]*0.1; // convert to [cm]
    m.cz = d[5]*0.1; // convert to [cm]
    m.rot_xu = d[ 6];
    m.rot_xv = d[ 7];
    m.rot_xw = d[ 8];
    m.rot_yu = d[ 9];
    m.rot_yv = d[10];
    m.rot_yw = d[11];
    m.rot_zu = d[12];
    m.rot_zv = d[13];
    m.rot_zw = d[14];
    m.module_t = d[15]*0.1; // convert to [cm]
    m.module_minhu = d[16]*0.1; // convert to [cm]
    m.module_maxhu = d[17]*0.1; // convert to [cm]
    m.module_hv = d[18]*0.1; // convert to [cm]
    m.pitch_u = d[19]*0.1; // convert to [cm]
    m.pitch_v = d[20]*0.1; // convert to [cm]

    m.module = m.module-1;
    m.layer = m.layer/2 -1;

    switch( m.volume )
      {
      case  8: m.volume = 0; break;
      case  7: m.volume = 1; m.layer = 6-m.layer; break;
      case  9: m.volume = 2; break;
      case 13: m.volume = 3; break;
      case 12: m.volume = 4; m.layer = 5-m.layer; break;
      case 14: m.volume = 5; break;
      case 17: m.volume = 6; break;
      case 16: m.volume = 7; m.layer = 5-m.layer; break;
      case 18: m.volume = 8; break;
      default:
	cout<<"Unknown detector volume: "<< m.volume << endl;
	exit(0);
      };

    if( m.layer<0 || m.layer>= Geo::volumes[m.volume].nLayers ){
      cout<<"Unknown detector layer: "<<m.layer<<endl;
      exit(0);
    }
    m.layerID = Geo::volumes[m.volume].layerIDs[m.layer];
    modules.push_back(m);

    d[ 0] = m.volume;
    d[ 1] = m.layer;
    d[ 2] = m.module;
    d[3]*=0.1; // convert to [cm]
    d[4]*=0.1;
    d[5]*=0.1;
    d[15]*=0.1;
    d[16]*=0.1;
    d[17]*=0.1;
    d[18]*=0.1;
    d[19]*=0.1;
    d[20]*=0.1;
    d[21] = m.layerID;
    ntModules->Fill(d);
  }
  ntModules->Write();
}

static void writeModules()
{
  ofstream out("geoLayerModules.txt");
  if( !out.is_open() ){
    cout<<"analyseGeometry:: Can not open output file"<<endl;
    exit(0);
  }

  out<<"layerID,volume,layer,module,alpha,pitchU,pitchV,thick05"<<endl;

  out << std::setprecision(10);

  for( int il=0; il<Geo::NLayers; il++){
    Layer &layer = Geo::layers[il];

    int im0=-1;
    for( uint im=0; im<modules.size(); im++ ){
      StatModule &m = modules[im];
      if( m.layerID == il && m.module==0 ){
	im0 = im;
	break;
      }
    }
    if( im0<0 ){ cout<<"analyseGeometry:: wrong layer indexing "<<endl; exit(-1); }

    if( layer.type!=0 ){ // vertical layer: write all the modules
      for( uint im=im0; im<modules.size(); im++ ){
	StatModule &m = modules[im];
	if( m.layerID != il ) break;
	if( m.module!= (int)im-im0 ){
	  cout<<"analyseGeometry:: wrong module indexing "<<endl; exit(-1);
	}
	double alpha = atan2(m.rot_yv,m.rot_xv)/TMath::Pi()*180.;
	out<<il<<" "<<layer.volume<<" "<<layer.layer<<" "
	   <<m.module<<" "<< alpha<<" "<<m.pitch_u<<" "<<m.pitch_v<<" "<<m.module_t<<endl;
      }
    } // layer type 1

    if( layer.type==0 ){ // radial layer: write one z slice

      StatModule &m0 = modules[im0];
      StatModule &m1 = modules[im0+1];
      if( m1.layerID != il || m1.module!=1 ) { cout<<"analyseGeometry:: wrong layer indexing "<<endl; exit(-1); }

      double Alpha0 = atan2(m0.rot_yw,m0.rot_xw)/TMath::Pi()*180.;
      double Alpha1 = atan2(m1.rot_yw,m1.rot_xw)/TMath::Pi()*180.;
      double dAlpha = Alpha1 - Alpha0;
      if( dAlpha>180 ) dAlpha-=360.;
      if( dAlpha<-180 ) dAlpha+=360.;

      int nSectors = 361./dAlpha;

      for( uint im=im0; im<modules.size(); im++ ){
	StatModule &m = modules[im];
	if( m.layerID != il ) break;
	if( m.module!= (int)im-im0 ){
	  cout<<"analyseGeometry:: wrong module indexing "<<endl; exit(-1);
	}
	double alphaOrig = atan2(m.rot_yw,m.rot_xw)/TMath::Pi()*180.;
	double alphaCalc = Alpha0 + dAlpha*( m.module % nSectors );
	if( alphaCalc>180. ) alphaCalc-=360.;
	if( fabs( (alphaOrig-alphaCalc) )>1.e-2 ){
	  cout<<"analyseGeometry:: something wrong with the calculation of sector angles "<<endl;
	  cout<<"orig angle: "<<alphaOrig<<", calculated angle "<<alphaCalc
	  <<" diff "<<alphaOrig - alphaCalc<<endl;
	  cout<<" vol "<<layer.volume<<" layer "<<layer.layer<<" module "
	      <<m.module<<" alpha "<< alphaOrig<<endl;
	  exit(-1);
	}
	if( m.module < nSectors ){
	  out<<il<<" "<<layer.volume<<" "<<layer.layer<<" "
	     <<m.module<<" "<< alphaOrig<<" "<<m.pitch_u<<" "<<m.pitch_v<<" "<<m.module_t<<endl;
	}
      }
    } // layer type 0
  }
  out.close();
}
#endif



void Tracker::analyzeGeometry( bool endOfData )
{
#ifdef USEROOT
  cout<<"\n-----------------\n analyse geometry  .."<<endl;

  static int nEvents = 0;
  static TFile *statFile = 0;
  static TNtuple *ntGeo = 0;

  if( nEvents == 0 ){
    statFile = new TFile("statGeometry.root", "RECREATE");
    ntGeo = new TNtuple ("geo","geo",
			 "vol:layer:layID");
    ntGeo->SetMarkerStyle(8);
    ntGeo->SetMarkerSize(0.3);
    initLayers();
    readLayerSizes();
    readModules();
  }
  nEvents++;

  if( endOfData ){
    statFile->Write();
    statFile->Close();
    nEvents=0;
    //writeLayerSizes();
    writeModules();
    return;
  }

  updateLayerSizesT(this);


  return;


  // analyse hit density
  /*
  for( int iv=0; iv<Geo::NVolumes; iv++){
    Volume &vol = Geo::volumes[iv];
    for( int il=0; il<vol.nLayers; il++ ){
      int layerId = vol.layerIDs[il];
      //Layer &layer = Geo::layers[layerId];
      FitLayer &flayer = mFitLayers[layerId];

      int nHits = 0;
      double sigmaXY=0, sigmaZ=0;
      int nHitsDense=0;
      double rDense = flayer.mTmin + 0.1*(flayer.mTmax - flayer.mTmin);

      for( int ih=0; ih<mHits.size(); ih++ ){
        Hit &h = mHits[ih];
	HitMC &mc = mHitsMC[ih];
	if( h.layerID!= layerId ) continue;
	if( mc.partID < 0) continue;
	double dx = h.x - mc.x;
	double dy = h.y - mc.y;
	double dz = h.z - mc.z;
	sigmaXY+=dx*dx + dy*dy;
 	sigmaZ+=dz*dz;
	nHits++;
	if( vol.type==0 ){
	  if( fabs(mc.z) <= 0.1*flayer.mTmax ){
	    nHitsDense++;
	  }
	} else {
	  double r = sqrt(mc.x*mc.x+mc.y*mc.y);
	  if( fabs(r) <= rDense ){
	    nHitsDense++;
	  }
	}
      }
      sigmaXY = 3.5*sqrt(sigmaXY/nHits/2);
      sigmaZ = 3.5*sqrt(sigmaZ/nHits);
      double S=0;
      if( vol.type==0 ){
	S = 2.*TMath::Pi()*flayer.mR*0.1*flayer.mTmax;
      } else {
	S = TMath::Pi()*(rDense*rDense - flayer.mTmin*flayer.mTmin);
      }
      double areaS = S/nHitsDense;
      double areaR = sqrt(areaS/TMath::Pi());
      cout<<"vol "<<iv<<" layer "<<il<<": dens "<< 2*areaR<<" xy=+-"<<sigmaXY<<" z=+-"<<sigmaZ<<endl;
    }
  }
  */
#endif
}
