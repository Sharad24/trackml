#ifndef CUTS_H
#define CUTS_H

#include "Geo.h"
#include <string.h>

using namespace std;


struct LayerCuts
{
  double cutPhi;
  double cutT;
  double cutV;
  double cutUZ;
  double cutVms;
  double cutUZms;

  static constexpr int nPar(){ return 6; }

  double &par( int i ){ return *(&cutPhi + i); }

  static double parInitial( int i ){
    if( i==0 ) return 0.01;
    else return 0.5;
  }

  static double parRound( int i ){
    if( i==0 ) return 1.e-5;
    else return 1.e-4;
  }

  static double parMaxStep( int i ){
    if( i==0 ) return 0.01;
    else return 0.2;
  }

};

struct RecoPassParameters
{
  int ID;
  string name;

  // fixed parameters
  int baseLayer1;
  int baseLayer2;
  int baseLayer3; // optional
  int model2; // track model at layer 2
  int model; // track model at layers 3+
  bool doBranching3; // branching at layer 3 (if baseLayer3>=0)
  bool doBranching; // branching at other layers

  int gap; // max gap between layers
  int gapTotal; // max total gap
  bool keepStopped; // keep stopped tracks

  int minTrackNLayers;
  int firstNonStopVolume;

  // adjusted parameters

  double l1AbsTMin;
  double l1AbsTMax;

  double l2PhiMin;
  double l2Phi;
  double l2T;
  double l2Pt;
  double l3Pt;

  LayerCuts searchCuts[Geo::NLayers];
  LayerCuts pickUpCuts[Geo::NLayers];

  static constexpr int nGlobalPar(){ return 7; }
  static constexpr int nPar(){ return nGlobalPar() + 2*Geo::NLayers*LayerCuts::nPar(); }

  double &par( int i ){
    constexpr int nLayerParam = LayerCuts::nPar();
    constexpr int nLayerParamTot = Geo::NLayers*nLayerParam;
    //if( i < nGlobalPar() ) return *( (&l1AbsTMin) + i);
    if( i == 0 )  return l1AbsTMin;
    if( i == 1 )  return l1AbsTMax;
    if( i == 2 )  return l2PhiMin;
    if( i == 3 )  return l2Phi;
    if( i == 4 )  return l2T;
    if( i == 5 )  return l2Pt;
    if( i == 6 )  return l3Pt;

    i-= nGlobalPar();
    if( i< nLayerParamTot ){
      int il = i/nLayerParam;
      int ip = i % nLayerParam;
      return searchCuts[il].par( ip );
    }
    i-=nLayerParamTot;
    int il = i/nLayerParam;
    int ip = i % nLayerParam;
    return pickUpCuts[il].par( ip );
  }


  static double parInitial( int i ){
    constexpr int nLayerParam = LayerCuts::nPar();
    if( i==0 ) return 0.;
    if( i==1 ) return 50.;
    if( i==2 ) return 0.;
    if( i==3 ) return 0.1;
    if( i==4 ) return 50;
    if( i==5 ) return 0.;
    if( i==6 ) return 0.;
    i-= nGlobalPar();
    int ip = i % nLayerParam;
    return LayerCuts::parInitial( ip );
  }

  static double parRound( int i ){
    constexpr int nLayerParam = LayerCuts::nPar();
    if( i==0 ) return 1.e-2;
    if( i==1 ) return 1.e-2;
    if( i==2 ) return 1.e-4;
    if( i==3 ) return 1.e-4;
    if( i==4 ) return 1.e-2;
    if( i==5 ) return 0.001;
    if( i==6 ) return 0.001;
    i-= nGlobalPar();
    int ip = i % nLayerParam;
    return LayerCuts::parRound( ip );
  }

  static double parMaxStep( int i ){
    constexpr int nLayerParam = LayerCuts::nPar();
    if( i==0 ) return 1.;
    if( i==1 ) return 1.;
    if( i==2 ) return 0.001;
    if( i==3 ) return 0.001;
    if( i==4 ) return 1.;
    if( i==5 ) return 0.1;
    if( i==6 ) return 0.1;
     i-= nGlobalPar();
    int ip = i % nLayerParam;
    return LayerCuts::parMaxStep( ip );
  }


};

istream & operator >> ( istream & in,  RecoPassParameters &cuts );
ostream & operator << ( ostream& out,  const RecoPassParameters &cuts );




struct Cuts
{
  static void init( const char *file = "/Users/sharadchitlangia/Desktop/trackml/code/MikadoTracker/cuts.txt");

  static int mirrorVolume( int v );
  static int mirrorLayer( int layerID );


  static RecoPassParameters mirror( const RecoPassParameters  &fwd, int ID ) ;


  static vector<RecoPassParameters> sliceCuts;


  static void ReadCuts( const char *file="/Users/sharadchitlangia/Desktop/trackml/code/MikadoTracker/cuts.txt" );
  static void WriteCuts( const char *file="cuts.txt" );

};


#endif
