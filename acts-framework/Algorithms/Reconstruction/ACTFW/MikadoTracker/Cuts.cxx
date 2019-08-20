#include "Cuts.h"

#include <sstream>
#include <iostream>
#include <iomanip>

vector<RecoPassParameters> Cuts::sliceCuts;




istream & operator>>( istream & in, RecoPassParameters &cuts)
{
  char tmpLine[256];
  in >> cuts.ID;
  in.ignore(2);
  in.getline(tmpLine,256);
  cuts.name = tmpLine;
  //cout<<"name:"<<cuts.name.c_str()<<endl;
  in.getline(tmpLine,256);
  in >> cuts.baseLayer1 >> cuts.baseLayer2 >> cuts.baseLayer3
     >> cuts.model2 >> cuts.model
     >> cuts.doBranching3 >> cuts.doBranching
     >> cuts.gap >> cuts.gapTotal >> cuts.keepStopped
     >> cuts.minTrackNLayers  >> cuts.firstNonStopVolume
     >> cuts.l1AbsTMin >> cuts.l1AbsTMax
     >> cuts.l2PhiMin >> cuts.l2Phi >> cuts.l2T >> cuts.l2Pt
     >> cuts.l3Pt;

  for( int ic=0; ic<2; ic++ ){
    LayerCuts* lCuts = (ic==0) ?cuts.searchCuts : cuts.pickUpCuts;
    in.getline(tmpLine,256);
    in.getline(tmpLine,256);
    for( int il=0; il<Geo::NLayers; il++ ){
      const Layer &l = Geo::layers[il];
      LayerCuts &c = lCuts[il];
      int id = -1, vol=-1, lay=-1;
      in >> id >> vol >> lay;
      if( id!=il || vol!=l.volume || lay!=l.layer ){
	cout<<"Cuts:: ReadCuts:: wrong indexing 2!! : "<<il<<" -> "<<id<<" "<<vol<<" "<<lay<<endl;
	exit(-1);
      }
      in >> c.cutPhi >> c.cutT >> c.cutV >> c.cutUZ >> c.cutVms >> c.cutUZms;
      if( ic==1 ){
	//c.cutPhi*=5;
	//c.cutT*=5;
      }
    }
  }

  return in;
}

ostream & operator<<( ostream & out, const RecoPassParameters &cuts)
{
  out<< endl<<cuts.ID<<"  "<<cuts.name.c_str()<<endl;
  out<<"  L1 L2 L3 mod2 mod br3 br gap totGap stopped nL stopV      t1      T1    phi2    Phi2      T2    Pt2    Pt3"<<endl;
  out
    << std::fixed
    << " "
    << " " << std::setw(2) << cuts.baseLayer1
    << " " << std::setw(2) << cuts.baseLayer2
    << " " << std::setw(2) << cuts.baseLayer3
    << " " << std::setw(4) << cuts.model2
    << " " << std::setw(3) << cuts.model
    << " " << std::setw(3) << cuts.doBranching3
    << " " << std::setw(2) << cuts.doBranching
    << " " << std::setw(3) << cuts.gap
    << " " << std::setw(6) << cuts.gapTotal
    << " " << std::setw(7) << cuts.keepStopped
    << " " << std::setw(2) << cuts.minTrackNLayers
    << " " << std::setw(5) << cuts.firstNonStopVolume

    << std::setprecision(2)
    << " " << std::setw(7) << cuts.l1AbsTMin
    << " " << std::setw(7) << cuts.l1AbsTMax
    << std::setprecision(4)
    <<"  " << std::setw(6) << cuts.l2PhiMin
    <<"  " << std::setw(6) << cuts.l2Phi
    << std::setprecision(2)
    <<"  " << std::setw(6) << cuts.l2T
    << std::setprecision(3)
    <<"  " << std::setw(5) << cuts.l2Pt
    << std::setprecision(3)
    <<"  " << std::setw(5) << cuts.l3Pt
    <<endl;
  for( int ic=0; ic<2; ic++ ){
    const LayerCuts* lCuts = (ic==0) ?cuts.searchCuts : cuts.pickUpCuts;
    out<<" id V L        phi        T        V       UZ      Vms     UZms"<<endl;
    for( int il=0; il<Geo::NLayers; il++ ){
      const Layer &l = Geo::layers[il];
      const LayerCuts &c = lCuts[il];
      out<<std::fixed<<std::setprecision(6);
      out<< std::setw(3)<<il<<" "<<l.volume<<" "<<l.layer
	 << std::setprecision(6)
	 <<"  "<< std::setw(9)<<c.cutPhi
	 << std::setprecision(4)
	 <<"  "<< std::setw(7)<<c.cutT
	 <<"  "<< std::setw(7)<<c.cutV
	 <<"  "<< std::setw(7)<<c.cutUZ
	 <<"  "<< std::setw(7)<<c.cutVms
	 <<"  "<< std::setw(7)<<c.cutUZms
	 <<endl;
    }
  }
  return out;
}


void Cuts::ReadCuts( const char *file ){

  sliceCuts.clear();

  // char filePath[256];
  // // sprintf( filePath, "%s/%s",path,"geoLayerSizes.txt");
  // sprintf( filePath, "%s/%s","/Users/sharadchitlangia/Desktop/trackml/code/MikadoTracker","geoLayerSizes.txt");

  ifstream in(file);
  if( !in.is_open() ){
    cout<<"Cuts:: Can not open input file"<<endl;
    exit(-1);
  }

  int nCuts = -1;
  in >> nCuts;

  for( int icut=0; icut<nCuts; icut++){
    RecoPassParameters cuts;
    in >> cuts;
    if( cuts.ID!= icut ){
      cout<<"Cuts:: ReadCuts:: wrong indexing!! index = "<<icut<<" id = "<<cuts.ID<<endl;
      exit(-1);
    }
    sliceCuts.push_back(cuts);
  }
  in.close();


  //sliceCuts[9] = mirror( sliceCuts[10], 9);
  //sliceCuts.push_back( mirror( sliceCuts[5], 10 ));
  // sliceCuts.push_back( convert( Cuts::CrCutsPrimV1Pt1, Cuts::primV1Pt1() , 3 ));
  //sliceCuts.push_back( convert( Cuts::CrCutsPrimV2Pt15, Cuts::primV2Pt15() , 11 ));
}


void Cuts::WriteCuts( const char *file )
{
  ofstream out(file);
  if( !out.is_open() ){
    cout<<"Cuts:: Can not open output file"<<endl;
    exit(-1);
  }

  out<<sliceCuts.size()<<endl;

  for( int icut=0; icut<(int)sliceCuts.size(); icut++){
    RecoPassParameters &cuts = sliceCuts[icut];
    if( cuts.ID != icut ){
      cout<<"Cuts:: WriteCuts:: wrong indexing!! index = "<<icut<<" id = "<<cuts.ID<<endl;
      exit(-1);
    }
    out<< cuts;
  }
  out.close();
}





int Cuts::mirrorVolume( int v )
{
  switch( v ){
  case 1:
  case 4:
  case 7:
    return v + 1;
  case 2:
  case 5:
  case 8:
    return v-1;
  case 0:
  case 3:
  case 6:
    return v;
  default:
    cout<<"cuts mirroring: something wrong!!"<<endl;
    exit(1);
  }
}

int Cuts::mirrorLayer( int layerID )
{
  if( layerID < 0 ) return layerID;
  int v = Geo::getVolume( layerID );
  int l = Geo::getVolumeLayer( layerID );
  return Geo::getLayerID( mirrorVolume(v), l);
}


RecoPassParameters Cuts::mirror( const RecoPassParameters  &fwd, int ID )
{
  RecoPassParameters bck = fwd;
  bck.ID = ID;
  bck.baseLayer1 = mirrorLayer( fwd.baseLayer1 );
  bck.baseLayer2 = mirrorLayer( fwd.baseLayer2 );
  bck.baseLayer3 = mirrorLayer( fwd.baseLayer3 );

  for( int i=0; i< Geo::NLayers; i++){
    int j = mirrorLayer( i );
    bck.searchCuts[i] = fwd.searchCuts[j];
    bck.pickUpCuts[i] = fwd.pickUpCuts[j];
  }
  return bck;
}





void Cuts::init( const char *file)
{
  ReadCuts( file );
}
