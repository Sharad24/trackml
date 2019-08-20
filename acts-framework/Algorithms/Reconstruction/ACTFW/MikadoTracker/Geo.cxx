#include "Geo.h"
#include <iostream>
#include <iomanip>
#include "Cuts.h"

 constexpr int Geo::NVolumes;
 constexpr int Geo::NLayers; // 48

 Volume Geo::volumes[NVolumes];
 Layer  Geo::layers[NLayers];



void Volume::init( int s, int t, int nL, int &nLayersTot)
{
  size = s;
  type = t;
  nLayers = nL;
  for( int i=0; i<nLayers; i++ ) layerIDs[i] = nLayersTot++;
}


void Geo::init( const char *path )
{
  cout<<"Geo:: init..."<<endl;
  // Volume:
  //  size; // radial position: 0 = inner, 1 = middle, 2 = outer
  //  type; // 0 = radial modules, 1 = vertical modules
  //  nLayers; // number of layers [2..7]

  for( int i=0; i<NLayers; i++ ){
    Layer &l = layers[i];
    l.id = i; // general layer id
  }

  int nLayersTot = 0;
  volumes[0].init( 0,  0, 4, nLayersTot);
  volumes[1].init( 0, -1, 7, nLayersTot);
  volumes[2].init( 0,  1, 7, nLayersTot);

  volumes[3].init( 1,  0, 4, nLayersTot);
  volumes[4].init( 1, -1, 6, nLayersTot);
  volumes[5].init( 1,  1, 6, nLayersTot);

  volumes[6].init( 2,  0, 2, nLayersTot);
  volumes[7].init( 2, -1, 6, nLayersTot);
  volumes[8].init( 2,  1, 6, nLayersTot);

  if( nLayersTot != NLayers ){
    cout<<"wrong n layers!!"<<endl;
    exit(0);
  }

  for( int iv=0; iv<NVolumes; iv++ ){
    Volume &v = volumes[iv];
    for( int il=0; il<v.nLayers; il++ ){
      Layer &l = layers[v.layerIDs[il]];
      l.volume = iv;
      l.layer = il;
      l.size = v.size;
      l.type = v.type;
    }
  }

  char filePath[256];


  // read layer modules
  {
    // sprintf( filePath, "%s/%s",path,"geoLayerModules.txt");
    sprintf( filePath, "%s/%s","/Users/sharadchitlangia/Desktop/trackml/code/MikadoTracker","geoLayerModules.txt");
    ifstream in(filePath);
    if( !in.is_open() ){
      cout<<"Geo:: Can not open geoLayerModules.txt file"<<endl;
      exit(0);
    }
    char tmpLine[256];
    in.getline(tmpLine,256);
    //cout<<tmpLine<<endl;
    while (1) {
      // layerID,volume,layer,module,alpha,pitchU,pitchV,thick05
      int ilayer;
      int tmp;
      double alpha, pitchU, pitchV, thick05;
      in >> ilayer >> tmp >> tmp >> tmp >> alpha >> pitchU >> pitchV>>thick05;
      if( !in.good() ) break;
      if( ilayer<0 || ilayer>=NLayers ){
	cout<<"wrong content of geoLayerModules.txt file"<<endl;
	cout<<"layerID "<<ilayer<<endl;
	exit(-1);
      }
      constexpr double PI = 3.14159265359;
      alpha = alpha*PI/180.;

      Layer &l = layers[ilayer];
      l.pitchU = pitchU;
      l.pitchV = pitchV;
      l.thick05 = thick05;
      l.sigmaUZ2 = l.pitchV*l.pitchV/12.;

      Module m;
      m.alpha = alpha;
      m.sinA = sin(m.alpha);
      m.cosA = cos(m.alpha);
      l.modules.push_back(m);
      l.nModules=l.modules.size();
    }
    in.close();
  }


  // read layer sizes

  ReadLayerGeometry(path);

   // read layer field
  {
    static constexpr double CLight = 0.000299792458; // speed of light
    // sprintf( filePath, "%s/%s",path,"geoLayerField.txt");
    sprintf( filePath, "%s/%s","/Users/sharadchitlangia/Desktop/trackml/code/MikadoTracker","geoLayerField.txt");
    ifstream in(filePath);
    if( !in.is_open() ){
      cout<<"Geo:: Can not open geoLayerField.txt file"<<endl;
     exit(0);
    }

    for( int il=0; il<NLayers; il++){
      Layer &l = layers[il];
      double sigma[3]={0,0,0};
      for( int ir=0; ir<3; ir++){
	int jl;
	in>>jl;
	if( jl!=il ){
	  cout<<"Geo:: geo field file broken"<<endl;
	  exit(1);
	}
	for( int i=0; i<Layer::NFieldPar; i++ ){
	  in>>l.field[ir][i];
	  l.field[ir][i]*=CLight;
	}
	in>>sigma[ir];
      }
      for( int ir=0; ir<2; ir++){
	if( sigma[ir]<=1.e-4 ){
	  for( int i=0; i<Layer::NFieldPar; i++ ) l.field[ir][i]=l.field[2][i];
	}
      }
    }
    in.close();
  }


  for( int i=0; i<NLayers; i++ ){
    //Layer &l = layers[i];
    //cout<<"layer id "<<i<<" vol "<<l.volume<<" il "<<l.layer<<endl;
  }

}

void Geo::ReadLayerGeometry( const char *path )
{
  char filePath[256];
  // sprintf( filePath, "%s/%s",path,"geoLayerSizes.txt");
  sprintf( filePath, "%s/%s","/Users/sharadchitlangia/Desktop/trackml/code/MikadoTracker","geoLayerSizes.txt");
  ifstream in(filePath);
  if( !in.is_open() ){
    cout<<"Geo:: Can not open geoLayerSizes.txt file"<<endl;
    exit(0);
  }
  char tmpLine[256];
  in.getline(tmpLine,256);

  for( int i=0; i<NLayers; i++){
    Layer &l = layers[i];
    l.ReadGeometry( in );
    if( l.id!=i ){
      cout<<"Geo: geoLayerSizes.txt file broken"<<endl;
      exit(1);
    }
  }
  in.close();

  // set neighbours for +Z side

  for( int id=0; id<NLayers; id++){
    Layer &l = layers[id];
    if( l.type == 1 ){ // +Z
      Layer &lm = layers[ Cuts::mirrorLayer(id) ];
      for( int k=0; k<5; k++ ){
	l.nextLayers[0][k].id = -1;
	l.nextLayers[1][k] =  lm.nextLayers[0][k];
	l.nextLayers[1][k].id = Cuts::mirrorLayer( lm.nextLayers[0][k].id );
	l.prevLayers[0][k].id = -1;
	l.prevLayers[1][k] =  lm.prevLayers[0][k];
	l.prevLayers[1][k].id = Cuts::mirrorLayer( lm.prevLayers[0][k].id );
      }
    } else if ( l.type == 0 ){ // radial layer
      for( int k=0; k<5; k++ ){
	l.nextLayers[1][k] = l.nextLayers[0][k];
 	l.nextLayers[1][k].id = Cuts::mirrorLayer( l.nextLayers[0][k].id );
	l.prevLayers[1][k] = l.prevLayers[0][k];
 	l.prevLayers[1][k].id = Cuts::mirrorLayer( l.prevLayers[0][k].id );
      }
    } else { // -Z side
      for( int k=0; k<5; k++ ){
	l.nextLayers[1][k].id = -1;
	l.prevLayers[1][k].id = -1;
      }
    }
  }

}

void Geo::WriteLayerGeometry( const char *path )
{
  char filePath[256];
  sprintf( filePath, "%s/%s",path,"geoLayerSizes.txt");
  ofstream out(filePath);
  if( !out.is_open() ){
    cout<<"Geo:: Can not open geoLayerSizes.txt file for writing"<<endl;
    exit(0);
  }

  out<<"id v l r z t T neigh: v l t T,  v l t T, v l t T" <<endl;

  for( int i=0; i<NLayers; i++){
    Layer &l = layers[i];
    if( l.id!=i ){
      cout<<"Geo: geometry broken"<<endl;
      exit(1);
    }
    l.WriteGeometry( out );
  }
  out.close();
}



void Layer::ReadGeometry( istream & in )
{
  int v=-1, l = -1;
  in >> id >> v >> l >> r >> z >> tMin >> tMax;
  if( id != Geo::getLayerID(v,l) ){
    cout<<"Geo: geometry broken"<<endl;
    exit(1);
  }
  for( int k=0; k<3; k++ ){
    LayerNeighbour &nb= nextLayers[0][k];
    nb.tMin = nb.tMax = 0.;
    nb.id = -1;
    v = -1;
    l = -1;
    in >> v >> l;// >> nb.tMin >> nb.tMax;
    if( v>=0 && l>=0 ) nb.id = Geo::getLayerID(v,l);
  }
  for( int k=3; k<5; k++ ){
    LayerNeighbour &nb= nextLayers[0][k];
    nb.tMin = nb.tMax = 0.;
    nb.id = -1;
  }
  for( int k=0; k<5; k++ ){
    LayerNeighbour &nb= prevLayers[0][k];
    nb.tMin = nb.tMax = 0.;
    nb.id = -1;
    v = -1;
    l = -1;
    in >> v >> l;// >> nb.tMin >> nb.tMax;
     if( v>=0 && l>=0 ) nb.id = Geo::getLayerID(v,l);
  }
}



void Layer::WriteGeometry( ostream & out )
{
  out
    << std::fixed
    << " "
    << " " << std::setw(2) << id
    << " " << std::setw(1) << Geo::getVolume(id)
    << " " << std::setw(1) << Geo::getVolumeLayer(id)

    << std::setprecision(4)
    << " " << std::setw(10) << r
    << " " << std::setw(10) << z
    << " " << std::setw(10) << tMin
    << " " << std::setw(10) << tMax
    ;
  for( int k=0; k<3; k++ ){
    LayerNeighbour &nb= nextLayers[0][k];
    int v=-1, l = -1;
    if( nb.id >= 0 ){
      v = Geo::getVolume(nb.id);
      l = Geo::getVolumeLayer(nb.id);
    }
    out
      << " " << std::setw(2) << v
      << " " << std::setw(2) << l
      //<< std::setprecision(4)
      //<< " " << std::setw(8) << nb.tMin
      //<< " " << std::setw(8) << nb.tMax
      ;
  }
  for( int k=0; k<5; k++ ){
    LayerNeighbour &nb= prevLayers[0][k];
    int v=-1, l = -1;
    if( nb.id >= 0 ){
      v = Geo::getVolume(nb.id);
      l = Geo::getVolumeLayer(nb.id);
    }
    out
      << " " << std::setw(2) << v
      << " " << std::setw(2) << l
      //<< std::setprecision(4)
      //<< " " << std::setw(8) << nb.tMin
      //<< " " << std::setw(8) << nb.tMax
      ;
  }
  out<<std::endl;
}
