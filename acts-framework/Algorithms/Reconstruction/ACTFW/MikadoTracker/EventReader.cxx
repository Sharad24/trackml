#include "EventReader.h"
#include "Tracker.h"
#include "util.h"
#include <iostream>
#include "TrackModelPhysical.h"
#include "Cuts.h"
#include <map>
#include <algorithm>
#include <random>


bool EventReader::checkLayers( const int *baseLayers, int nBaseLayers, Particle &p, double &fitpt, double &phi )
{
  fitpt = -1.;
  phi = -1;
  int icross=0;
  int lNHits[nBaseLayers];
  for( int i=0; i<nBaseLayers; i++) lNHits[i] = 0;
  for( uint i=0; i<p.hits.size(); i++ ){
    int l = mTracker->mHits[p.hits[i]].layerID;
    for( int j=0; j<nBaseLayers; j++ ){
      if( l==baseLayers[j] ) lNHits[j]++;
    }
  }
  for( int il=0; il<nBaseLayers; il++ ){
    if( lNHits[il]>0 ) icross++;    
  }
  
  if( icross!=nBaseLayers ) return 0;

  int layID[3];
  if( nBaseLayers==2 ){
    layID[0] = -1;
    layID[1] = baseLayers[0];
    layID[2] = baseLayers[1];
  } else if( nBaseLayers==3 ){
    layID[0] = baseLayers[0];
    layID[1] = baseLayers[1];
    layID[2] = baseLayers[2];
  } else {
    cout<<"wrong number of base layers!!"<<endl;
    exit(1);
  }  

  for( uint i0=0; i0<p.hits.size(); i0++ ){ 
    Hit h0;
    if( layID[0] < 0 ){
      h0.fitHit.x = 0.;
      h0.fitHit.y = 0.;
      h0.fitHit.z = 0.;
      h0.phi=0.;
      i0 = p.hits.size();
    } else {
      if( mTracker->mHits[p.hits[i0]].layerID != layID[0] ) continue;
      h0 = mTracker->mHits[p.hits[i0]];
    }
    
    for( uint i1=0; i1<p.hits.size(); i1++ ){
     if( mTracker->mHits[p.hits[i1]].layerID != layID[1] ) continue;
     Hit h1 = mTracker->mHits[p.hits[i1]];
           
     if( nBaseLayers==3){
       double dphi = fabs(h1.phi-h0.phi);
       if( dphi > M_PI ) dphi = fabs(2*M_PI - dphi);
       if( phi<0. || phi < dphi )  phi = dphi;
     }

     for( uint i2=0; i2<p.hits.size(); i2++ ){
       if( mTracker->mHits[p.hits[i2]].layerID != layID[2] ) continue;
       Hit h2 = mTracker->mHits[p.hits[i2]];

       if( nBaseLayers==2){
	 double dphi = fabs(h2.phi-h1.phi);	 
	 if( dphi > M_PI ) dphi = fabs(2*M_PI - dphi);
	 if( phi<0. || phi < dphi )  phi = dphi;       
       }

       TrackModelPhysical t;
       double Bz = Geo::OriginBzkG;    
       int err = t.createXYLast( h0.x(), h0.y(), h0.z(),
				 h1.x(), h1.y(), h1.z(),
				 h2.x(), h2.y(), h2.z() , Bz  );	    
       if( err!=0 ) continue;
       if( fitpt<0 || t.getPtkG() < fitpt ) fitpt = t.getPtkG();
     }
    }
  }
  
  return 1;
}



int EventReader::readEvent( const char *directory, int event, bool loadMC )
{

  char filePrefix[256];
  sprintf( filePrefix, "%sevent%09d",directory,event);

  mTracker->mcFlag = loadMC;

  mTracker->mHits.clear();
  //mHitRecInfo.clear();
  //mTracks.clear();
  mTracker->mHitsMC.clear();
  mTracker->mParticles.clear();
  mTracker->mHits.reserve(150000); 

  // ===== load hits
  { 
    char fname[256];
    sprintf( fname, "%s-hits.csv", filePrefix);

    ifstream in(fname);    
    if( !in.is_open() ){
      cout<<"Event "<<event<<" does not exist!!"<<endl;
      return -1;
      //exit(0);
    }
    char tmpLine[256];
    in.getline(tmpLine,256);
    cout<<tmpLine<<endl;
    while (1) {
      double h[7]; //file line: id:x:y:z:volume:layer:module
      if( !readLine(in,h,7) ) break;
      if( h[0]-1 != mTracker->mHits.size() ){
	cout<<"Hit index is wrong: "<<h[0]<<endl;
	exit(-1);
      }      
      Hit hit;
      hit.fitHit.x = h[1]*0.1; // convert to [cm]
      hit.fitHit.y = h[2]*0.1; // convert to [cm]
      hit.fitHit.z = h[3]*0.1; // convert to [cm]      
      hit.r = sqrt( hit.x()*hit.x() + hit.y()*hit.y() );
      hit.phi = atan2( hit.y(), hit.x() );
      //hit.module = h[6]-1;
      int layer = ( (int) h[5])/2 -1;
      int volume = 0;
      switch( (int) h[4] )
	{
	case  8: volume = 0; break;
	case  7: volume = 1; layer = 6-layer; break;
	case  9: volume = 2; break;
	case 13: volume = 3; break;
	case 12: volume = 4; layer = 5-layer; break;
	case 14: volume = 5; break;
	case 17: volume = 6; break;
	case 16: volume = 7; layer = 5-layer; break;
	case 18: volume = 8; break;
	default:	  
	  cout<<"Unknown detector volume: "<< (int) h[4] << endl;
	  exit(-1);	
	};     

      if( layer<0 || layer>= Geo::volumes[volume].nLayers ){
	cout<<"Unknown detector layer: "<<layer<<endl;
	exit(-1);	
      }      
      hit.layerID = Geo::volumes[volume].layerIDs[layer];
      Layer &glayer = Geo::layers[hit.layerID];
      if( Geo::volumes[volume].type==0 ){
	hit.t = hit.z() / hit.r * glayer.r;
      } else {
	hit.t = hit.r / hit.z() * glayer.z ;
      }
    
      hit.isUsed=0;
      hit.trackID=0;
      //hit.nCellsU=1;
      //hit.nCellsV=1;
      //hit.layerHitID = -1;
      FitHit &fh = hit.fitHit;
      fh.BzFwd = Geo::OriginBzkG;
      fh.BzMid = Geo::OriginBzkG;
      fh.BzBck = Geo::OriginBzkG;
      int vol = Geo::getVolume(hit.layerID);
      if( 1 || ( vol==4 || vol==5 || vol==7 || vol==8 ) ){
	fh.BzFwd = Geo::layers[hit.layerID].getFieldFwd(hit.phi, hit.t);
 	fh.BzMid = Geo::layers[hit.layerID].getFieldMid(hit.phi, hit.t);
	fh.BzBck = Geo::layers[hit.layerID].getFieldBck(hit.phi, hit.t);
      }

      mTracker->mHits.push_back(hit);
    }
    cout<<" loaded "<<mTracker->mHits.size()<<" hits "<<endl;
    in.close();    

  } // load hits


  // ===== load cells
  { 
    int nHits = (int) mTracker->mHits.size();
    vector<int> cMinU( nHits );
    vector<int> cMaxU( nHits );
    vector<int> cMinV( nHits );
    vector<int> cMaxV( nHits );
    for( int i=0; i<nHits; i++ ){
      cMinU[i] = 1000000;
      cMaxU[i] = -1;
      cMinV[i] = 1000000;
      cMaxV[i] = -1;
    }

    char fname[256];
    sprintf( fname, "%s-cells.csv", filePrefix);

    ifstream in(fname );    
    if( !in.is_open() ){
      cout<<"Event "<<event<<": cells file does not exist!!"<<endl;
      exit(-1);
    }
    char tmpLine[256];
    in.getline(tmpLine,256);
    cout<<tmpLine<<endl;
    int nCells=0;
    while (1) {    
      double c[4]; //file line: hit_id,ch0,ch1,value
      if( !readLine(in,c,4) ) break;       
      int ihit = (int) c[0] - 1;
      if( ihit<0 || ihit>=nHits ){
	cout<<"read cells: Hit index is wrong: "<<ihit<<endl;
	exit(-1);
      }
      int ch0 = (int) c[1];
      int ch1 = (int) c[2];
      if( cMinU[ihit] > ch0 ) cMinU[ihit] = ch0;
      if( cMaxU[ihit] < ch0 ) cMaxU[ihit] = ch0;
      if( cMinV[ihit] > ch1 ) cMinV[ihit] = ch1;
      if( cMaxV[ihit] < ch1 ) cMaxV[ihit] = ch1;
      nCells++;
    }
    cout<<" loaded "<<nCells<<" cells "<<endl;
    in.close();    
    for( int i=0; i<nHits; i++ ){
      //Hit &hit = mTracker->mHits[i];
      if( cMaxU[i]<0 || cMaxV[i]<0 ){
	cout<<"read cells: a hit has no cells, something wrong "<<endl;
	exit(-1);
      }
      //hit.nCellsU = cMaxU[i] - cMinU[i] + 1;
      //hit.nCellsV = cMaxV[i] - cMinV[i] + 1;
    }

  } // load cells


  if( mTracker->mcFlag ){
    int err = readMCEvent( directory, event );
    if( err!=0 ) return err;
  }

 
  // create layer hits
  {
    int *nLayerHits = new int[Geo::NLayers];
    for( int il=0; il<Geo::NLayers; il++ ){
      nLayerHits[il]=0;
    }
    for( uint ih=0; ih<mTracker->mHits.size(); ih++ ){
      Hit &hit = mTracker->mHits[ih];
      nLayerHits[hit.layerID]++;
    }
    for( int il=0; il<Geo::NLayers; il++ ){      
      mTracker->mLayerHits[il].clear();
      mTracker->mLayerHits[il].reserve(nLayerHits[il]);      
    }
    for( uint ih=0; ih<mTracker->mHits.size(); ih++ ){
      Hit &h = mTracker->mHits[ih];
      SearchHit fh;
      //fh.fitHit = h.fitHit;
      //fh.phi = h.phi;
      fh.phi = h.phi;
      fh.t = h.t;
      fh.id = ih;

      // SG!!! check MC info
      if( 0 && mTracker->mcFlag ){
	const HitMC &mc = mTracker->mHitsMC[ih];
	if( mc.partID<0 ) continue;    		
	//if( mc.partID!=1155 ) continue;    	
	if( 1 && mc.partID>=0 ){
	  Particle &p = mTracker->mParticles[mc.partID];
	  //cout<<" orig ID "<<p.origID<<endl;
	  if( p.phiV1 < 0.0   ) continue;	  
	  if( p.phiV1 >= 0.01   ) continue;
	  //if( p.phiV0L12 >0.1   ) continue;
	  //if( p.phiV1 >= 0.055 ) continue;
	}
      }     
      mTracker->mLayerHits[h.layerID].push_back(fh);
      //h.layerHitID = mLayerHits[h.layerID].size()-1;
    }
    delete[] nLayerHits;
  }     
 
  return 0;
}



int EventReader::readEvent( float *x, float *y, float *z, int *id, int *volumes, int *layers, int nHits )
{
  mTracker->mcFlag = 0;

  mTracker->mHits.clear();
  mTracker->mHitsMC.clear();
  mTracker->mParticles.clear();
  mTracker->mHits.reserve(nHits); 

  // ===== load hits
  
  for( int i=0; i<nHits; i++ ){     
    Hit hit;
    if( id[i]!=i+1 ){
      cout<<"unexpected hit IDs "<<endl;
      exit(-1);
    }

    hit.fitHit.x = x[i]*0.1; // convert to [cm]
    hit.fitHit.y = y[i]*0.1; // convert to [cm]
    hit.fitHit.z = z[i]*0.1; // convert to [cm]      
    hit.r = sqrt( hit.x()*hit.x() + hit.y()*hit.y() );
    hit.phi = atan2( hit.y(), hit.x() );
    int layer = ( (int) layers[i])/2 -1;
    int volume = 0;
    switch( (int) volumes[i] )
      {
      case  8: volume = 0; break;
      case  7: volume = 1; layer = 6-layer; break;
      case  9: volume = 2; break;
      case 13: volume = 3; break;
      case 12: volume = 4; layer = 5-layer; break;
      case 14: volume = 5; break;
      case 17: volume = 6; break;
      case 16: volume = 7; layer = 5-layer; break;
      case 18: volume = 8; break;
      default:	  
	cout<<"Unknown detector volume: "<< volumes[i] << endl;
	exit(-1);	
      };     

    if( layer<0 || layer>= Geo::volumes[volume].nLayers ){
      cout<<"Unknown detector layer: "<<layer<<endl;
      exit(-1);	
    }      
    hit.layerID = Geo::volumes[volume].layerIDs[layer];
    Layer &glayer = Geo::layers[hit.layerID];
    if( Geo::volumes[volume].type==0 ){
      hit.t = hit.z() / hit.r * glayer.r;
    } else {
      hit.t = hit.r / hit.z() * glayer.z ;
    }
    
    hit.isUsed=0;
    hit.trackID=0;
    //hit.nCellsU=1;
    //hit.nCellsV=1;
    //hit.layerHitID = -1;
    FitHit &fh = hit.fitHit;
    fh.BzFwd = Geo::OriginBzkG;
    fh.BzMid = Geo::OriginBzkG;
    fh.BzBck = Geo::OriginBzkG;
    int vol = Geo::getVolume(hit.layerID);
    if( 1 || ( vol==4 || vol==5 || vol==7 || vol==8 ) ){
      fh.BzFwd = Geo::layers[hit.layerID].getFieldFwd(hit.phi, hit.t);
      fh.BzMid = Geo::layers[hit.layerID].getFieldMid(hit.phi, hit.t);
      fh.BzBck = Geo::layers[hit.layerID].getFieldBck(hit.phi, hit.t);
    }

    mTracker->mHits.push_back(hit);
  }
  if( mTracker->doPrint ) cout<<" loaded "<<mTracker->mHits.size()<<" hits "<<endl;
  

  // create layer hits
  {
    int *nLayerHits = new int[Geo::NLayers];
    for( int il=0; il<Geo::NLayers; il++ ){
      nLayerHits[il]=0;
    }
    for( uint ih=0; ih<mTracker->mHits.size(); ih++ ){
      Hit &hit = mTracker->mHits[ih];
      nLayerHits[hit.layerID]++;
    }
    for( int il=0; il<Geo::NLayers; il++ ){      
      mTracker->mLayerHits[il].clear();
      mTracker->mLayerHits[il].reserve(nLayerHits[il]);
    }
    for( uint ih=0; ih<mTracker->mHits.size(); ih++ ){
      Hit &h = mTracker->mHits[ih];
      SearchHit fh;
      fh.phi = h.phi;
      fh.t = h.t;
      fh.id = ih;
      mTracker->mLayerHits[h.layerID].push_back(fh);
    }
    delete[] nLayerHits;
  }     
 
  return 0;
}






int EventReader::readMCEvent( const char *directory, int event )
{
  char filePrefix[256];
  sprintf( filePrefix, "%sevent%09d",directory,event);

  mTracker->mcFlag = 1;

  mTracker->mHitsMC.clear();
  mTracker->mParticles.clear();

  // create particle ID->index map
  std::map<long unsigned int,int> partIDmap;   

  // ========= load particles with reindexing
  {
    mTracker->mParticles.clear();          

    char fname[256];
    sprintf( fname, "%s-particles.csv", filePrefix);

    ifstream in(fname);
    if( !in.is_open() ){
      cout<<"Particle file for event "<<event<<" does not exist!!"<<endl;
      exit(0);
    }
    char tmpLine[256];
    in.getline(tmpLine,256);
    cout<<tmpLine<<endl;
    while (1) {    
      double f[10]; //  particle_id,particle_type,vx,vy,vz,px,py,pz,q,nhits
      if( !readLine(in,f,10) ) break;
      int nhits = (int) f[9];
      if( nhits==0 ) continue; // no hits belong to the particle		
      Particle p( nhits );
      p.countID = mTracker->mParticles.size();
      p.origID = ( (long unsigned int) f[0]);
      p.pid = (int) f[1];
      p.x = f[2]*0.1; // [cm]
      p.y = f[3]*0.1; // [cm]
      p.z = f[4]*0.1; // [cm]
      p.r = sqrt(p.x*p.x+p.y*p.y);
      p.px = f[5];
      p.py = f[6];
      p.pz = f[7];
      p.q  = f[8];
      p.xl = 0;
      p.yl = 0;
      p.zl = 0;
      p.rl = 0;
      p.pt = sqrt(p.px*p.px + p.py*p.py);
      p.p  = sqrt(p.px*p.px + p.py*p.py + p.pz*p.pz );
      //p.prim = fabs(p.z)<1.2 && p.r<0.05;
      p.prim = fabs(p.z)<=21. && p.r<=0.006;      
      p.w = 0;
      p.isSelected = 1;
      p.dcaR = fabs( -p.x*p.py + p.y*p.px )/p.pt;
      p.dcaZ = p.z - p.pz*( p.x*p.px + p.y*p.py )/p.pt/p.pt;

      partIDmap[ p.origID ] = mTracker->mParticles.size();
      //if( mParticles.size()==293 ) cout<<"origID = "<<p.origID<<endl;
      mTracker->mParticles.push_back(p);      
    }
    cout <<" loaded "<<mTracker->mParticles.size() <<" particles in event "<<event<<endl;    
    in.close();
  } // particles



  { // ============  read  mc truth

    mTracker->mHitsMC.clear();
    mTracker->mHitsMC.reserve(150000);

    char fname[256];
    sprintf( fname, "%s-truth.csv", filePrefix);

    ifstream in(fname);    
    if( !in.is_open() ){
      cout<<"Truth file for event "<<event<<" does not exist!!"<<endl;
      exit(0);
    }

    std::default_random_engine generator;
    generator.seed(1);
    std::uniform_real_distribution<double> distribution(0.,1.);

    char tmpLine[256];    
    in.getline(tmpLine,256);
    cout<<tmpLine<<endl;

    while (1) {    
      double mc[9]; //mc: hit_id,particle_id,tx,ty,tz,tpx,tpy,tpz,weight     
      if( !readLine(in,mc,9) ) break;
      if( mc[0]-1 != mTracker->mHitsMC.size() ){
	cout<<"MC hit index is wrong: "<<mc[0]<<endl;
	exit(0);
      } 
      HitMC hitmc;
      Hit &hit = mTracker->mHits[mTracker->mHitsMC.size()];
 
      hitmc.hitID = mTracker->mHitsMC.size();
      hitmc.x = mc[2]*0.1; // convert to [cm]
      hitmc.y = mc[3]*0.1; // convert to [cm]
      hitmc.z = mc[4]*0.1; // convert to [cm]      
      hitmc.partID = -1;
      hitmc.w = mc[8]; // weight 
      hitmc.px = mc[5];
      hitmc.py = mc[6];
      hitmc.pz = mc[7];
      hitmc.pt = sqrt(hitmc.px*hitmc.px + hitmc.py*hitmc.py); // pt
      hitmc.dOrigin2 = -1.; // squared distance to particle origin 
      hitmc.p = sqrt( hitmc.px*hitmc.px + hitmc.py*hitmc.py + hitmc.pz*hitmc.pz ); // momentum  for sorting
      hitmc.q = -1000;
      hitmc.r = sqrt( hitmc.x*hitmc.x + hitmc.y*hitmc.y );
      hitmc.phi = atan2( hitmc.y, hitmc.x );     
      if( Geo::volumes[Geo::getVolume(hit.layerID)].type==0 ) hitmc.t = hitmc.z;
      else hitmc.t = hitmc.r;

      hitmc.random = distribution(generator);//( (double) rand() ) / RAND_MAX;

      if(0){//SG!!! put MC truth in hit:
	hit.fitHit.x = hitmc.x;
	hit.fitHit.y = hitmc.y;
	hit.fitHit.z = hitmc.z;
	hit.r = sqrt( hit.x()*hit.x() + hit.y()*hit.y() );
	hit.phi = atan2( hit.y(), hit.x() );
 	Layer &layer = Geo::layers[hit.layerID];
	if( Geo::volumes[Geo::getVolume(hit.layerID)].type==0 ){
	  hit.t = hit.z() / hit.r * layer.r;
	} else {
	  hit.t = hit.r / hit.z() * layer.z ;
	}
      }

      // find mapped particle id
      long unsigned int id = mc[1];
      int newID = 0;
      if( id==0 ){ // hit is not associated to any particle 
	mTracker->mHitsMC.push_back(hitmc);
	continue;     	
      }
      std::map<long unsigned int, int>::iterator it = partIDmap.find(id);
      if( it==partIDmap.end() || it->first!=id){
	cout<<"Mapped particle ID is wrong!!!"<<endl;
	cout<<"ID= "<<id<<" hit "<<id<<" iterator at ID "<<it->first<<endl;
	exit(0);
      }
      newID = it->second;
      if( newID < 0 || newID >= (int) mTracker->mParticles.size() ){
	cout<<"Mapped particle ID is wrong!!!"<<endl;
	cout<<"ID= "<<id<<" new ID "<<newID<<endl;
	exit(0);
      }
      hitmc.partID = newID;     
      Particle &p = mTracker->mParticles[newID];

      double dx = hitmc.x - p.x;
      double dy = hitmc.y - p.y;
      double dz = hitmc.z - p.z;
      hitmc.dOrigin2 = dx*dx + dy*dy + dz*dz;

      p.w += hitmc.w;
      hitmc.q = p.q;

      p.hits.push_back( hitmc.hitID );
      mTracker->mHitsMC.push_back(hitmc);           
    }
    cout << " read "<<mTracker->mHitsMC.size() << " mc hits for event "<<event<<endl;
    in.close();      
    if( mTracker->mHitsMC.size() != mTracker->mHits.size() ){
      cout<<"number of MC hits is not equal to the number of hits"<<endl;
      exit(0);
    }
  } // read mc hits
    
    
  // sort particle hits

  for( unsigned int ipart=0; ipart<mTracker->mParticles.size(); ipart++ ){      
    Particle &p = mTracker->mParticles[ipart];
    if( p.hits.size()<1 ) continue;
    std::vector<HitMC> v;
    for( uint i=0; i<p.hits.size(); i++ ){
      v.push_back( mTracker->mHitsMC[p.hits[i]] );
    }
    std::sort(v.begin(), v.end() );
    for( uint i=0; i<p.hits.size(); i++ ){
      p.hits[i] = v[i].hitID;
    }
    p.hitClusterIds.resize(p.hits.size());
    int il=-1;
    int lastLid=-1;
    for( uint i=0; i<p.hits.size(); i++ ){
      if( mTracker->mHits[p.hits[i]].layerID != lastLid ){
	il++;
	lastLid = mTracker->mHits[p.hits[i]].layerID;
      }
      p.hitClusterIds[i] = il;
    } 

    {
      int ihlast = p.hits[p.hits.size()-1];
      HitMC &mc = mTracker->mHitsMC[ihlast];      
      p.xl = mc.x;
      p.yl = mc.y;
      p.zl = mc.z;
      p.rl = sqrt(mc.x*mc.x+mc.y*mc.y);
    }
    if( 0 && ipart==467 ){
      for( uint i=0; i<p.hits.size(); i++ ){
	Hit &h = mTracker->mHits[p.hits[i]];      
	h.Print();
      }
      for( uint i=0; i<p.hits.size(); i++ ){
	HitMC &mc = mTracker->mHitsMC[p.hits[i]];      
	mc.Print();	
      }
    }

    // exclude some particles    
    if( 0 && 
       ( fabs(p.z)>16.5 || p.r>0.2 || p.hits.size()<8 ) ){
      p.w = 0;
      for( uint i=0; i<p.hits.size(); i++ ){
	HitMC &mc = mTracker->mHitsMC[p.hits[i]];      
	mc.w = 0;
      }
    }
  }
  //exit(0);
  

  for( unsigned int ipart=0; ipart<mTracker->mParticles.size(); ipart++ ){      
    Particle &p = mTracker->mParticles[ipart];
    p.isSelected = ( p.w > 0.f );
    for( int i=0; i<Geo::NLayers; i++ ) layerNHits[i]=0;

    for( uint i=0; i<p.hits.size(); i++ ){
      layerNHits[mTracker->mHits[p.hits[i]].layerID]++;      
    }
 
   p.nLayers = 0;
    for( int i=0; i<Geo::NLayers; i++ ){
      if( layerNHits[i]>0 ) p.nLayers ++;
    }

    p.baseV0 = 0;
    p.baseV1 = 0;
    p.baseV2 = 0;
    p.baseV  = 0;
    p.ptV = -1;
    p.ptV0 = -1;
    p.ptV1 = -1;
    p.ptV2 = -1;

    if( 1 || p.prim ){ //SG!!
      int baseV0[2] ={ Geo::getLayerID(0,0), Geo::getLayerID(0,1) };
      int baseV1[2] ={ Geo::getLayerID(1,0), Geo::getLayerID(1,1) };
      int baseV2[2] ={ Geo::getLayerID(2,0), Geo::getLayerID(2,1) };
      p.baseV0 = checkLayers(baseV0,2, p, p.ptV0, p.phiV0 );
      p.baseV1 = checkLayers(baseV1,2, p, p.ptV1, p.phiV1 );
      p.baseV2 = checkLayers(baseV2,2, p, p.ptV2, p.phiV2 );
      p.baseV  =  ( p.baseV0 || p.baseV1 || p.baseV2 );      
      p.ptV = 1.e20;
      if( p.ptV0>=0 && p.ptV > p.ptV0 ) p.ptV = p.ptV0;
      if( p.ptV1>=0 && p.ptV > p.ptV1 ) p.ptV = p.ptV1;
      if( p.ptV2>=0 && p.ptV > p.ptV2 ) p.ptV = p.ptV2;
      if( p.ptV > 1.e19 ) p.ptV = -1;

      p.phiV = -1;
      if( p.phiV < p.phiV0 ) p.phiV = p.phiV0;
      if( p.phiV < p.phiV1 ) p.phiV = p.phiV1;
      if( p.phiV < p.phiV2 ) p.phiV = p.phiV2;
    }
    
    int baseA0[3] ={ Geo::getLayerID(0,0), Geo::getLayerID(0,1), Geo::getLayerID(0,2) };
    p.baseA0 = checkLayers(baseA0,3,p, p.ptA0, p.phiA0 );

    int baseA1[3] ={ Geo::getLayerID(1,0), Geo::getLayerID(1,1), Geo::getLayerID(1,2) };
    p.baseA1 = checkLayers(baseA1,3,p, p.ptA1, p.phiA1 );

    int baseA2[3] ={ Geo::getLayerID(2,2), Geo::getLayerID(2,3), Geo::getLayerID(2,4) };
    p.baseA2 = checkLayers(baseA2,3,p, p.ptA2, p.phiA2);

    int baseA3[3] ={ Geo::getLayerID(3,0), Geo::getLayerID(3,1), Geo::getLayerID(3,2) };
    p.baseA3 = checkLayers(baseA3,3,p, p.ptA3, p.phiA3);
    
    double pt;
    
    int baseV0L12[2] ={ Geo::getLayerID(0,1), Geo::getLayerID(0,2) };
    checkLayers( baseV0L12, 2, p, pt, p.phiV0L12 );

    int baseV0L23[2] ={ Geo::getLayerID(0,2), Geo::getLayerID(0,3) };
    checkLayers( baseV0L23, 2, p, pt, p.phiV0L23 );

    int baseV0V3[2] ={ Geo::getLayerID(0,3), Geo::getLayerID(3,0) };
    checkLayers( baseV0V3, 2, p, pt, p.phiV0V3 );


    int baseV0V1[2] ={ Geo::getLayerID(0,0), Geo::getLayerID(1,0) };
    checkLayers(baseV0V1,2,p, pt, p.phiV0V1 );
 
    int baseV3L01[2] ={ Geo::getLayerID(3,0), Geo::getLayerID(3,1) };
    checkLayers(baseV3L01,2,p, pt, p.phiV3L01 );
 
    int baseV3L12[2] ={ Geo::getLayerID(3,1), Geo::getLayerID(3,2) };
    checkLayers(baseV3L12,2,p, pt, p.phiV3L12 );

    int baseV1L23[2] ={ Geo::getLayerID(1,2), Geo::getLayerID(1,3) };
    checkLayers(baseV1L23,2,p, pt, p.phiV1L23 );

    int baseV2L23[2] ={ Geo::getLayerID(2,2), Geo::getLayerID(2,3) };
    checkLayers(baseV2L23,2,p, pt, p.phiV2L23 );

    int baseV1L45[2] ={ Geo::getLayerID(1,4), Geo::getLayerID(1,5) };
    checkLayers( baseV1L45, 2, p, pt, p.phiV1L45 );

    int baseV2L45[2] ={ Geo::getLayerID(2,4), Geo::getLayerID(2,5) };
    checkLayers( baseV2L45, 2, p, pt, p.phiV2L45 );

    int baseV1L56[2] ={ Geo::getLayerID(1,5), Geo::getLayerID(1,6) };
    checkLayers( baseV1L56, 2, p, pt, p.phiV1L56 );

    int baseV2L56[2] ={ Geo::getLayerID(2,5), Geo::getLayerID(2,6) };
    checkLayers( baseV2L56, 2, p, pt, p.phiV2L56 );
 
    if( p.phiV0 >= 0 && ipart%20 !=0 ){
      p.phiV0L12 = -1;
      p.phiV0L23 = -1;
      p.phiV0V3 = -1;
      p.phiV3L01 = -1;
      p.phiV3L12 = -1;
    }

  } // particles

  return 0;
}


