#include <stdio.h>
#include "Tracker.h"
#include "EventReader.h"
#include <iostream>
#include <thread>

using namespace std;

Tracker tracker;


extern "C" void SGTrackerInit( const char *path )
{
  cout<<"init tracker.. path = "<<path<<endl;


  Geo::init( path );  

  char filePath[256];
  sprintf( filePath, "%s/%s",path,"cuts.txt");
  Cuts::ReadCuts(filePath);

  tracker.mNThreads = std::thread::hardware_concurrency();  
  cout<<"System N threads = " <<tracker.mNThreads<<endl;
  if( tracker.mNThreads < 2 ) tracker.mNThreads = 2;  
  tracker.mNThreads = 2;

  cout<<"init tracker ok"<<endl;
  tracker.doPrint = 1;
}


extern "C" void SGTrackerProcessEvent( float *x, float *y, float *z, int *id, int *vol, int *layer, int *labels, int nHits)
{
  static int iEvent=-1;
  //cout<<"tracker: processing event "<< ++iEvent <<endl;
  EventReader reader( &tracker );

  int err = reader.readEvent( x, y, z, id, vol, layer, nHits );
  if( err!=0 ){
    cout<<"Can not read event "<<iEvent<<endl;
    exit(-1);
  }
 
  tracker.reconstruct( 0 ); // don't stop at the learning point   

  //cout<<"Event "<<iEvent<<": reconstructed with "<<tracker.proc.vTracks.size()<<" tracks"<<endl;
      
  for( uint ih=0; ih<tracker.mHits.size(); ih++ ){
    Hit &h = tracker.mHits[ih];
    h.trackID=0;
  }
  
  int currentID=1;
      
  for( uint itr=0; itr<tracker.mTracks.size(); itr++ ){
    Track &t = tracker.mTracks[itr];
    if( t.nHits <=0 ) continue;
    for( int ih=0; ih<t.nHits; ih++ ){
      Hit &h = tracker.mHits[t.hitIDs[ih]];
      if( h.trackID>0 ){
	cout<<"reconstruction.cxx: Wrong hit indexing!"<<endl;
	cout<<"track "<<itr+1<<" nhits "<<t.nHits<<endl;
	cout<<"hit "<<ih<<" trackID "<<h.trackID<<endl;
	h.Print();
	exit(1);
      }
      h.trackID=currentID;
    }
    currentID++;
  }
      
  for( uint ih=0; ih<tracker.mHits.size(); ih++ ){
    Hit &h = tracker.mHits[ih];
    labels[ih] = h.trackID;
  }  
}

