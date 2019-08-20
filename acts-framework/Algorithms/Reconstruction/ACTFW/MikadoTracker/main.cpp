// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "main.hpp"

#include <stdexcept>

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/GeometryID.hpp>

#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

#include "Tracker.h"
#include "EventReader.h"

#ifdef USEROOT
#endif

#include <iostream>
#include <thread>

using namespace std;

FW::MikadoTracker::MikadoTracker(
    const FW::MikadoTracker::Config& cfg,
    Acts::Logging::Level                            logLevel)
  : BareAlgorithm("MikadoTracker", logLevel), m_cfg(cfg)
{
  if (m_cfg.spacePointCollection.empty()) {
    throw std::invalid_argument("Missing input space points collection");
  }
}

FW::ProcessCode
FW::MikadoTracker::execute(FW::AlgorithmContext ctx) const
{
  ACTS_INFO("empty reconstruction on event " << ctx.eventNumber);
  bool analyseTruth = true;

  const char *dir = "/Users/sharadchitlangia/Desktop/trackml/train_sample/";
  //const char *dir = "../data/sample_1000/";
  int nEvents = 80;
  int firstEvent=21100;

  // big test event

  if(1){
    nEvents = 1;
    firstEvent=21139;
    //firstEvent=21147;
  }

  //TString dir = "data/train_100_events/";

  /*
  const int nEvents = 125;
  const int firstEvent=0;
  TString dir = "data/test/";
  */

  /*
  const int nEvents = 1000;
  const int firstEvent=1000;
  TString dir = "../task/train_1/";
  */

  cout<<"init tracker.."<<endl;
  Geo::init();
  Cuts::init();
  //Geo::WriteLayerGeometry();

  Tracker tracker;

  tracker.mRecoPassForLearning = -1; //18;

  cout<<"init tracker ok"<<endl;

  tracker.mNThreads = std::thread::hardware_concurrency();
  //if( tracker.mNThreads < 2 ) tracker.mNThreads = 2;

  //tracker.mNThreads = 1;

  //return 0;
  //tracker.doPrint = 0;


  Cuts::ReadCuts("/Users/sharadchitlangia/Desktop/trackml/code/MikadoTracker/cuts.txt");
  //Cuts::ReadCuts("doc/cutsEv21139.txt");
  //Cuts::sliceCuts[29] = Cuts::mirror( Cuts::sliceCuts[41], 29 );
  //Cuts::WriteCuts("cutsNew.txt");


  tracker.mUseMC = 0; //SG!!
  tracker.wrongHitMix = 1.;

  if( tracker.mRecoPassForLearning >= 0 ){ // learn some parameters

    Learning ML;
    ML.init();

    int nThreadsTotal = std::thread::hardware_concurrency();
    if( nThreadsTotal < 1 ) nThreadsTotal = 1;

    //nThreadsTotal = 10;

    int nTrackers = nThreadsTotal;
    if( nTrackers  > nEvents ) nTrackers = nEvents;

    // set up trackers

    Tracker *trackers = new Tracker[nTrackers];
    std::thread *trackerThreads[nTrackers];
    {
      int nThreadsLeft = nThreadsTotal - nTrackers;
      for( int i=0; i<nTrackers; i++ ){
	trackers[i].mRecoPassForLearning = tracker.mRecoPassForLearning;
	trackers[i].mNThreads = 1;
	trackers[i].doPrint = 0;
	trackers[i].mUseMC = 0; //SG!!
 	trackers[i].wrongHitMix = 0; //SG!!

	trackerThreads[i] = NULL;
      }
      while( nThreadsLeft > 0 ){
	for( int i=0; i<nTrackers && nThreadsLeft>0; i++ ){
	  trackers[i].mNThreads++;
	  nThreadsLeft--;
	}
      }
    }


    for( int startEvent = 0; startEvent<nEvents;  startEvent+=nTrackers ){

      int err[nTrackers];

      for( int itr=0; itr<nTrackers; itr++ ){
	int &terr = err[itr];
	terr = -1;
	int event = startEvent + itr;
	if( event >= nEvents ) break;
	Tracker &tr = trackers[itr];

	trackerThreads[itr] = new std::thread
	  ( [ &tr, &terr, event, firstEvent, dir ]
	    {
	      cout<<"read event N "<<event<<" id "<<firstEvent+event<<" for learning.. "<<endl;
	      EventReader reader( &tr );
	      terr = reader.readEvent( dir,  firstEvent+event, true );
	      if( terr!=0 ) return;
	      tr.reconstruct( 1 ); // stop at the learning point
	    }
	    );
      }
      for( int itr=0; itr<nTrackers; itr++ ){
	if( !trackerThreads[itr] ) break;
	trackerThreads[itr]->join();
	if( err[itr]==0 ) ML.StoreEvent( trackers[itr] );
	delete trackerThreads[itr];
	trackerThreads[itr] = NULL;
      }
    }



    cout<<" Start learning.. "<<endl<<endl;;

    tracker.doPrint = 0;

    RecoPassParameters &cuts = Cuts::sliceCuts[tracker.mRecoPassForLearning];

    /*
    //ML.extendToMin( cuts );
    ML.resetPhi( cuts, 1.5 );
    ML.resetT  ( cuts, 2.5 );
    ML.resetUV( cuts, 2.5 );
    //ML.extendToMin( cuts );
    */
    ML.selectAll(0);
    //ML.selectDead( cuts, 1., 1., 1. );

    ML.SetCalibrationObject( cuts, 1);

    ML.setDefaultCutsClosed();
    //ML.setDefaultCutsOpen( cuts, 1.1, 1.1, 1.1 );

    //ML.setDefaultCuts( cuts, 1. );
    //ML.writeDefaultCuts("cutsDefault.txt");
    //ML.readDefaultCuts("cutsDefault.txt");
    //ML.setDefaultCutsClosed();


    //ML.selectParameters(0);
    ML.startSample=3;


    //ML.selectParameter(0,1); // t1
    ML.selectParameter(1,1); // T1
    ML.selectParameter(2,0); // phi min
    ML.selectParameter(3,0); // phi max
    ML.selectParameter(4,1); // T2
    //ML.selectParameter(5,1); // Pt2
    //ML.selectParameter(6,1); // Pt3


    /*
    ML.selectPhiT(1);
    ML.selectUV(1);
    ML.selectUVms(1);
    */

    ML.selectSearchLayers(1);
    ML.selectPickUpLayers(1);


    //ML.selectSearchLayer(28,1);
    //ML.selectPickUpLayer(11,1);


    /*
    ML.selectVolume(3,1);
    ML.selectVolume(5,1);
    */
    //ML.selectVolume(6,1);
    //ML.selectVolume(8,1);

   //ML.selectPhiT(0);


    //ML.selectParameter(19,0);
    //ML.selectParameter(20,0);
    //ML.selectParameter(21,0);
    //ML.selectParameter(22,0);
   //ML.selectParameter(27,0);
    //ML.selectParameter(28,0);

    //ML.selectParameter(135,0);
    //ML.selectParameter(136,0);


    ML.StartLearning();

    while( ML.startNewTry() ){

      double qa[4] = {0,0,0,0};

      for( uint startEvent = 0; startEvent<ML.mEvents.size();  startEvent+=nTrackers ){
	int err[nTrackers];
	for( int itr=0; itr<nTrackers; itr++ ){
	  int &terr = err[itr];
	  terr = -1;
	  uint event = startEvent + itr;
	  if( event >= ML.mEvents.size() ) break;
	  Tracker &tr = trackers[itr];

	  trackerThreads[itr] = new std::thread
	    ( [ &tr, &terr, event, &ML ]
	      {
		ML.RestoreEvent( event, tr );
		tr.reconstruct( 2 ); // continue from the learning point
		terr=0;
	      }
	      );
	}
	for( int itr=0; itr<nTrackers; itr++ ){
	  if( !trackerThreads[itr] ) break;
	  trackerThreads[itr]->join();
	  if( err[itr]==0 ){
	    Tracker &tr = trackers[itr];
	    qa[0]+=tr.mQA.mTotalEfficiency;//mTotalRemoved - tr.mQA.mTotalEfficiency;
	    qa[1]+=tr.mQA.mEfficiency;
	    qa[2]+=tr.mQA.mPurity;
	    //qa[1]+=100.-tr.mQA.mHitsMissed;
	    qa[3]+=100.-tr.mQA.mHitsFake;
	  }
	  delete trackerThreads[itr];
	  trackerThreads[itr] = NULL;
	}
      }

      for( int i=0; i<4; i++ ) qa[i] /= ML.mEvents.size();
      double eff = (qa[1] + 0.1*qa[2] + 0.*qa[3] );
      ML.PrintStatus( eff );
      ML.setTryResult( eff, qa, 4);
      Cuts::WriteCuts("cutsNewTmp1.txt");
      Cuts::WriteCuts("cutsNewTmp2.txt");
      //break;
    }
    Cuts::WriteCuts("cutsNew.txt");

    delete[] trackers;

    // return 0;

  } // learning


  Cuts::WriteCuts("cutsNew.txt");
  // real reconstruction

  ofstream out("mysubmission.csv");
  if( !out.is_open() ){
    cout<<"Can not open output file"<<endl;
    exit(0);
  }

  out<<"event_id,hit_id,track_id"<<endl;


  long int currentID=1;

  EventReader reader( &tracker );

  for( int event = firstEvent; event<firstEvent+nEvents; event++){
    cout<<"read event N "<<event-firstEvent<<" id "<<event<<endl;
    int err = reader.readEvent( dir,  event, analyseTruth );
    if( err!=0 ) break;

    tracker.reconstruct( 0 ); // don't stop at the learning point

    cout<<"Event "<<event<<": reconstructed "<<tracker.mTracks.size()<<" tracks"<<endl;


    for( uint ih=0; ih<tracker.mHits.size(); ih++ ){
      Hit &h = tracker.mHits[ih];
      h.trackID=0;
    }
    currentID=1;

    for( uint itr=0; itr<tracker.mTracks.size(); itr++ ){
      Track &t = tracker.mTracks[itr];
      if( t.nHits <= 0 ) continue;
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
      out<<event<<","<<ih+1<<","<<h.trackID<<endl;
    }
  } // events

  //tracker.analyzeGeometry(1);

  out.close();
/*
      for (auto it : tracks) {
          if (it.first<=0) continue; // Holds unassigned points
          auto track = it.second;
          if (track.size() == 0) continue;
          if (it.first<MAXTRACK || it.first>nt-MAXTRACK) {
              cout << "Track " << it.first << ": ";
              //for (auto it : track) cout << it << "/" << tracker->shortid(it) << "(" << tracker->truthPart(it)%nhits << ") ";
              tracker->printshort(track);
          }
          if (it.first == MAXTRACK) cout << endl << "..." << endl;
      }
*/
  // ACTS_INFO("Loaded input points for event " << ctx.eventNumber);
  // for (auto& volumeData : *spacePoints) {
  //   for (auto& layerData : volumeData.second) {
  //     for (auto& moduleData : layerData.second) {
  //       ACTS_INFO("volume " << volumeData.first << " layer " << layerData.first
  //                           << " module " << moduleData.first << " "
  //                           << moduleData.second.size() << " space points");
  //       for (auto& cluster : moduleData.second) {
  //         // TODO do something w/ the space points and write out track
  //         // candidates NOTE space points are grouped by the geometry id
  //       }
  //     }
  //   }
  // }

  return ProcessCode::SUCCESS;
}
