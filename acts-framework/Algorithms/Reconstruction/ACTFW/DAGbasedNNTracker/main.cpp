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

#define VERBOSE false
#define MAXTRACK 10
#define FILENUM 21101

FW::DAGbasedNNTracker::DAGbasedNNTracker(
    const FW::DAGbasedNNTracker::Config& cfg,
    Acts::Logging::Level                            logLevel)
  : BareAlgorithm("DAGbasedNNTracker", logLevel), m_cfg(cfg)
{
  if (m_cfg.spacePointCollection.empty()) {
    throw std::invalid_argument("Missing input space points collection");
  }
}

FW::ProcessCode
FW::DAGbasedNNTracker::execute(FW::AlgorithmContext ctx) const
{
  ACTS_INFO("empty reconstruction on event " << ctx.eventNumber);

  // initialize point struct
  int filenum = FILENUM;
  int number = 1;
  int step = 1;
  for (int n=0;n<number;n++,filenum++) {

      ACTS_INFO("Running on event #" << filenum << "\n");
      ACTS_INFO("Step " << step << "\n");

      Tracker *tracker = new Tracker(filenum,DATAPATH,WORKPATH);

      if (!tracker->evaluation()) {
          //tracker->readBlacklist();
          //tracker->readWhitelist();
          tracker->readStarts();
          tracker->readTruth();
          tracker->sortTracks();
      }

      tracker->readHits();
      tracker->readCells();

      tracker->setDebug(VERBOSE);
      int *assignment = tracker->findTracks(step,WORKPATH);

      // Assemble tracks
      ACTS_INFO("Assembling tracks..." << "\n");
      std::map<int,std::vector<int> > tracks;

      long nhits = tracker->numberHits();
      int is = VERBOSE ? 0 : 1; // track 0 holds the unassigned points
      for(int i=is;i<nhits;i++) {
          int track = assignment[i];
          tracks[track].push_back(i);
      }

      long nt = tracks.size();
      ACTS_INFO("\n" << "Number of tracks: " << nt << "\n");
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
      delete tracker;
  }
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
