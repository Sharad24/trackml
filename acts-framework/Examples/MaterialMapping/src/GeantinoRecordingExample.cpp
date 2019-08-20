// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>
#include "ACTFW/Common/CommonOptions.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/MaterialMapping/GeantinoRecording.hpp"
#include "ACTFW/Plugins/DD4hep/DD4hepDetectorOptions.hpp"
#include "ACTFW/Plugins/DD4hep/GeometryService.hpp"
#include "ACTFW/Plugins/DD4hepG4/DD4hepToG4Svc.hpp"
#include "ACTFW/Plugins/Root/RootMaterialTrackWriter.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "ACTFW/Writers/IWriterT.hpp"
#include "Acts/Detector/TrackingGeometry.hpp"

namespace po = boost::program_options;

int
main(int argc, char* argv[])
{
  // Declare the supported program options.
  po::options_description desc("Allowed options");
  // add the standard options
  FW::Options::addStandardOptions<po::options_description>(desc, 100, 2);
  // add the detector options
  FW::Options::addDD4hepOptions<po::options_description>(desc);
  po::variables_map vm;
  // Get all options from contain line and store it into the map
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  // print help if requested
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }
  // now read the standard options
  auto standardOptions
      = FW::Options::readStandardOptions<po::variables_map>(vm);
  // @todo update - make program options in separate MR
  auto   nEvents     = standardOptions.first;
  size_t nTracks     = 100;
  int    randomSeed1 = 536235167;
  int    randomSeed2 = 729237523;

  // DETECTOR:
  // --------------------------------------------------------------------------------
  // DD4Hep detector definition
  // read the detector config & dd4hep detector
  auto dd4HepDetectorConfig
      = FW::Options::readDD4hepConfig<po::variables_map>(vm);
  auto geometrySvc
      = std::make_shared<FW::DD4hep::GeometryService>(dd4HepDetectorConfig);
  std::shared_ptr<const Acts::TrackingGeometry> tGeometry
      = geometrySvc->trackingGeometry();

  // DD4Hep to Geant4 conversion
  //
  FW::DD4hepG4::DD4hepToG4Svc::Config dgConfig("DD4hepToG4",
                                               Acts::Logging::INFO);
  dgConfig.dd4hepService = geometrySvc;
  auto dd4hepToG4Svc = std::make_shared<FW::DD4hepG4::DD4hepToG4Svc>(dgConfig);

  // --------------------------------------------------------------------------------
  // Geant4 JOB:
  // --------------------------------------------------------------------------------
  // set up the writer for
  FW::Root::RootMaterialTrackWriter::Config g4WriterConfig(
      "MaterialTrackWriter", Acts::Logging::INFO);
  g4WriterConfig.fileName = "GeantMaterialTracks.root";
  g4WriterConfig.treeName = "GeantMaterialTracks";
  auto g4TrackRecWriter
      = std::make_shared<FW::Root::RootMaterialTrackWriter>(g4WriterConfig);

  // set up the algorithm writing out the material map
  FW::GeantinoRecording::Config g4rConfig;
  g4rConfig.materialTrackWriter = g4TrackRecWriter;
  g4rConfig.geant4Service       = dd4hepToG4Svc;
  g4rConfig.tracksPerEvent      = nTracks;
  g4rConfig.seed1               = randomSeed1;
  g4rConfig.seed2               = randomSeed2;
  // create the geant4 algorithm
  auto g4rAlgorithm
      = std::make_shared<FW::GeantinoRecording>(g4rConfig, Acts::Logging::INFO);

  // Geant4 job - these can be many Geant4 jobs, indeed
  //
  // create the config object for the sequencer
  FW::Sequencer::Config g4SeqConfig;
  // now create the sequencer
  FW::Sequencer g4Sequencer(g4SeqConfig);
  // the writer is a service as it needs initialize, finalize
  g4Sequencer.addServices({g4TrackRecWriter});
  g4Sequencer.appendEventAlgorithms({g4rAlgorithm});
  g4Sequencer.run(nEvents);
}
