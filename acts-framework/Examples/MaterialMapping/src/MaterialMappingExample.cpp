// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>
#include "ACTFW/Common/CommonOptions.hpp"
#include "ACTFW/Extrapolation/ExtrapolationUtils.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/MaterialMapping/MaterialMapping.hpp"
#include "ACTFW/Plugins/DD4hep/DD4hepDetectorOptions.hpp"
#include "ACTFW/Plugins/DD4hep/GeometryService.hpp"
#include "ACTFW/Plugins/DD4hepG4/DD4hepToG4Svc.hpp"
#include "ACTFW/Plugins/Root/RootIndexedMaterialWriter.hpp"
#include "ACTFW/Plugins/Root/RootMaterialTrackReader.hpp"
#include "ACTFW/Plugins/Root/RootMaterialTrackWriter.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Plugins/MaterialMapping/MaterialMapper.hpp"

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
  // read the detector config & dd4hep detector
  auto dd4HepDetectorConfig
      = FW::Options::readDD4hepConfig<po::variables_map>(vm);
  auto geometrySvc
      = std::make_shared<FW::DD4hep::GeometryService>(dd4HepDetectorConfig);

  std::shared_ptr<const Acts::TrackingGeometry> tGeometry
      = geometrySvc->trackingGeometry();

  auto nEvents = standardOptions.first;

  // DD4Hep to Geant4 conversion
  //
  FW::DD4hepG4::DD4hepToG4Svc::Config dgConfig("DD4hepToG4",
                                               Acts::Logging::INFO);
  dgConfig.dd4hepService = geometrySvc;
  auto dd4hepToG4Svc = std::make_shared<FW::DD4hepG4::DD4hepToG4Svc>(dgConfig);

  // @todo update - make program options in separate MR
  // --------------------------------------------------------------------------------
  // MaterialMapping Algorithm configruation:
  //
  // set up the writer for the surface material maps
  FW::Root::RootMaterialTrackReader::Config mtrReaderConfig(
      "MaterialTrackReader", Acts::Logging::DEBUG);
  mtrReaderConfig.fileList
      = {"GeantMaterialTracks0.root", "GeantMaterialTracks1.root"};
  //                               "GeantMaterialTracks2.root",
  //                               "GeantMaterialTracks3.root",
  //                               "GeantMaterialTracks4.root"};
  mtrReaderConfig.treeName = "GeantMaterialTracks";
  auto mtrReader
      = std::make_shared<FW::Root::RootMaterialTrackReader>(mtrReaderConfig);

  // EXTRAPOLATOR - set up the extrapolator
  // set up the magnetic field
  std::shared_ptr<Acts::ConstantBField> magFieldSvc(
      new Acts::ConstantBField{{0., 0., 0.002}});  // field is given in kT
  // EXTRAPOLATOR - set up the extrapolator
  std::shared_ptr<Acts::IExtrapolationEngine> extrapolationEngine
      = FW::initExtrapolator(tGeometry, magFieldSvc, Acts::Logging::INFO);

  // create material mapping
  Acts::MaterialMapper::Config mapperConf;
  mapperConf.extrapolationEngine = extrapolationEngine;
  auto mtrMapper                 = std::make_shared<Acts::MaterialMapper>(
      mapperConf,
      Acts::getDefaultLogger("MaterialMapper", Acts::Logging::DEBUG));

  // create the mapped material writer
  // set up the writer for
  FW::Root::RootMaterialTrackWriter::Config mtrWriterConfig(
      "MappedMaterialTrackWriter", Acts::Logging::INFO);
  mtrWriterConfig.fileName = "MappedMaterialTracks.root";
  mtrWriterConfig.treeName = "MappedMaterialTracks";
  auto mtrWriter
      = std::make_shared<FW::Root::RootMaterialTrackWriter>(mtrWriterConfig);

  // create the material writer
  FW::Root::RootIndexedMaterialWriter::Config imatWriterConfig(
      "MaterialWriter", Acts::Logging::INFO);

  imatWriterConfig.fileName = "$PWD/LayerMaterialMaps.root";
  auto imaterialWriter
      = std::make_shared<FW::Root::RootIndexedMaterialWriter>(imatWriterConfig);

  // set up the algorithm reading in the material map and mapping the material
  // onto the tracking geometry
  FW::MaterialMapping::Config mmConfig;
  mmConfig.materialTrackReader   = mtrReader;
  mmConfig.materialTrackWriter   = mtrWriter;
  mmConfig.materialMapper        = mtrMapper;
  mmConfig.indexedMaterialWriter = imaterialWriter;
  mmConfig.trackingGeometry      = tGeometry;
  mmConfig.maximumTrackRecords   = 10000;
  auto materialMappingAlg
      = std::make_shared<FW::MaterialMapping>(mmConfig, Acts::Logging::INFO);

  // --------------------------------------------------------------------------------
  // Mapping job configruation
  //
  // create the config object for the sequencer
  FW::Sequencer::Config mapSeqConfig;
  // now create the sequencer
  FW::Sequencer mappingSequencer(mapSeqConfig);
  mappingSequencer.addServices({mtrReader, mtrWriter, imaterialWriter});
  mappingSequencer.appendEventAlgorithms({materialMappingAlg});
  mappingSequencer.run(nEvents);
}
