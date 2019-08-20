// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include <boost/program_options.hpp>

#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/Common/CommonOptions.hpp"
#include "ACTFW/Common/GeometryOptions.hpp"
#include "ACTFW/Common/OutputOptions.hpp"
#include "ACTFW/Digitization/DigitizationOptions.hpp"
#include "ACTFW/Fatras/FatrasOptions.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Options/ParticleGunOptions.hpp"
#include "ACTFW/Options/Pythia8Options.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Plugins/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Random/RandomNumbersOptions.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "FatrasDigitizationBase.hpp"
#include "FatrasEvgenBase.hpp"
#include "FatrasSimulationBase.hpp"

namespace po = boost::program_options;

/// @brief The Fatras example
///
/// @tparam geometry_getter_t Type of the geometry getter struct
///
/// @param argc the number of argumetns of the call
/// @param aegv the argument list
template <typename geometry_options_t, typename geometry_getter_t>
int
fatrasExample(int                argc,
              char*              argv[],
              geometry_options_t geometryOptions,
              geometry_getter_t  trackingGeometry)
{
  // Create the config object for the sequencer
  FW::Sequencer::Config seqConfig;
  // Now create the sequencer
  FW::Sequencer sequencer(seqConfig);
  // Declare the supported program options.
  po::options_description desc("Allowed options");
  // Add the Common options
  FW::Options::addCommonOptions<po::options_description>(desc);
  // Add the geometry options
  FW::Options::addGeometryOptions<po::options_description>(desc);
  // Add the particle gun options
  FW::Options::addParticleGunOptions(desc);
  // Add the Pythia 8 options
  FW::Options::addPythia8Options(desc);
  // Add the random number options
  FW::Options::addRandomNumbersOptions<po::options_description>(desc);
  // Add the bfield options
  FW::Options::addBFieldOptions<po::options_description>(desc);
  // Add the fatras options
  FW::Options::addFatrasOptions<po::options_description>(desc);
  // Add the digization options
  FW::Options::addDigitizationOptions<po::options_description>(desc);
  // Add the output options
  FW::Options::addOutputOptions<po::options_description>(desc);
  // Add program specific options: input / output
  desc.add_options()("evg-input-type",
                     po::value<std::string>()->default_value("pythia8"),
                     "Type of evgen input 'gun', 'pythia8'");

  // Add specific options for this geometry
  geometryOptions(desc);

  // Map to store the given program options
  po::variables_map vm;
  // Get all options from contain line and store it into the map
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  // Print help if requested
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }
  // Read the common options : number of events and log level
  auto nEvents  = FW::Options::readNumberOfEvents<po::variables_map>(vm);
  auto logLevel = FW::Options::readLogLevel<po::variables_map>(vm);

  // Create the random number engine
  auto randomNumberSvcCfg
      = FW::Options::readRandomNumbersConfig<po::variables_map>(vm);
  auto randomNumberSvc
      = std::make_shared<FW::RandomNumbersSvc>(randomNumberSvcCfg);

  // Add it to the sequencer
  sequencer.addServices({randomNumberSvc});
  // Create the barcode service
  FW::BarcodeSvc::Config barcodeSvcCfg;
  auto                   barcodeSvc = std::make_shared<FW::BarcodeSvc>(
      barcodeSvcCfg, Acts::getDefaultLogger("BarcodeSvc", logLevel));
  // Add it to the sequencer
  sequencer.addServices({barcodeSvc});

  // Get the tracking geometry
  auto tGeometry = trackingGeometry(vm);

  // (A) EVGEN
  // Setup the evgen input to the simulation
  setupEvgenInput<po::variables_map>(
      vm, sequencer, barcodeSvc, randomNumberSvc);

  // (B) SIMULATION
  // Setup the simulation
  setupSimulation<po::variables_map>(
      vm, sequencer, tGeometry, barcodeSvc, randomNumberSvc);

  // (C) DIGITIZATION
  // Setup the digitization
  setupDigitization<po::variables_map>(
      vm, sequencer, barcodeSvc, randomNumberSvc);

  // (D) TRUTH TRACKING

  // (E) PATTERN RECOGNITION

  // Initiate the run
  sequencer.run(nEvents);
  // Return 0 for success
  return 0;
}
