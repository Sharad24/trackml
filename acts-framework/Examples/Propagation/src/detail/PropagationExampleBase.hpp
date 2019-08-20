// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/program_options.hpp>
#include <memory>
#include "ACTFW/Common/CommonOptions.hpp"
#include "ACTFW/Common/GeometryOptions.hpp"
#include "ACTFW/Common/OutputOptions.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Plugins/Obj/ObjPropagationStepsWriter.hpp"
#include "ACTFW/Plugins/Root/RootPropagationStepsWriter.hpp"
#include "ACTFW/Propagation/PropagationAlgorithm.hpp"
#include "ACTFW/Propagation/PropagationOptions.hpp"
#include "ACTFW/Random/RandomNumbersOptions.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"

namespace po = boost::program_options;

/// @brief Propgation setup
///
/// @tparam sequencer_t Type of the sequencer of the framework
/// @tparam bfield_t Type of the magnetic field
///
/// @param sequencer The framework sequencer, Propgation algorithm to be added
/// @param bfield The bfield object needed for the Stepper & propagagor
/// @param vm The program options for the log file
/// @param randomNumberSvc The framework random number engine
/// @param tGeometry The TrackingGeometry object
///
/// @return a process code
template <typename sequencer_t, typename bfield_t>
FW::ProcessCode
setupPropgation(sequencer_t&                                  sequencer,
                bfield_t                                      bfield,
                po::variables_map&                            vm,
                std::shared_ptr<FW::RandomNumbersSvc>         randomNumberSvc,
                std::shared_ptr<const Acts::TrackingGeometry> tGeometry)
{
  // Get the log level
  auto logLevel = FW::Options::readLogLevel<po::variables_map>(vm);

  // Get a Navigator
  Acts::Navigator navigator(tGeometry);

  // Resolve the bfield map template and create the propgator
  using Stepper    = Acts::EigenStepper<bfield_t>;
  using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
  Stepper    stepper(std::move(bfield));
  Propagator propagator(std::move(stepper), std::move(navigator));

  // Read the propagation config and create the algorithms
  auto pAlgConfig = FW::Options::readPropagationConfig(vm, propagator);
  pAlgConfig.randomNumberSvc = randomNumberSvc;
  auto propagationAlg = std::make_shared<FW::PropagationAlgorithm<Propagator>>(
      pAlgConfig, logLevel);

  // Add the propagation algorithm
  sequencer.appendEventAlgorithms({propagationAlg});

  return FW::ProcessCode::SUCCESS;
}

/// @brief Straight Line Propgation setup
///
/// @tparam sequencer_t Type of the sequencer of the framework
///
/// @param sequencer The framework sequencer, Propgation algorithm to be added
/// @param vm The program options for the log file
/// @param randomNumberSvc The framework random number engine
/// @param tGeometry The TrackingGeometry object
///
/// @return a process code
template <typename sequencer_t>
FW::ProcessCode
setupStraightLinePropgation(
    sequencer_t&                                  sequencer,
    po::variables_map&                            vm,
    std::shared_ptr<FW::RandomNumbersSvc>         randomNumberSvc,
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry)
{
  // Get the log level
  auto logLevel = FW::Options::readLogLevel<po::variables_map>(vm);

  // Get a Navigator
  Acts::Navigator navigator(tGeometry);

  // Straight line stepper
  using SlStepper  = Acts::StraightLineStepper;
  using Propagator = Acts::Propagator<SlStepper, Acts::Navigator>;
  // Make stepper and propagator
  SlStepper  stepper;
  Propagator propagator(std::move(stepper), std::move(navigator));

  // Read the propagation config and create the algorithms
  auto pAlgConfig = FW::Options::readPropagationConfig(vm, propagator);
  pAlgConfig.randomNumberSvc = randomNumberSvc;
  auto propagationAlg = std::make_shared<FW::PropagationAlgorithm<Propagator>>(
      pAlgConfig, logLevel);

  // Add the propagation algorithm
  sequencer.appendEventAlgorithms({propagationAlg});

  return FW::ProcessCode::SUCCESS;
}

/// The Propagation example
///
/// @tparam geometry_getter_t Type of the geometry getter struct
///
/// @param argc the number of argumetns of the call
/// @param atgv the argument list
/// @param trackingGeometry is the access struct for the trackingGeometry
///
template <typename geometry_options_t, typename geometry_getter_t>
int
propagationExample(int                argc,
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
  // Add the common options
  FW::Options::addCommonOptions<po::options_description>(desc);
  // Add the geometry options
  FW::Options::addGeometryOptions<po::options_description>(desc);
  // Add the bfield options
  FW::Options::addBFieldOptions<po::options_description>(desc);
  // Add the random number options
  FW::Options::addRandomNumbersOptions<po::options_description>(desc);
  // Add the fatras options
  FW::Options::addPropagationOptions<po::options_description>(desc);
  // Add the output options
  FW::Options::addOutputOptions<po::options_description>(desc);

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

  // The Log level
  auto nEvents  = FW::Options::readNumberOfEvents<po::variables_map>(vm);
  auto logLevel = FW::Options::readLogLevel<po::variables_map>(vm);

  // Create the random number engine
  auto randomNumberSvcCfg
      = FW::Options::readRandomNumbersConfig<po::variables_map>(vm);
  auto randomNumberSvc
      = std::make_shared<FW::RandomNumbersSvc>(randomNumberSvcCfg);
  // Add it to the sequencer
  sequencer.addServices({randomNumberSvc});

  // Get the tracking geometry
  auto tGeometry = trackingGeometry(vm);

  // Create BField service
  auto bField = FW::Options::readBField<po::variables_map>(vm);

  if (vm["prop-stepper"].template as<int>() == 0) {
    // Straight line stepper was chosen
    setupStraightLinePropgation(sequencer, vm, randomNumberSvc, tGeometry);
  } else if (bField.first) {
    // Define the interpolated b-field
    using BField = Acts::SharedBField<Acts::InterpolatedBFieldMap>;
    BField fieldMap(bField.first);
    setupPropgation(sequencer, fieldMap, vm, randomNumberSvc, tGeometry);
  } else {
    // Create the constant  field
    using CField = Acts::ConstantBField;
    CField fieldMap(*bField.second);
    setupPropgation(sequencer, fieldMap, vm, randomNumberSvc, tGeometry);
  }

  // ---------------------------------------------------------------------------------
  // Output directory
  std::string outputDir    = vm["output-dir"].template as<std::string>();
  auto        psCollection = vm["prop-step-collection"].as<std::string>();

  if (vm["output-root"].template as<bool>()) {
    // Write the propagation steps as ROOT TTree
    FW::Root::RootPropagationStepsWriter::Config pstepWriterRootConfig;
    pstepWriterRootConfig.collection = psCollection;
    pstepWriterRootConfig.filePath
        = FW::joinPaths(outputDir, psCollection + ".root");
    auto pstepWriterRoot
        = std::make_shared<FW::Root::RootPropagationStepsWriter>(
            pstepWriterRootConfig);
    if (sequencer.addWriters({pstepWriterRoot}) != FW::ProcessCode::SUCCESS)
      return -1;
  }

  if (vm["output-obj"].template as<bool>()) {

    using PropagationSteps = Acts::detail::Step;
    using ObjPropagationStepsWriter
        = FW::Obj::ObjPropagationStepsWriter<PropagationSteps>;

    // Write the propagation steps as Obj TTree
    ObjPropagationStepsWriter::Config pstepWriterObjConfig;
    pstepWriterObjConfig.collection = psCollection;
    pstepWriterObjConfig.outputDir  = outputDir;
    auto pstepWriteObj
        = std::make_shared<ObjPropagationStepsWriter>(pstepWriterObjConfig);
    if (sequencer.addWriters({pstepWriteObj}) != FW::ProcessCode::SUCCESS)
      return -1;
  }

  // Initiate the run
  sequencer.run(nEvents);
  // Return success code
  return 0;
}
