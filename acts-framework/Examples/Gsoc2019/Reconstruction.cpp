// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>
#include <memory>
#include <vector>

#include <boost/program_options.hpp>

#include <Acts/Detector/TrackingGeometry.hpp>
#include <Acts/Extrapolator/Navigator.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/Plugins/Digitization/PlanarModuleStepper.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/StraightLineStepper.hpp>
#include <Acts/Utilities/Units.hpp>
#include <Fatras/Kernel/Interactor.hpp>
#include <Fatras/Kernel/SelectorList.hpp>
#include <Fatras/Kernel/Simulator.hpp>
#include <Fatras/Selectors/ChargeSelectors.hpp>
#include <Fatras/Selectors/KinematicCasts.hpp>
#include <Fatras/Selectors/SelectorHelpers.hpp>

#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/Common/CommonOptions.hpp"
#include "ACTFW/Digitization/DigitizationAlgorithm.hpp"
#include "ACTFW/Digitization/DigitizationOptions.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Fatras/FatrasAlgorithm.hpp"
#include "ACTFW/Fatras/FatrasOptions.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/GenericDetector/BuildGenericDetector.hpp"
#include "ACTFW/Options/ParticleGunOptions.hpp"
#include "ACTFW/Random/RandomNumbersOptions.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "../../Algorithms/Reconstruction/ACTFW/DAGbasedNNTracker/main.hpp"
#include "../../Algorithms/Reconstruction/ACTFW/MikadoTracker/main.hpp"
#include "../../Algorithms/Reconstruction/ACTFW/Top-quarks/main.hpp"
#include "../../Algorithms/Reconstruction/ACTFW/pytorch/main.hpp"


using namespace boost::program_options;
using namespace FW;

// data type definitions
using TruthParticle   = Data::SimParticle;
using TruthHit        = Data::SimHit<Data::SimParticle>;
using TruthHitCreator = Data::SimHitCreator;
using TruthEvent      = std::vector<Data::SimVertex<TruthParticle>>;

// type definitions for the propagation and simulation
// for charge particles
using MagneticField     = Acts::ConstantBField;
using ChargedStepper    = Acts::EigenStepper<MagneticField>;
using ChargedPropagator = Acts::Propagator<ChargedStepper, Acts::Navigator>;
using ChargedInteractor = Fatras::
    Interactor<RandomEngine, TruthParticle, TruthHit, TruthHitCreator>;
// for neutral particles
using NeutralStepper    = Acts::StraightLineStepper;
using NeutralPropagator = Acts::Propagator<NeutralStepper, Acts::Navigator>;
using NeutralInteractor = Fatras::
    Interactor<RandomEngine, TruthParticle, TruthHit, TruthHitCreator>;
// the simulator kernel
using Simulator = Fatras::Simulator<ChargedPropagator,
                                    Fatras::SelectorListOR<>,
                                    ChargedInteractor,
                                    NeutralPropagator,
                                    Fatras::SelectorListOR<>,
                                    NeutralInteractor>;
// the combined simulation algorithm
using SimulationAlgorithm = FatrasAlgorithm<Simulator, TruthEvent, TruthHit>;

SimulationAlgorithm::Config
buildSimulationConfig(std::shared_ptr<const Acts::TrackingGeometry>& geo,
                      Acts::ConstantBField&& magneticField)
{
  using namespace Acts;

  ChargedStepper    chargedStepper(std::move(magneticField));
  ChargedPropagator chargedProp(std::move(chargedStepper), Navigator(geo));
  NeutralPropagator neutralProp(NeutralStepper{}, Navigator(geo));
  Simulator         simulator(std::move(chargedProp), std::move(neutralProp));

  return SimulationAlgorithm::Config(std::move(simulator));
}

int
main(int argc, char* argv[])
{
  // Setup command line arguments and options
  options_description desc("Gsoc2019 reconstruction tool");
  Options::addCommonOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addParticleGunOptions(desc);

  // Process command line arguments and options
  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);
  // Print help if reqested
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return EXIT_SUCCESS;
  }

  // Extract common options
  auto nEvents  = Options::readNumberOfEvents(vm);
  auto logLevel = Options::readLogLevel(vm);

  // Create a sequencer w/ the default config
  Sequencer seq({});

  // Setup basic services
  RandomNumbersSvc::Config rndCfg;
  rndCfg.seed  = 123;
  auto rng     = std::make_shared<RandomNumbersSvc>(rndCfg);
  auto barcode = std::make_shared<BarcodeSvc>(
      BarcodeSvc::Config{}, Acts::getDefaultLogger("BarcodeSvc", logLevel));

  // Always use the generic (TrackML) detector
  std::shared_ptr<const Acts::TrackingGeometry> geo(
      Generic::buildGenericDetector(logLevel));
  // Always use a constant magnetic field
  Acts::ConstantBField bfield(0.0, 0.0, 2.0 * Acts::units::_T);

  // Event generation w/ particle gun
  EventGenerator::Config evgen = Options::readParticleGunOptions(vm);
  evgen.randomNumbers          = rng;
  evgen.barcodeSvc             = barcode;
  evgen.output                 = "particles";
  seq.addReaders({std::make_shared<EventGenerator>(evgen, logLevel)});

  // Setup fast simulation
  auto sim                     = buildSimulationConfig(geo, std::move(bfield));
  sim.randomNumberSvc          = rng;
  sim.inputEventCollection     = evgen.output;
  sim.simulatedEventCollection = "simulated_particles";
  sim.simulatedHitCollection   = "simulated_hits";
  seq.appendEventAlgorithms(
      {std::make_shared<SimulationAlgorithm>(sim, logLevel)});

  // Setup digitization to convert truth hits into digitized hits/ space points
  DigitizationAlgorithm::Config digi;
  digi.planarModuleStepper = std::make_shared<Acts::PlanarModuleStepper>(
      Acts::getDefaultLogger("PlanarModuleStepper", logLevel));
  digi.randomNumberSvc        = rng;
  digi.simulatedHitCollection = sim.simulatedHitCollection;
  digi.spacePointCollection   = "spacepoints";
  digi.clusterCollection      = "clusters";
  seq.appendEventAlgorithms(
      {std::make_shared<DigitizationAlgorithm>(digi, logLevel)});

  // Add the empty reconstruction algorithm
  // DAGbasedNNTracker::Config Reco1;
  // Reco1.spacePointCollection = digi.spacePointCollection;
  // seq.appendEventAlgorithms(
  //     {std::make_shared<DAGbasedNNTracker>(Reco1, logLevel)});
  //
  MikadoTracker::Config Reco2;
  Reco2.spacePointCollection = digi.spacePointCollection;
  seq.appendEventAlgorithms(
      {std::make_shared<MikadoTracker>(Reco2, logLevel)});

  // Topquarks::Config Reco3;
  // Reco3.spacePointCollection = digi.spacePointCollection;
  // seq.appendEventAlgorithms(
  //     {std::make_shared<Topquarks>(Reco3, logLevel)});
  // PytorchReconstruction::Config Reco3;
  // Reco3.spacePointCollection = digi.spacePointCollection;
  // seq.appendEventAlgorithms(
  //     {std::make_shared<PytorchReconstruction>(Reco3, logLevel)});

  return (seq.run(nEvents) == ProcessCode::SUCCESS) ? EXIT_SUCCESS
                                                    : EXIT_FAILURE;
}
