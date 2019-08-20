// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Fatras/FatrasAlgorithm.hpp"
#include "ACTFW/Fatras/FatrasOptions.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Plugins/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Plugins/Root/RootParticleWriter.hpp"
#include "ACTFW/Plugins/Root/RootSimHitWriter.hpp"
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
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/GeometryID.hpp"
#include "Fatras/Kernel/Interactor.hpp"
#include "Fatras/Kernel/SelectorList.hpp"
#include "Fatras/Kernel/Simulator.hpp"
#include "Fatras/Selectors/ChargeSelectors.hpp"
#include "Fatras/Selectors/KinematicCasts.hpp"
#include "Fatras/Selectors/SelectorHelpers.hpp"

typedef FW::Data::SimHit<FW::Data::SimParticle> FatrasHit;
typedef std::vector<FW::Data::SimVertex<>>      FatrasEvent;

/// Simple struct to select sensitive surfaces
struct SurfaceSelector
{
  bool selectSensitive = true;
  bool selectMaterial  = false;
  bool selectPassive   = false;

  /// SurfaceSelector with options
  ///
  /// @param sSensitive is the directive to select sensitive surfaces
  /// @param sMaterial is the directive to select material surfaces
  /// @param sPassive is the directive to select passive surfaces
  SurfaceSelector(bool sSensitive = true,
                  bool sMaterial  = false,
                  bool sPassive   = false)
    : selectSensitive(sSensitive)
    , selectMaterial(sMaterial)
    , selectPassive(sPassive)
  {
  }

  /// Call operator to check if a surface should be selected
  ///
  /// @param surface is the test surface
  bool
  operator()(const Acts::Surface& surface) const
  {
    if (selectSensitive && surface.associatedDetectorElement()) {
      return true;
    }
    if (selectMaterial && surface.associatedMaterial()) {
      return true;
    }
    if (selectPassive) {
      return true;
    }
    return false;
  }
};

/// @brief Simulation setup for the FatrasAlgorithm
///
/// @tparam bfield_t Type of the bfield for the simulation to be set up
///
/// @param fieldMap The field map for the simulation setup
/// @param sequencer The framework sequencer
/// @param vm The boost variable map to resolve
/// @param tGeometry The TrackingGeometry for the tracking setup
/// @param barcodesSvc The barcode service to be used for the simulation
/// @param randomNumberSvc The random number service to be used for the
/// simulation
template <typename bfield_t>
void
setupSimulationAlgorithm(
    bfield_t                                      fieldMap,
    FW::Sequencer&                                sequencer,
    po::variables_map&                            vm,
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
    std::shared_ptr<FW::BarcodeSvc>               barcodeSvc,
    std::shared_ptr<FW::RandomNumbersSvc>         randomNumberSvc)
{
  // Read the log level
  Acts::Logging::Level logLevel = FW::Options::readLogLevel(vm);

  /// Read the evgen particle collection
  std::string evgenCollection = "particles";

  // Create a navigator for this tracking geometry
  Acts::Navigator cNavigator(tGeometry);
  Acts::Navigator nNavigator(tGeometry);

  // using ChargedStepper     = Acts::AtlasStepper<bfield_t>;
  using ChargedStepper    = Acts::EigenStepper<bfield_t>;
  using ChargedPropagator = Acts::Propagator<ChargedStepper, Acts::Navigator>;
  using NeutralStepper    = Acts::StraightLineStepper;
  using NeutralPropagator = Acts::Propagator<NeutralStepper, Acts::Navigator>;

  ChargedStepper    cStepper(std::move(fieldMap));
  ChargedPropagator cPropagator(std::move(cStepper), std::move(cNavigator));
  NeutralStepper    nStepper;
  NeutralPropagator nPropagator(std::move(nStepper), std::move(nNavigator));

  // The Selector for charged particles, including kinematic cuts
  typedef Fatras::ChargedSelector            CSelector;
  typedef Fatras::Max<Fatras::casts::absEta> CMaxEtaAbs;
  typedef Fatras::Min<Fatras::casts::pT>     CMinPt;
  typedef Fatras::SelectorListAND<CSelector, CMinPt, CMaxEtaAbs>
      ChargedSelector;

  typedef Fatras::NeutralSelector            NSelector;
  typedef Fatras::Max<Fatras::casts::absEta> NMaxEtaAbs;
  typedef Fatras::Min<Fatras::casts::E>      NMinE;
  typedef Fatras::SelectorListAND<NSelector, NMinE, NMaxEtaAbs> NeutralSelector;

  typedef Fatras::PhysicsList<> PhysicsList;

  typedef Fatras::Interactor<FW::RandomEngine,
                             FW::Data::SimParticle,
                             FW::Data::SimHit<FW::Data::SimParticle>,
                             FW::Data::SimHitCreator,
                             SurfaceSelector,
                             PhysicsList>
      ChargedInteractor;

  typedef Fatras::Interactor<FW::RandomEngine,
                             FW::Data::SimParticle,
                             FW::Data::SimHit<FW::Data::SimParticle>,
                             FW::Data::SimHitCreator>
      NeutralInteractor;

  typedef Fatras::Simulator<ChargedPropagator,
                            ChargedSelector,
                            ChargedInteractor,
                            NeutralPropagator,
                            NeutralSelector,
                            NeutralInteractor>
      FatrasSimulator;

  FatrasSimulator fatrasSimulator(cPropagator, nPropagator);
  fatrasSimulator.debug = vm["fatras-debug-output"].template as<bool>();

  using FatrasAlgorithm
      = FW::FatrasAlgorithm<FatrasSimulator, FatrasEvent, FatrasHit>;

  typename FatrasAlgorithm::Config fatrasConfig
      = FW::Options::readFatrasConfig<po::variables_map,
                                      FatrasSimulator,
                                      FatrasEvent,
                                      FatrasHit>(vm, fatrasSimulator);
  fatrasConfig.randomNumberSvc      = randomNumberSvc;
  fatrasConfig.inputEventCollection = evgenCollection;

  // Finally the fatras algorithm
  auto fatrasAlgorithm
      = std::make_shared<FatrasAlgorithm>(fatrasConfig, logLevel);

  // Finalize the squencer setting and run
  sequencer.appendEventAlgorithms({fatrasAlgorithm});

  // Output directory
  std::string outputDir = vm["output-dir"].template as<std::string>();

  // Write simulation information as CSV files
  std::shared_ptr<FW::Csv::CsvParticleWriter> pWriterCsv = nullptr;
  if (vm["output-csv"].template as<bool>()) {
    FW::Csv::CsvParticleWriter::Config pWriterCsvConfig;
    pWriterCsvConfig.collection = fatrasConfig.simulatedEventCollection;
    pWriterCsvConfig.outputDir  = outputDir;
    pWriterCsvConfig.outputFileName
        = fatrasConfig.simulatedEventCollection + ".csv";
    auto pWriterCsv
        = std::make_shared<FW::Csv::CsvParticleWriter>(pWriterCsvConfig);

    sequencer.addWriters({pWriterCsv});
  }

  // Write simulation information as ROOT files
  std::shared_ptr<FW::Root::RootParticleWriter> pWriterRoot = nullptr;
  if (vm["output-root"].template as<bool>()) {
    // Write particles as ROOT TTree
    FW::Root::RootParticleWriter::Config pWriterRootConfig;
    pWriterRootConfig.collection = fatrasConfig.simulatedEventCollection;
    pWriterRootConfig.filePath   = FW::joinPaths(
        outputDir, fatrasConfig.simulatedEventCollection + ".root");
    pWriterRootConfig.treeName   = fatrasConfig.simulatedEventCollection;
    pWriterRootConfig.barcodeSvc = barcodeSvc;
    auto pWriterRoot
        = std::make_shared<FW::Root::RootParticleWriter>(pWriterRootConfig);

    // Write simulated hits as ROOT TTree
    FW::Root::RootSimHitWriter::Config fhitWriterRootConfig;
    fhitWriterRootConfig.collection = fatrasConfig.simulatedHitCollection;
    fhitWriterRootConfig.filePath   = FW::joinPaths(
        outputDir, fatrasConfig.simulatedHitCollection + ".root");
    fhitWriterRootConfig.treeName = fatrasConfig.simulatedHitCollection;
    auto fhitWriterRoot
        = std::make_shared<FW::Root::RootSimHitWriter>(fhitWriterRootConfig);

    sequencer.addWriters({pWriterRoot, fhitWriterRoot});
  }
}

/// @brief Simulation setup
///
/// @tparam bfield_t Type of the bfield for the simulation to be set up
///
/// @param fieldMap The field map for the simulation setup
/// @param sequencer The framework sequencer
/// @param vm The boost variable map to resolve
/// @param tGeometry The TrackingGeometry for the tracking setup
/// @param barcodesSvc The barcode service to be used for the simulation
/// @param randomNumberSvc The random number service to be used for the
/// simulation
template <typename vmap_t>
void
setupSimulation(vmap_t&                                       vm,
                FW::Sequencer&                                sequencer,
                std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
                std::shared_ptr<FW::BarcodeSvc>               barcodeSvc,
                std::shared_ptr<FW::RandomNumbersSvc>         randomNumberSvc)
{
  // create BField service
  auto bField = FW::Options::readBField<vmap_t>(vm);
  // a charged propagator
  if (bField.first) {
    // create the shared field
    using BField = Acts::SharedBField<Acts::InterpolatedBFieldMap>;
    BField fieldMap(bField.first);
    // now setup of the simulation and append it to the sequencer
    setupSimulationAlgorithm(std::move(fieldMap),
                             sequencer,
                             vm,
                             tGeometry,
                             barcodeSvc,
                             randomNumberSvc);
  } else {
    // create the shared field
    using CField = Acts::ConstantBField;
    CField fieldMap(*bField.second);
    // now setup of the simulation and append it to the sequencer
    setupSimulationAlgorithm(std::move(fieldMap),
                             sequencer,
                             vm,
                             tGeometry,
                             barcodeSvc,
                             randomNumberSvc);
  }
}
