// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/Common/CommonOptions.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Generators/EventGenerator.hpp"
#include "ACTFW/Options/ParticleGunOptions.hpp"
#include "ACTFW/Options/Pythia8Options.hpp"
#include "ACTFW/Plugins/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Plugins/Root/RootParticleWriter.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "ACTFW/Utilities/Paths.hpp"

template <typename vmap_t>
void
setupEvgenInput(vmap_t&                               vm,
                FW::Sequencer&                        sequencer,
                std::shared_ptr<FW::BarcodeSvc>       barcodeSvc,
                std::shared_ptr<FW::RandomNumbersSvc> randomNumberSvc)
{
  // Read the standard options
  auto logLevel = FW::Options::readLogLevel(vm);
  auto nEvents  = FW::Options::readNumberOfEvents(vm);

  // Add requeste event generator
  auto evgenInput = vm["evg-input-type"].template as<std::string>();
  if (evgenInput == "gun") {
    auto evgCfg          = FW::Options::readParticleGunOptions(vm);
    evgCfg.output        = "particles";
    evgCfg.randomNumbers = randomNumberSvc;
    evgCfg.barcodeSvc    = barcodeSvc;
    sequencer.addReaders(
        {std::make_shared<FW::EventGenerator>(evgCfg, logLevel)});

  } else if (evgenInput == "pythia8") {
    auto evgCfg          = FW::Options::readPythia8Options(vm);
    evgCfg.output        = "particles";
    evgCfg.randomNumbers = randomNumberSvc;
    evgCfg.barcodeSvc    = barcodeSvc;
    sequencer.addReaders(
        {std::make_shared<FW::EventGenerator>(evgCfg, logLevel)});

  } else {
    throw std::runtime_error("unknown event generator input: " + evgenInput);
  }

  // Output directory
  std::string outputDir = vm["output-dir"].template as<std::string>();

  // Write particles as CSV files
  if (vm["output-csv"].template as<bool>()) {
    FW::Csv::CsvParticleWriter::Config pWriterCsvConfig;
    pWriterCsvConfig.collection     = "particles";
    pWriterCsvConfig.outputDir      = outputDir;
    pWriterCsvConfig.outputFileName = "particles.csv";
    sequencer.addWriters(
        {std::make_shared<FW::Csv::CsvParticleWriter>(pWriterCsvConfig)});
  }

  // Write particles as ROOT file
  if (vm["output-root"].template as<bool>()) {
    // Write particles as ROOT TTree
    FW::Root::RootParticleWriter::Config pWriterRootConfig;
    pWriterRootConfig.collection = "particles";
    pWriterRootConfig.filePath   = FW::joinPaths(outputDir, "particles.root");
    pWriterRootConfig.treeName   = "particles";
    pWriterRootConfig.barcodeSvc = barcodeSvc;
    sequencer.addWriters(
        {std::make_shared<FW::Root::RootParticleWriter>(pWriterRootConfig)});
  }
}
