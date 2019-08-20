// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>
#include <memory>

#include <Acts/Utilities/Units.hpp>
#include <boost/program_options.hpp>

#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/Common/CommonOptions.hpp"
#include "ACTFW/Common/OutputOptions.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Generators/ParticleSelector.hpp"
#include "ACTFW/Options/Pythia8Options.hpp"
#include "ACTFW/Plugins/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Plugins/Root/RootParticleWriter.hpp"
#include "ACTFW/Random/RandomNumbersOptions.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "ACTFW/Utilities/Paths.hpp"

using namespace Acts::units;
using namespace FW;

int
main(int argc, char* argv[])
{
  namespace po = boost::program_options;

  // define command line options and parse them
  po::options_description desc("Allowed options");
  Options::addCommonOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addPythia8Options(desc);
  Options::addOutputOptions(desc);
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  // print help if requested
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return EXIT_FAILURE;
  }

  Acts::Logging::Level logLevel  = Options::readLogLevel(vm);
  size_t               numEvents = Options::readNumberOfEvents(vm);

  Sequencer::Config sequencerCfg;
  Sequencer         sequencer(sequencerCfg);

  // basic services
  RandomNumbersSvc::Config rndCfg;
  rndCfg.seed  = 123;
  auto rnd     = std::make_shared<RandomNumbersSvc>(rndCfg);
  auto barcode = std::make_shared<BarcodeSvc>(
      BarcodeSvc::Config(), Acts::getDefaultLogger("BarcodeSvc", logLevel));

  // event generation w/ process guns
  EventGenerator::Config evgenCfg = Options::readPythia8Options(vm);
  evgenCfg.output                 = "generated_particles";
  evgenCfg.randomNumbers          = rnd;
  evgenCfg.barcodeSvc             = barcode;
  sequencer.addReaders({std::make_shared<EventGenerator>(evgenCfg, logLevel)});

  // event selection
  ParticleSelector::Config selectorCfg;
  selectorCfg.input  = evgenCfg.output;
  selectorCfg.output = "selected_particles";
  //  selectorCfg.absEtaMax = 2.0;
  //  selectorCfg.ptMin     = 0.5 * _GeV;
  selectorCfg.keepNeutral = false;  // retain only charged particles
  sequencer.appendEventAlgorithms(
      {std::make_shared<ParticleSelector>(selectorCfg, logLevel)});

  // different output modes
  std::string outputDir = vm["output-dir"].as<std::string>();

  if (vm["output-csv"].as<bool>()) {
    Csv::CsvParticleWriter::Config csvWriterCfg;
    csvWriterCfg.collection     = selectorCfg.output;
    csvWriterCfg.outputDir      = outputDir;
    csvWriterCfg.outputFileName = "particles.csv";
    sequencer.addWriters(
        {std::make_shared<Csv::CsvParticleWriter>(csvWriterCfg, logLevel)});
  }
  if (vm["output-root"].as<bool>()) {
    Root::RootParticleWriter::Config rootWriterCfg;
    rootWriterCfg.collection = selectorCfg.output;
    rootWriterCfg.filePath   = joinPaths(outputDir, "particles.root");
    rootWriterCfg.barcodeSvc = barcode;
    sequencer.addWriters(
        {std::make_shared<Root::RootParticleWriter>(rootWriterCfg, logLevel)});
  }

  return (sequencer.run(numEvents) == ProcessCode::SUCCESS) ? EXIT_SUCCESS
                                                            : EXIT_FAILURE;
}
