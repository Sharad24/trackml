// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>
#include <cstdlib>
#include <memory>
#include "ACTFW/Common/CommonOptions.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Random/RandomNumbersOptions.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "RandomNumbersAlgorithm.hpp"

namespace po = boost::program_options;

// @brief An example to use the random number service
//
// @param argc the number of arguments provided
// @param argv the argument list
int
main(int argc, char* argv[])
{
  // Declare the supported program options.
  po::options_description desc("Allowed options");
  // Add the common options
  FW::Options::addCommonOptions<po::options_description>(desc);
  // Add the random number options
  FW::Options::addRandomNumbersOptions<po::options_description>(desc);
  // Map to store the given program options
  po::variables_map vm;
  // Get all options from contain line and store it into the map
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  // print help if requested
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  // Read the common options : number of events and log level
  auto nEvents  = FW::Options::readNumberOfEvents<po::variables_map>(vm);
  auto logLevel = FW::Options::readLogLevel<po::variables_map>(vm);

  // Create the random number engine
  auto randomNumbersCfg
      = FW::Options::readRandomNumbersConfig<po::variables_map>(vm);
  auto randomNumbers = std::make_shared<FW::RandomNumbersSvc>(randomNumbersCfg);

  // Create the config object for the hello world algorithm
  FWE::RandomNumbersAlgorithm::Config rNumbersConfig;
  rNumbersConfig.randomNumbers     = randomNumbers;
  rNumbersConfig.gaussParameters   = {{0., 1.}};
  rNumbersConfig.uniformParameters = {{0., 1.}};
  rNumbersConfig.landauParameters  = {{1., 7.}};
  rNumbersConfig.gammaParameters   = {{1., 1.}};
  rNumbersConfig.drawsPerEvent     = 5000;

  // And now the hello world algorithm
  std::shared_ptr<FW::IAlgorithm> rNumbers(
      new FWE::RandomNumbersAlgorithm(rNumbersConfig, logLevel));
  // Create the config object for the sequencer
  FW::Sequencer::Config seqConfig;
  // Create the sequencer and run the example
  FW::Sequencer sequencer(seqConfig);
  sequencer.addServices({randomNumbers});
  sequencer.appendEventAlgorithms({rNumbers});
  sequencer.run(nEvents);

  return 0;
}
