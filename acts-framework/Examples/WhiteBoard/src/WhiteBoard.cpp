// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Framework/WhiteBoard.hpp"
#include <boost/program_options.hpp>
#include <cstdlib>
#include <memory>
#include "ACTFW/Common/CommonOptions.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "WhiteBoardAlgorithm.hpp"

namespace po = boost::program_options;

/// Main read evgen executable
///
/// @param argc The argument count
/// @param argv The argument list
int
main(int argc, char* argv[])
{
  // Declare the supported program options.
  po::options_description desc("Allowed options");
  // Add the standard/common options
  FW::Options::addCommonOptions<po::options_description>(desc);
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
  // Read the common options
  auto nEvents  = FW::Options::readNumberOfEvents<po::variables_map>(vm);
  auto logLevel = FW::Options::readLogLevel<po::variables_map>(vm);

  // Create an algorithm that writes to the event store
  FWE::WhiteBoardAlgorithm::Config wBoardConfigWrite;
  wBoardConfigWrite.outputClassOneCollection = "ClassOneCollection";
  wBoardConfigWrite.outputClassTwoCollection = "ClassTwoCollection";
  auto wBoardWrite
      = std::make_shared<FWE::WhiteBoardAlgorithm>(wBoardConfigWrite, logLevel);

  // Create an algorithm that reads from the event store
  FWE::WhiteBoardAlgorithm::Config wBoardConfigRead;
  wBoardConfigRead.inputClassOneCollection = "ClassOneCollection";
  wBoardConfigRead.inputClassTwoCollection = "ClassTwoCollection";
  auto wBoardRead
      = std::make_shared<FWE::WhiteBoardAlgorithm>(wBoardConfigRead, logLevel);

  // Create the event loop
  FW::Sequencer::Config seqConfig;
  seqConfig.eventStoreLogLevel = logLevel;
  FW::Sequencer sequencer(seqConfig);
  sequencer.appendEventAlgorithms({wBoardWrite, wBoardRead});

  // Run the event loop
  sequencer.run(nEvents);

  return 0;
}
