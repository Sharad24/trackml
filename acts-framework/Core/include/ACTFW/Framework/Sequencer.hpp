// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2016-05-11 Initial version
/// @date 2017-07-27 Clean up with simplified interfaces
/// @author Andreas Salzburger
/// @author Moritz Kiehn <msmk@cern.ch>

#ifndef ACTFW_SEQUENCER_H
#define ACTFW_SEQUENCER_H

#include <boost/optional.hpp>
#include <memory>
#include <string>
#include <tbb/task_scheduler_init.h>
#include <vector>
#include "ACTFW/Framework/IAlgorithm.hpp"
#include "ACTFW/Framework/IReader.hpp"
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/IWriter.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace FW {

/// @class  Sequencer
///
/// This is the backbone of the mini framework, it initializes all algorithms,
/// calls execute per event and deals with the event store */
///
class Sequencer
{
public:
  struct Config
  {
    /// job store logging level
    Acts::Logging::Level jobStoreLogLevel = Acts::Logging::INFO;
    /// event store logging level
    Acts::Logging::Level eventStoreLogLevel = Acts::Logging::INFO;
  };

  /// Constructor
  ///
  /// @param cfg is the configuration object
  Sequencer(const Config&                       cfg,
            std::unique_ptr<const Acts::Logger> logger
            = Acts::getDefaultLogger("Sequencer", Acts::Logging::INFO));

  /// Add services
  ///
  /// @param services is the vector of services to be added
  ProcessCode
  addServices(std::vector<std::shared_ptr<IService>> services);

  /// Add algorithms for reading
  ///
  /// @param readers is the vector of reader algorithms to be added
  ProcessCode
  addReaders(std::vector<std::shared_ptr<IReader>> readers);

  /// Add algorithms for writing
  ///
  /// @param writers is the vector of writer algorithms to be added
  ProcessCode
  addWriters(std::vector<std::shared_ptr<IWriter>> writers);

  /// Prepend algorithms
  ///
  /// @param algorithms is the vector of algorithms to be prepended
  ProcessCode
  prependEventAlgorithms(std::vector<std::shared_ptr<IAlgorithm>> algorithms);

  /// Append algorithms
  ///
  /// @param algorithms is the vector of algorithms to be appended
  ProcessCode
  appendEventAlgorithms(std::vector<std::shared_ptr<IAlgorithm>> algorithms);

  /// Run the event loop over the given number of events.
  ///
  /// @param events (optional) Number of events to process
  /// @note This parameter is optional when input is read from a file. In this
  /// scenario, leaving it unset will process events until the end of the file,
  /// and setting it will put an upper bound on the number of events to be
  /// processed.
  /// @param skip Number of events to skip before processing
  ///
  /// This will run all configured algorithms for each event, potentially in
  /// parallel, then invoke the endRun hook of writers and services.
  ProcessCode
  run(boost::optional<size_t> events, size_t skip = 0);

private:
  std::vector<std::shared_ptr<IService>>   m_services;
  std::vector<std::shared_ptr<IReader>>    m_readers;
  std::vector<std::shared_ptr<IWriter>>    m_writers;
  std::vector<std::shared_ptr<IAlgorithm>> m_algorithms;
  Config                                   m_cfg;
  std::unique_ptr<const Acts::Logger>      m_logger;
  tbb::task_scheduler_init                 m_tbb_init;

  const Acts::Logger&
  logger() const
  {
    return *m_logger;
  }
};

}  // namespace FW

#endif  // ACTFW_SEQUENCER_H
