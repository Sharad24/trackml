// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <exception>

#include <tbb/tbb.h>

#include <TROOT.h>

#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

FW::Sequencer::Sequencer(const Sequencer::Config&            cfg,
                         std::unique_ptr<const Acts::Logger> logger)
  : m_cfg(cfg)
  , m_logger(std::move(logger))
  , m_tbb_init(tbb::task_scheduler_init::deferred)
{
  ROOT::EnableThreadSafety();

  const char* num_threads_str = getenv("ACTSFW_NUM_THREADS");
  int         num_threads;
  if (num_threads_str) {
    num_threads = std::stoi(num_threads_str);
  } else {
    num_threads = tbb::task_scheduler_init::automatic;
  }
  m_tbb_init.initialize(num_threads);
}

FW::ProcessCode
FW::Sequencer::addServices(std::vector<std::shared_ptr<FW::IService>> services)
{
  for (auto& svc : services) {
    if (!svc) {
      ACTS_FATAL("Trying to add empty service to sequencer");
      return ProcessCode::ABORT;
    }
    m_services.push_back(std::move(svc));
    ACTS_INFO("Added service " << m_services.back()->name());
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::Sequencer::addReaders(std::vector<std::shared_ptr<FW::IReader>> readers)
{
  for (auto& rdr : readers) {
    if (!rdr) {
      ACTS_FATAL("Trying to add empty reader to sequencer");
      return ProcessCode::ABORT;
    }
    m_readers.push_back(std::move(rdr));
    ACTS_INFO("Added reader " << m_readers.back()->name());
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::Sequencer::addWriters(std::vector<std::shared_ptr<FW::IWriter>> writers)
{
  for (auto& wrt : writers) {
    if (!wrt) {
      ACTS_FATAL("Trying to add empty writer to sequencer");
      return ProcessCode::ABORT;
    }
    m_writers.push_back(std::move(wrt));
    ACTS_INFO("Added writer " << m_writers.back()->name());
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::Sequencer::prependEventAlgorithms(
    std::vector<std::shared_ptr<FW::IAlgorithm>> algorithms)
{
  for (auto& alg : algorithms) {
    if (!alg) {
      ACTS_FATAL("Trying to prepend empty algorithm");
      return ProcessCode::ABORT;
    }
    m_algorithms.insert(m_algorithms.begin(), std::move(alg));
    ACTS_INFO("Prepended algorithm " << m_algorithms.front()->name());
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::Sequencer::appendEventAlgorithms(
    std::vector<std::shared_ptr<FW::IAlgorithm>> algorithms)
{
  for (auto& alg : algorithms) {
    if (!alg) {
      ACTS_FATAL("Trying to append empty algorithm.");
      return ProcessCode::ABORT;
    }
    m_algorithms.push_back(std::move(alg));
    ACTS_INFO("Appended algorithm " << m_algorithms.back()->name());
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::Sequencer::run(boost::optional<size_t> events, size_t skip)
{
  // Print some introduction
  ACTS_INFO("Starting event loop for");
  ACTS_INFO("  " << m_services.size() << " services");
  ACTS_INFO("  " << m_readers.size() << " readers");
  ACTS_INFO("  " << m_writers.size() << " writers");
  ACTS_INFO("  " << m_algorithms.size() << " algorithms");

  // The number of events to be processed
  size_t numEvents = 0;

  // There are two possibilities how the event loop can be steered
  // 1) By the number of given events
  // 2) By the number of events given by the readers
  // Calulate minimum and maximum of events to be read in
  auto min = std::min_element(
      m_readers.begin(), m_readers.end(), [](const auto& a, const auto& b) {
        return (a->numEvents() < b->numEvents());
      });
  // Check if number of events is given by the reader(s)
  if (min == m_readers.end()) {
    // 1) In case there are no readers, no event should be skipped
    if (skip != 0) {
      ACTS_ERROR(
          "Number of skipped events given although no readers present. Abort");
      return ProcessCode::ABORT;
    }
    // Number of events is not given by readers, in this case the parameter
    // 'events' must be specified - Abort, if this is not the case
    if (!events) {
      ACTS_ERROR("Number of events not specified. Abort");
      return ProcessCode::ABORT;
    }
    // 'events' is specified, set 'numEvents'
    numEvents = *events;
  } else {
    // 2) Number of events given by reader(s)
    numEvents = ((*min)->numEvents());
    // Check if the number of skipped events is smaller then the overall number
    // if events
    if (skip > numEvents) {
      ACTS_ERROR("Number of events to be skipped > than total number of "
                 "events. Abort");
      return ProcessCode::ABORT;
    }
    // The total number of events is the maximum number of events minus the
    // number of skipped evebts
    numEvents -= skip;
    // Check if user wants to process less events than given by the reader
    if (events && (*events) < numEvents) numEvents = *events;
  }

  // Execute the event loop
  ACTS_INFO("Run the event loop");
  tbb::parallel_for(
      tbb::blocked_range<size_t>(skip, numEvents + skip),
      [&](const tbb::blocked_range<size_t>& r) {
        for (size_t event = r.begin(); event != r.end(); ++event) {
          ACTS_INFO("start event " << event);

          // Setup the event and algorithm context
          WhiteBoard eventStore(Acts::getDefaultLogger(
              "EventStore#" + std::to_string(event), m_cfg.eventStoreLogLevel));
          size_t ialg = 0;

          // read everything in
          for (auto& rdr : m_readers) {
            if (rdr->read({ialg++, event, eventStore}) != ProcessCode::SUCCESS)
              throw std::runtime_error("Failed to read input data");
          }
          // process all algorithms
          for (auto& alg : m_algorithms) {
            if (alg->execute({ialg++, event, eventStore})
                != ProcessCode::SUCCESS)
              throw std::runtime_error("Failed to process event data");
          }
          // write out results
          for (auto& wrt : m_writers) {
            if (wrt->write({ialg++, event, eventStore}) != ProcessCode::SUCCESS)
              throw std::runtime_error("Failed to write output data");
          }

          ACTS_INFO("event " << event << " done");
        }
      });

  // Call endRun() for writers and services
  ACTS_INFO("Running end-of-run hooks of writers and services");
  for (auto& wrt : m_writers)
    if (wrt->endRun() != ProcessCode::SUCCESS) return ProcessCode::ABORT;
  for (auto& svc : m_services)
    if (svc->endRun() != ProcessCode::SUCCESS) return ProcessCode::ABORT;
  return ProcessCode::SUCCESS;
}
