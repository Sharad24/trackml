// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename simulator_t, typename event_collection_t, typename hit_t>
FW::FatrasAlgorithm<simulator_t, event_collection_t, hit_t>::FatrasAlgorithm(
    const Config&        cfg,
    Acts::Logging::Level loglevel)
  : FW::BareAlgorithm("FatrasAlgorithm", loglevel), m_cfg(cfg)
{
}

template <typename simulator_t, typename event_collection_t, typename hit_t>
FW::ProcessCode
FW::FatrasAlgorithm<simulator_t, event_collection_t, hit_t>::execute(
    const FW::AlgorithmContext context) const
{

  // Create an algorithm local random number generator
  RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator(context);

  // Read Particles from input collection
  const event_collection_t* inputEvent = nullptr;
  if (context.eventStore.get(m_cfg.inputEventCollection, inputEvent)
      == FW::ProcessCode::ABORT)
    return FW::ProcessCode::ABORT;

  ACTS_DEBUG("Read collection '" << m_cfg.inputEventCollection << "' with "
                                 << inputEvent->size()
                                 << " vertices");

  // Output: simulated particles attached to their process vertices
  // we start with a copy of the current event
  event_collection_t simulatedEvent(*inputEvent);

  // Create the hit collection
  HitCollection<hit_t> simulatedHits;

  // The simulation call
  m_cfg.simulator(rng, simulatedEvent, simulatedHits);

  // Write simulated data to the event store
  // - the simulated particles
  if (context.eventStore.add(m_cfg.simulatedEventCollection,
                             std::move(simulatedEvent))
      == FW::ProcessCode::ABORT) {
    return FW::ProcessCode::ABORT;
  }

  // The simulated hits
  if (context.eventStore.add(m_cfg.simulatedHitCollection,
                             std::move(simulatedHits.hits))
      == FW::ProcessCode::ABORT) {
    return FW::ProcessCode::ABORT;
  }
  return FW::ProcessCode::SUCCESS;
}
