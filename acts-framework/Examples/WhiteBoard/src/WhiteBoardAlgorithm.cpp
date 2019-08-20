// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "WhiteBoardAlgorithm.hpp"
#include <iostream>
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "DataClassOne.hpp"
#include "DataClassTwo.hpp"

FWE::WhiteBoardAlgorithm::WhiteBoardAlgorithm(const Config&        cfg,
                                              Acts::Logging::Level level)
  : FW::BareAlgorithm("WhiteBoardAlgorithm", level), m_cfg(cfg)
{
}

FW::ProcessCode
FWE::WhiteBoardAlgorithm::execute(FW::AlgorithmContext ctx) const
{
  // -------- Reading -----------------------
  // Reading Class One
  if (!m_cfg.inputClassOneCollection.empty()) {
    ACTS_INFO("Reading ClassOneCollection " << m_cfg.inputClassOneCollection);
    // read in the collection
    const FWE::DataClassOneCollection* dcoCollIn = nullptr;
    // write to the EventStore
    if (ctx.eventStore.get(m_cfg.inputClassOneCollection, dcoCollIn)
        == FW::ProcessCode::ABORT)
      return FW::ProcessCode::ABORT;
    // screen output
    ACTS_VERBOSE("Read DataClassOneCollection with size " << dcoCollIn->size());
    for (auto& idco : (*dcoCollIn))
      ACTS_VERBOSE("Read in  DataClassOne object as " << idco.data());
  }

  // Reading Class Two
  if (!m_cfg.inputClassTwoCollection.empty()) {
    ACTS_INFO("Reading ClassTwoCollection " << m_cfg.inputClassTwoCollection);
    // read in the collection
    const FWE::DataClassTwoCollection* dctCollIn = nullptr;
    // write to the EventStore
    if (ctx.eventStore.get(m_cfg.inputClassTwoCollection, dctCollIn)
        == FW::ProcessCode::ABORT)
      return FW::ProcessCode::ABORT;
    // screen output
    ACTS_VERBOSE("Read DataClassTwoCollection with size " << dctCollIn->size());
    for (auto& idct : (*dctCollIn))
      ACTS_VERBOSE("Read in  DataClassTwo object as " << idct.data());
  }

  // ---------- Writing -----------------------
  // Writing Class One
  if (!m_cfg.outputClassOneCollection.empty()) {
    ACTS_INFO("Writing ClassOneCollection " << m_cfg.outputClassOneCollection);
    // create a new collection
    DataClassOneCollection dcoCollOut = {{"One", ctx.eventNumber}};
    ACTS_VERBOSE("Written out DataClassOne object as "
                 << dcoCollOut.back().data());
    // write to the EventStore
    if (ctx.eventStore.add(m_cfg.outputClassOneCollection,
                           std::move(dcoCollOut))
        == FW::ProcessCode::ABORT)
      return FW::ProcessCode::ABORT;
  }

  // Writing Class Two
  if (!m_cfg.outputClassTwoCollection.empty()) {
    ACTS_INFO("Writing ClassTwoCollection " << m_cfg.outputClassTwoCollection);
    // create a new collection
    DataClassTwoCollection dctCollOut = {{"Two", double(ctx.eventNumber)}};
    ACTS_VERBOSE("Written out DataClassTwo object as "
                 << dctCollOut.back().data());
    // write to the EventStore
    if (ctx.eventStore.add(m_cfg.outputClassTwoCollection,
                           std::move(dctCollOut))
        == FW::ProcessCode::ABORT)
      return FW::ProcessCode::ABORT;
  }
  // Return with success
  return FW::ProcessCode::SUCCESS;
}
