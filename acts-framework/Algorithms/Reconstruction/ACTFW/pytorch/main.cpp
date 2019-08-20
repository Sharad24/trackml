// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "main.hpp"
#include <torch/torch.h>
#include <torch/script.h>
#include <iostream>

#include <stdexcept>

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/GeometryID.hpp>

#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

FW::PytorchReconstruction::PytorchReconstruction(
    const FW::PytorchReconstruction::Config& cfg,
    Acts::Logging::Level                            logLevel)
  : BareAlgorithm("PytorchReconstruction", logLevel), m_cfg(cfg)
{
  if (m_cfg.spacePointCollection.empty()) {
    throw std::invalid_argument("Missing input space points collection");
  }
  std::vector<torch::jit::IValue> inputs;
  inputs.push_back(torch::ones({2, 8}));

  // initialize script Module. The model should've been saved through torchscript.
  torch::jit::script::Module module;
  module = torch::jit::load("/tmp/tmp.pb");

  at::Tensor output = module.forward(inputs).toTensor();
  ACTS_INFO(output);

  // Load model
  // module = torch::jit::load("/tmp/tmp.pb");
  // ACTS_INFO("Loaded pretrained model");
}

FW::ProcessCode
FW::PytorchReconstruction::execute(FW::AlgorithmContext ctx) const
{
  ACTS_INFO("empty reconstruction on event " << ctx.eventNumber);

  const DetectorData<geo_id_value, Acts::Vector3D>* spacePoints;
  if (ctx.eventStore.get(m_cfg.spacePointCollection, spacePoints)
      == FW::ProcessCode::ABORT) {
    ACTS_WARNING("missing space point input");
    return FW::ProcessCode::ABORT;
  }
  // Create input data to pass into the network
  std::vector<torch::jit::IValue> inputs;
  inputs.push_back(torch::ones({2, 8}));
  ACTS_INFO("Initialized Input Data");
  //
  // // initialize script Module. The model should've been saved through torchscript.
  torch::jit::script::Module module2 = module;
  //
  // // Load model
  // module2 = torch::jit::load("/tmp/tmp.pb");

  // Get outputs and convert to tensor
  at::Tensor output = module2.forward(inputs).toTensor();
  ACTS_INFO(output);

  return ProcessCode::SUCCESS;
}
