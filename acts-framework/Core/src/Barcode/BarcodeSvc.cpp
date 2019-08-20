// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//
//  BarcodeSvc.cpp
//  ACTFW
//
//  Created by Andreas Salzburger on 17/05/16.
//
//
#include "ACTFW/Barcode/BarcodeSvc.hpp"

FW::BarcodeSvc::BarcodeSvc(const FW::BarcodeSvc::Config&       cfg,
                           std::unique_ptr<const Acts::Logger> logger)
  : m_cfg(cfg), m_logger(std::move(logger))
{
}

std::string
FW::BarcodeSvc::name() const
{
  return "BarcodeSvc";
}
