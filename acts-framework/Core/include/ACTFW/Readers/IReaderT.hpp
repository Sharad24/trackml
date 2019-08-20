// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// IReaderT.h
///////////////////////////////////////////////////////////////////
#pragma once

#include <string>
#include <vector>
#include "ACTFW/Framework/AlgorithmContext.hpp"
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace FW {

/// @class IWriterT
///
/// Interface to read in an templated object
template <typename T>
class IReaderT : public IService
{
public:
  /// Virtual destructor
  virtual ~IReaderT() = default;

  /// The Read interface method
  /// @param object is the thing to be read back in
  /// @param skip The events to be skipped
  /// @param context The algorithm context in case it is needed (e.g. random
  /// number creation)
  /// @return is a ProcessCode indicating return/failure
  virtual ProcessCode
  read(T&                          object,
       size_t                      skip    = 0,
       const FW::AlgorithmContext* context = nullptr)
      = 0;
};
}
