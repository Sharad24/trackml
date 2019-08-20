// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LayerArrayCreator.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#ifndef ACTS_TOOLS_TAKESMALLERBIGGER
#define ACTS_TOOLS_TAKESMALLERBIGGER
#define takeSmaller(current, test) current = current < test ? current : test
#define takeBigger(current, test) current  = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test)                           \
  takeSmaller(cSmallest, test);                                                \
  takeBigger(cBiggest, test)
#endif

#include <algorithm>
#include "Acts/Tools/ILayerArrayCreator.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

class Surface;
class Layer;

/// @class LayerArrayCreator
///  The LayerArrayCreator is a simple Tool that helps to construct
///  LayerArrays from std::vector of Acts::CylinderLayer, Acts::DiscLayer,
/// Acts::PlaneLayer.
///
///  It fills the gaps automatically with Acts::NavigationLayer to be processed
/// easily in the
///  Navigation of the Extrapolation process.
///

class LayerArrayCreator : public ILayerArrayCreator
{
public:
  /// Constructor
  ///
  /// @param logger logging instance
  LayerArrayCreator(std::unique_ptr<const Logger> logger
                    = getDefaultLogger("LayerArrayCreator", Logging::INFO))
    : m_logger(std::move(logger))
  {
  }

  /// Destructor
  ~LayerArrayCreator() override = default;

  /// LayerArrayCreator interface method
  ///
  /// @param layersInput are the layers to be moved into an array
  /// @param min is the minimum value for binning
  /// @param max is the maximum value for binning
  /// @param bType is the binning type
  /// @param bValue is the value in which the binning should be done
  ///
  /// @return unique pointer to a newly created LayerArray
  std::unique_ptr<const LayerArray>
  layerArray(const LayerVector& layersInput,
             double             min,
             double             max,
             BinningType        bType  = arbitrary,
             BinningValue       bValue = binX) const override;

  /// set logging instance
  void
  setLogger(std::unique_ptr<const Logger> logger)
  {
    m_logger = std::move(logger);
  }

private:
  /// Private access method to the logging instance
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private helper method for creating a surface for
  /// the NavigationLayer
  std::shared_ptr<Surface>
  createNavigationSurface(const Layer& layer,
                          BinningValue bValue,
                          double       offset) const;
};

}  // namespace