// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LayerCreator.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <boost/optional.hpp>
#include "Acts/Tools/SurfaceArrayCreator.hpp"
#include "Acts/Utilities/ApproachDescriptor.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

#ifndef ACTS_LAYERCREATOR_TAKESMALLERBIGGER
#define ACTS_LAYERCREATOR_TAKESMALLERBIGGER
#define takeSmaller(current, test) current = current < test ? current : test
#define takeBigger(current, test) current  = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test)                           \
  takeSmaller(cSmallest, test);                                                \
  takeBigger(cBiggest, test)
#endif

namespace Acts {

namespace Test {
  struct LayerCreatorFixture;
}
class Surface;
class Layer;
using MutableLayerPtr = std::shared_ptr<Layer>;

/// @class LayerCreator
///
/// The LayerCreator is able to build cylinde,r disc layers or plane layers from
/// detector elements
///
class LayerCreator
{
public:
  friend Acts::Test::LayerCreatorFixture;
  ///  @struct Config
  ///  Configuration for the LayerCreator
  ///  This is the nexted configuration struct for the LayerCreator class
  struct Config
  {
    /// surface array helper
    std::shared_ptr<const SurfaceArrayCreator> surfaceArrayCreator = nullptr;
    /// cylinder module z tolerance : it counts at same z, if ...
    double cylinderZtolerance{10.};
    /// cylinder module phi tolerance : it counts at same phi, if ...
    double cylinderPhiTolerance{0.1};

    // standard constructor
    Config() = default;
  };

  /// Constructor
  ///
  /// @param lcConfig is the configuration object
  /// @param logger logging instance
  LayerCreator(const Config&                 lcConfig,
               std::unique_ptr<const Logger> logger
               = getDefaultLogger("LayerCreator", Logging::INFO));

  /// Destructor
  ~LayerCreator() = default;

  /// returning a cylindrical layer
  ///
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// represented by this layer
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param binsPhi is number of bins the sensitive surfaces are ordered in phi
  /// @param binsZ is number of bins the sensitive surfaces are ordered in Z
  /// @param _protoLayer (optional) proto layer specifying the dimensions and
  /// envelopes
  /// @param transform is the (optional) transform of the layer
  /// @param ad possibility to hand over a specific ApproachDescriptor, which is
  /// needed for material mapping. Otherwise the default ApproachDescriptor will
  /// be taken used for this layer
  ///
  /// @return shared pointer to a newly created layer
  MutableLayerPtr
  cylinderLayer(std::vector<std::shared_ptr<const Surface>> surfaces,
                size_t                                      binsPhi,
                size_t                                      binsZ,
                boost::optional<ProtoLayer>         _protoLayer = boost::none,
                std::shared_ptr<const Transform3D>  transform   = nullptr,
                std::unique_ptr<ApproachDescriptor> ad = nullptr) const;

  /// returning a cylindrical layer
  ///
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// represented by this layer
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param bTypePhi binning type in phi (equidistant/arbitrary)
  /// @param bTypeZ binning type in z (equidistant/arbitrary)
  /// @param _protoLayer (optional) proto layer specifying the dimensions and
  /// envelopes
  /// @param transform is the (optional) transform of the layer
  /// @param ad possibility to hand over a specific ApproachDescriptor, which is
  /// needed for material mapping. Otherwise the default ApproachDescriptor will
  /// be taken used for this layer
  ///
  /// @return shared pointer to a newly created layer
  MutableLayerPtr
  cylinderLayer(std::vector<std::shared_ptr<const Surface>> surfaces,
                BinningType                                 bTypePhi,
                BinningType                                 bTypeZ,
                boost::optional<ProtoLayer>         _protoLayer = boost::none,
                std::shared_ptr<const Transform3D>  transform   = nullptr,
                std::unique_ptr<ApproachDescriptor> ad = nullptr) const;

  /// returning a disc layer
  ///
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// represented by this layer
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param binsR is number of bins the sensitive surfaces are ordered in R
  /// @param binsPhi is number of bins the sensitive surfaces are ordered in Phi
  /// @param transform is the (optional) transform of the layer
  /// @param _protoLayer (optional) proto layer specifying the dimensions and
  /// envelopes
  /// @param ad possibility to hand over a specific ApproachDescriptor, which is
  /// needed for material mapping. Otherwise the default ApproachDescriptor will
  /// be taken used for this layer
  ///
  /// @return shared pointer to a newly created layer
  MutableLayerPtr
  discLayer(std::vector<std::shared_ptr<const Surface>> surfaces,
            size_t                                      binsR,
            size_t                                      binsPhi,
            boost::optional<ProtoLayer>         _protoLayer = boost::none,
            std::shared_ptr<const Transform3D>  transform   = nullptr,
            std::unique_ptr<ApproachDescriptor> ad          = nullptr) const;

  /// returning a disc layer
  ///
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// represented by this layer
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param bTypeR binning type in r (equidistant/arbitrary)
  /// @param bTypePhi binning type in phi (equidistant/arbitrary)
  /// @param transform is the (optional) transform of the layer
  /// @param _protoLayer (optional) proto layer specifying the dimensions and
  /// envelopes
  /// @param ad possibility to hand over a specific ApproachDescriptor, which is
  /// needed for material mapping. Otherwise the default ApproachDescriptor will
  /// be taken used for this layer
  ///
  /// @return shared pointer to a newly created layer
  MutableLayerPtr
  discLayer(std::vector<std::shared_ptr<const Surface>> surfaces,
            BinningType                                 bTypeR,
            BinningType                                 bTypePhi,
            boost::optional<ProtoLayer>         _protoLayer = boost::none,
            std::shared_ptr<const Transform3D>  transform   = nullptr,
            std::unique_ptr<ApproachDescriptor> ad          = nullptr) const;

  /// returning a plane layer
  ///
  /// @param [in] surfaces is the vector of pointers to sensitive surfaces
  /// represented by this layer
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param [in] bins1 is the number of bins in the orthogonal direction to @p
  /// bValue
  /// @param [in] bins2 is the number of bins in the orthogonal direction to @p
  /// bValue
  /// @param [in] bValue Direction of the aligned surfaces
  /// @param [in] transform is the (optional) transform of the layer
  /// @param [in] _protoLayer (optional) proto layer specifying the dimensions
  /// and
  /// envelopes
  /// @param [in] ad possibility to hand over a specific ApproachDescriptor,
  /// which is needed for material mapping. Otherwise the default
  /// ApproachDescriptor will be taken used for this layer
  ///
  /// @return shared pointer to a newly created layer
  MutableLayerPtr
  planeLayer(std::vector<std::shared_ptr<const Surface>> surfaces,
             size_t                                      bins1,
             size_t                                      bins2,
             BinningValue                        bValue = BinningValue::binX,
             boost::optional<ProtoLayer>         _protoLayer = boost::none,
             std::shared_ptr<const Transform3D>  transform   = nullptr,
             std::unique_ptr<ApproachDescriptor> ad          = nullptr) const;

  /// Set the configuration object
  /// @param lcConfig is the configuration struct
  void
  setConfiguration(const Config& lcConfig);

  /// Access th configuration object
  Config
  getConfiguration() const;

  /// set logging instance
  /// @param newLogger the logger instance
  void
  setLogger(std::unique_ptr<const Logger> newLogger);

  // associate surfaces contained by this layer to this layer
  void
  associateSurfacesToLayer(Layer& layer) const;

private:
  /// Validates that all the sensitive surfaces are actually accessible through
  /// the binning
  ///
  /// @param surfGrid is the object grid from the surface array
  /// @para surfaces is the vector of sensitive surfaces
  bool
  checkBinning(const SurfaceArray& sArray) const;

  /// configuration object
  Config m_cfg;

  /// Private acces method to the logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// logging instance
  std::unique_ptr<const Logger> m_logger;
};

inline LayerCreator::Config
LayerCreator::getConfiguration() const
{
  return m_cfg;
}

}  // namespace Acts
