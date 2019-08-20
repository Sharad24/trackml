// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LayerArrayCreator.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Tools/LayerArrayCreator.hpp"
#include <cmath>
#include "Acts/Layers/Layer.hpp"
#include "Acts/Layers/NavigationLayer.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryObjectSorter.hpp"
#include "Acts/Utilities/GeometryStatics.hpp"

std::unique_ptr<const Acts::LayerArray>
Acts::LayerArrayCreator::layerArray(const LayerVector& layersInput,
                                    double             min,
                                    double             max,
                                    BinningType        bType,
                                    BinningValue       bValue) const
{
  ACTS_VERBOSE("Build LayerArray with " << layersInput.size()
                                        << " layers at input.");
  ACTS_VERBOSE("       min/max provided : " << min << " / " << max);
  ACTS_VERBOSE("       binning type     : " << bType);
  ACTS_VERBOSE("       binning value    : " << bValue);

  // create a local copy of the layer vector
  LayerVector layers(layersInput);

  // sort it accordingly to the binning value
  GeometryObjectSorterT<std::shared_ptr<const Layer>> layerSorter(bValue);
  std::sort(layers.begin(), layers.end(), layerSorter);
  // useful typedef
  using LayerOrderPosition = std::pair<std::shared_ptr<const Layer>, Vector3D>;
  // needed for all cases
  std::shared_ptr<const Layer>      layer      = nullptr;
  std::unique_ptr<const BinUtility> binUtility = nullptr;
  std::vector<LayerOrderPosition>   layerOrderVector;

  // switch the binning type
  switch (bType) {
  // equidistant binning - no navigation layers built - only equdistant layers
  case equidistant: {
    // loop over layers and put them in
    for (auto& layIter : layers) {
      ACTS_VERBOSE("equidistant : registering a Layer at binning position : "
                   << (layIter->binningPosition(bValue)));
      layerOrderVector.push_back(
          LayerOrderPosition(layIter, layIter->binningPosition(bValue)));
    }
    // create the binUitlity
    binUtility = std::make_unique<const BinUtility>(
        layers.size(), min, max, open, bValue);
    ACTS_VERBOSE("equidistant : created a BinUtility as " << *binUtility);
  } break;

  // arbitrary binning
  case arbitrary: {
    std::vector<float> boundaries;
    // initial step
    boundaries.push_back(min);
    double                       layerValue     = 0.;
    double                       layerThickness = 0.;
    std::shared_ptr<const Layer> navLayer       = nullptr;
    std::shared_ptr<const Layer> lastLayer      = nullptr;
    // loop over layers
    for (auto& layIter : layers) {
      // estimate the offset
      layerThickness = layIter->thickness();
      layerValue     = layIter->binningPositionValue(bValue);
      // register the new boundaries in the step vector
      boundaries.push_back(layerValue - 0.5 * layerThickness);
      boundaries.push_back(layerValue + 0.5 * layerThickness);
      // calculate the layer value for the offset
      double navigationValue = 0.5 * ((layerValue - 0.5 * layerThickness)
                                      + boundaries.at(boundaries.size() - 3));
      // if layers are attached to each other bail out - navigation will not
      // work anymore
      if (navigationValue == (layerValue - 0.5 * layerThickness)) {
        ACTS_ERROR("Layers are attached to each other at: "
                   << layerValue - 0.5 * layerThickness
                   << ", which corrupts "
                      "navigation. This should never happen. Please detach the "
                      "layers in your geometry description.");
      }
      // if layers are overlapping bail out
      if (navigationValue > (layerValue - 0.5 * layerThickness)) {
        ACTS_ERROR("Layers are overlapping at: "
                   << layerValue - 0.5 * layerThickness
                   << ". This should never happen. "
                      "Please check your geometry description.");
      }

      // create the navigation layer surface from the layer
      std::shared_ptr<const Surface> navLayerSurface = createNavigationSurface(
          *layIter, bValue, -std::abs(layerValue - navigationValue));
      ACTS_VERBOSE("arbitrary : creating a  NavigationLayer at "
                   << (navLayerSurface->binningPosition(bValue)).x()
                   << ", "
                   << (navLayerSurface->binningPosition(bValue)).y()
                   << ", "
                   << (navLayerSurface->binningPosition(bValue)).z());
      navLayer = NavigationLayer::create(std::move(navLayerSurface));
      // push the navigation layer in
      layerOrderVector.push_back(
          LayerOrderPosition(navLayer, navLayer->binningPosition(bValue)));

      // push the original layer in
      layerOrderVector.push_back(
          LayerOrderPosition(layIter, layIter->binningPosition(bValue)));
      ACTS_VERBOSE("arbitrary : registering MaterialLayer at  "
                   << (layIter->binningPosition(bValue)).x()
                   << ", "
                   << (layIter->binningPosition(bValue)).y()
                   << ", "
                   << (layIter->binningPosition(bValue)).z());
      // remember the last
      lastLayer = layIter;
    }
    // a final navigation layer
    // calculate the layer value for the offset
    double navigationValue = 0.5 * (boundaries.at(boundaries.size() - 1) + max);
    // create navigation layer only when necessary
    if (navigationValue != max) {
      // create the navigation layer surface from the layer
      std::shared_ptr<const Surface> navLayerSurface = createNavigationSurface(
          *lastLayer, bValue, navigationValue - layerValue);
      ACTS_VERBOSE("arbitrary : creating a  NavigationLayer at "
                   << (navLayerSurface->binningPosition(bValue)).x()
                   << ", "
                   << (navLayerSurface->binningPosition(bValue)).y()
                   << ", "
                   << (navLayerSurface->binningPosition(bValue)).z());
      navLayer = NavigationLayer::create(std::move(navLayerSurface));
      // push the navigation layer in
      layerOrderVector.push_back(
          LayerOrderPosition(navLayer, navLayer->binningPosition(bValue)));
    }
    // now close the boundaries
    boundaries.push_back(max);
    // some screen output
    ACTS_VERBOSE(layerOrderVector.size()
                 << " Layers (material + navigation) built. ");
    // create the BinUtility
    binUtility = std::make_unique<const BinUtility>(boundaries, open, bValue);
    ACTS_VERBOSE("arbitrary : created a BinUtility as " << *binUtility);

  } break;
  // default return nullptr
  default: {
    return nullptr;
  }
  }
  // return the binned array
  return std::make_unique<const BinnedArrayXD<LayerPtr>>(layerOrderVector,
                                                         std::move(binUtility));
}

std::shared_ptr<Acts::Surface>
Acts::LayerArrayCreator::createNavigationSurface(const Layer& layer,
                                                 BinningValue bValue,
                                                 double       offset) const
{
  // surface reference
  const Surface& layerSurface = layer.surfaceRepresentation();
  // translation to be applied
  Vector3D translation(0., 0., 0.);
  // switching he binnig values
  switch (bValue) {
  // case x
  case binX: {
    translation = Vector3D(offset, 0., 0.);
  } break;
  // case y
  case binY: {
    translation = Vector3D(0., offset, 0.);
  } break;
  // case z
  case binZ: {
    translation = Vector3D(0., 0., offset);
  } break;
  // case R
  case binR: {
    // binning in R and cylinder surface means something different
    if (layerSurface.type() == Surface::Cylinder) {
      break;
    }
    translation = Vector3D(offset, 0., 0.);
  } break;
  // do nothing for the default
  default: {
    ACTS_WARNING("Not yet implemented.");
  }
  }
  // navigation surface
  std::shared_ptr<Surface> navigationSurface;
  // for everything else than a cylinder it's a copy with shift
  if (layerSurface.type() != Surface::Cylinder) {
    // create a transform that does the shift
    const Transform3D* shift = new Transform3D(Translation3D(translation));
    navigationSurface        = layerSurface.clone(shift);
    // delete the shift again
    delete shift;
  } else {
    // get the bounds
    const CylinderBounds* cBounds
        = dynamic_cast<const CylinderBounds*>(&(layerSurface.bounds()));
    double navigationR = cBounds->r() + offset;
    double halflengthZ = cBounds->halflengthZ();
    // create the new layer surface
    std::shared_ptr<const Transform3D> navTrasform
        = (!layerSurface.transform().isApprox(s_idTransform))
        ? std::make_shared<const Transform3D>(layerSurface.transform())
        : nullptr;
    // new navigation layer
    navigationSurface = Surface::makeShared<CylinderSurface>(
        navTrasform, navigationR, halflengthZ);
  }
  return navigationSurface;
}
