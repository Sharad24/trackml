// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <vector>

namespace FW {

/// Data containers designed to fit around the GeometryID structure
///
/// internal map structure is
/// { volume : layer : module , data }
template <typename data_t>
using ModuleData = std::vector<data_t>;
template <typename identifier_t, typename data_t>
using LayerData = std::map<identifier_t, ModuleData<data_t>>;
template <typename identifier_t, typename data_t>
using VolumeData = std::map<identifier_t, LayerData<identifier_t, data_t>>;
template <typename identifier_t, typename data_t>
using DetectorData = std::map<identifier_t, VolumeData<identifier_t, data_t>>;

namespace Data {

  /// @brief Insert (& create container if necessary)
  ///
  /// @tparam identifier_t Type of the identifier key
  /// @tparam data_t Type of the object to be filled in
  ///
  /// @param dData The detector data map
  /// @param volumeKey The identifier for the detector volume
  /// @param layerKey The identifier for the layer
  /// @param moduleKey The identifier for the module
  /// @param obj The data object to be stored
  template <typename identifier_t, typename data_t>
  void
  insert(DetectorData<identifier_t, data_t>& dData,
         identifier_t volumeKey,
         identifier_t layerKey,
         identifier_t moduleKey,
         data_t       obj)
  {
    // find if the volume has an entry
    auto volumeData = dData.find(volumeKey);
    if (volumeData == dData.end()) {
      // insert at the volumeKey
      dData[volumeKey] = FW::VolumeData<identifier_t, data_t>();
      volumeData       = dData.find(volumeKey);
    }
    // find the layer data
    auto layerData = (volumeData->second).find(layerKey);
    if (layerData == (volumeData->second).end()) {
      // insert a layer key for this
      (volumeData->second)[layerKey] = FW::LayerData<identifier_t, data_t>();
      layerData                      = (volumeData->second).find(layerKey);
    }
    // find the module data
    auto moduleData = (layerData->second).find(moduleKey);
    if (moduleData == (layerData->second).end()) {
      // insert the module for this
      (layerData->second)[moduleKey] = FW::ModuleData<data_t>();
      moduleData                     = (layerData->second).find(moduleKey);
    }
    // and now push back
    (moduleData->second).push_back(std::move(obj));
  };

  /// @brief Read and return the module data
  ///
  /// @tparam identifier_t Type of the identifier key
  /// @tparam data_t Type of the object to be filled in
  ///
  /// @param dData The detector data map
  /// @param volumeKey The identifier for the detector volume
  /// @param layerKey The identifier for the layer
  /// @param moduleKey The identifier for the module
  template <typename identifier_t, typename data_t>
  const ModuleData<identifier_t>*
  read(DetectorData<identifier_t, data_t>& dData,
       identifier_t volumeKey,
       identifier_t layerKey,
       identifier_t moduleKey)
  {
    // find if the volume has an entry
    auto volumeData = dData.find(volumeKey);
    if (volumeData == dData.end()) return nullptr;
    // find the layer data
    auto layerData = (volumeData->second).find(layerKey);
    if (layerData == (volumeData->second).end()) return nullptr;
    // find the module data
    auto moduleData = (layerData->second).find(moduleKey);
    if (moduleData == (layerData->second).end()) return nullptr;
    // and now return as a pointer
    return (&(moduleData->second));
  };

  /// @brief Read and return the layer data
  ///
  /// @tparam identifier_t Type of the identifier key
  /// @tparam data_t Type of the object to be filled in
  ///
  /// @param dData The detector data map
  /// @param volumeKey The identifier for the detector volume
  /// @param layerKey The identifier for the layer
  /// @param moduleKey The identifier for the module
  template <typename identifier_t, typename data_t>
  const LayerData<identifier_t, data_t>*
  read(DetectorData<identifier_t, data_t>& dData,
       identifier_t volumeKey,
       identifier_t layerKey)
  {
    // find if the volume has an entry
    auto volumeData = dData.find(volumeKey);
    if (volumeData == dData.end()) return nullptr;
    // find the layer data
    auto layerData = (volumeData->second).find(layerKey);
    if (layerData == (volumeData->second).end()) return nullptr;
    // and now return as a pointer
    return (&(layerData->second));
  };

}  // namespace Data
}  // namespace FW
