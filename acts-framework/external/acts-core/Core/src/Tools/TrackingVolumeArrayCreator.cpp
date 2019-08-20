// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingVolumeArrayCreator.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Tools/TrackingVolumeArrayCreator.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryObjectSorter.hpp"
#include "Acts/Volumes/VolumeBounds.hpp"

std::shared_ptr<const Acts::TrackingVolumeArray>
Acts::TrackingVolumeArrayCreator::trackingVolumeArray(
    const TrackingVolumeVector& tVolumes,
    BinningValue                bValue) const
{
  // MSG_VERBOSE("Create VolumeArray of "<< tVolumes.size() << " TrackingVolumes
  // with binning in : " << binningValueNames[bValue] );
  // let's copy and sort
  TrackingVolumeVector volumes(tVolumes);
  // sort it accordingly to the binning value
  GeometryObjectSorterT<std::shared_ptr<const TrackingVolume>> volumeSorter(
      bValue);
  std::sort(volumes.begin(), volumes.end(), volumeSorter);

  // prepare what we need :
  // (1) arbitrary binning for volumes is fast enough
  std::vector<float> boundaries;
  boundaries.reserve(tVolumes.size() + 1);
  // (2) the vector needed for the BinnedArray
  std::vector<TrackingVolumeOrderPosition> tVolumesOrdered;

  // let's loop over the (sorted) volumes
  for (auto& tVolume : volumes) {
    // get the binning position
    Vector3D binningPosition = tVolume->binningPosition(bValue);
    double   binningBorder   = tVolume->volumeBounds().binningBorder(bValue);
    // get the center value according to the bin
    double value = tVolume->binningPositionValue(bValue);
    // for the first one take low and high boundary
    if (boundaries.empty()) {
      boundaries.push_back(value - binningBorder);
    }
    // always take the high boundary
    boundaries.push_back(value + binningBorder);
    // record the volume to be ordered
    tVolumesOrdered.push_back(
        TrackingVolumeOrderPosition(tVolume, binningPosition));
  }

  // now create the bin utility
  auto binUtility
      = std::make_unique<const BinUtility>(boundaries, open, bValue);

  // and return the newly created binned array
  return std::make_shared<const BinnedArrayXD<TrackingVolumePtr>>(
      tVolumesOrdered, std::move(binUtility));
}
