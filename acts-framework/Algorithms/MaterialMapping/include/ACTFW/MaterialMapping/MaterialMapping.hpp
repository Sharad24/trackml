// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTFW_ALGORITHMS_MATERIALMAPPING_MATERIALMAPPING_H
#define ACTFW_ALGORITHMS_MATERIALMAPPING_MATERIALMAPPING_H

#include <climits>
#include <memory>

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Readers/IReaderT.hpp"
#include "ACTFW/Writers/IWriterT.hpp"
#include "Acts/Layers/Layer.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Plugins/MaterialMapping/MaterialMapper.hpp"
#include "Acts/Plugins/MaterialMapping/MaterialTrack.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace FW {
class WhiteBoard;
}

namespace Acts {
class TrackingGeometry;
}

namespace FW {

/// @class MaterialMapping
///
/// @brief Initiates material mapping
///
/// The MaterialMapping reads in the MaterialTrack with a dedicated
/// reader and uses the material mapper to project the material onto
/// the tracking geometry
///
/// In a final step, the material maps are written out for further usage

class MaterialMapping : public FW::BareAlgorithm
{
public:
  /// @class nested Config class
  /// of the MaterialMapping algorithm
  struct Config
  {
  public:
    /// The reader to read in the MaterialTrack entities
    std::shared_ptr<FW::IReaderT<Acts::MaterialTrack>> materialTrackReader
        = nullptr;
    /// The ACTS material mapper
    std::shared_ptr<Acts::MaterialMapper> materialMapper = nullptr;
    /// The validation writer of the material
    std::shared_ptr<FW::IWriterT<Acts::MaterialTrack>> materialTrackWriter
        = nullptr;
    /// The writer of the material
    std::shared_ptr<FW::IWriterT<Acts::IndexedSurfaceMaterial>>
        indexedMaterialWriter = nullptr;
    /// The TrackingGeometry to be mapped on
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
    /// mapping conditions
    size_t maximumTrackRecords = std::numeric_limits<size_t>::infinity();
  };

  /// Constructor
  MaterialMapping(const Config&        cfg,
                  Acts::Logging::Level level = Acts::Logging::INFO);

  /// Framework execute method
  FW::ProcessCode
  execute(FW::AlgorithmContext context) const final override;

private:
  Config m_cfg;
};

}  // namespace FW

#endif  // ACTFW_ALGORITHMS_MATERIALMAPPING_MATERIALMAPPING_H
