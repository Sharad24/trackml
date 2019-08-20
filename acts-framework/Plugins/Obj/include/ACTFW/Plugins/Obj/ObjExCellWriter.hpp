// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <TTree.h>
#include <fstream>
#include <mutex>
#include <sstream>
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/Plugins/Obj/ObjHelper.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Extrapolation/ExtrapolationCell.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace FW {

namespace Obj {

  /// @class ExtrapolationCellWriter
  ///
  /// An obj based implementation to write out extrapolation steps.
  ///
  /// Safe to use from multiple writer threads, each event gets
  /// its own stream and thus is thread safe per call
  template <typename charge_t>
  class ObjExCellWriter
      : public FW::WriterT<std::vector<Acts::ExtrapolationCell<charge_t>>>
  {
  public:
    // The nested configuration struct
    struct Config
    {
    public:
      /// The Input collection
      std::string collection = "";

      /// Output directory
      std::string outputDir = "";

      /// Output scalor
      double outputScalor = 1.;

      /// Precision for out
      size_t outputPrecision = 4;

      /// transverse momentum cut
      double outputPtCut = 150. * Acts::units::_MeV;
      double outputMaxVr = 1000. * Acts::units::_mm;

      /// the size of the bezier segment
      double outputBezierSegment = 0. * Acts::units::_mm;
    };

    /// Constructor
    ///
    /// @param cfg is the configuration object
    /// @parm level is the output logging level
    ObjExCellWriter(const Config&        cfg,
                    Acts::Logging::Level level = Acts::Logging::INFO);

    /// End-of-run hook
    ProcessCode
    endRun() final override
    {
      return ProcessCode::SUCCESS;
    }

  protected:
    /// The protected writeT method, called by the WriterT base
    /// @param [in] ctx is the algorithm context for event consistency
    /// @param [in] ecells are the celss to be written out
    ProcessCode
    writeT(const FW::AlgorithmContext&                           ctx,
           const std::vector<Acts::ExtrapolationCell<charge_t>>& ecells)
        final override;

  private:
    Config m_cfg;  ///< the configuration class of this writer
  };

  template <typename charge_t>
  ObjExCellWriter<charge_t>::ObjExCellWriter(
      const ObjExCellWriter<charge_t>::Config& cfg,
      Acts::Logging::Level                     level)
    : FW::WriterT<std::vector<Acts::ExtrapolationCell<charge_t>>>(
          cfg.collection,
          "ObjExCellWriter",
          level)
    , m_cfg(cfg)
  {
    // Validate the configuration
    if (m_cfg.collection.empty()) {
      throw std::invalid_argument("Missing input collection");
    }
  }

  template <class T>
  FW::ProcessCode
  ObjExCellWriter<T>::writeT(
      const FW::AlgorithmContext&                    ctx,
      const std::vector<Acts::ExtrapolationCell<T>>& ecells)
  {

    // Initialize the vertex counter
    int vCounter = 0;

    // event suffix
    std::string suffix = "tracks.obj";

    std::string path
        = FW::perEventFilepath(m_cfg.outputDir, suffix, ctx.eventNumber);
    std::ofstream os(path, std::ofstream::out | std::ofstream::trunc);
    if (!os) {
      throw std::ios_base::failure("Could not open '" + path + "' to write");
    }
    os << std::setprecision(m_cfg.outputPrecision);

    // Loop over the cells
    for (auto& eCell : ecells) {

      // The ecc start paramters
      auto sPosition = eCell.startParameters->position();
      auto sMomentum = eCell.startParameters->momentum();
      // The momentum & vertex cut
      if (Acts::VectorHelpers::perp(sMomentum) < m_cfg.outputPtCut) continue;
      if (Acts::VectorHelpers::perp(sPosition) > m_cfg.outputMaxVr) continue;

      // Remember the first counter - for obj lines
      size_t fCounter = vCounter;

      // Loop over extrapolation steps - add bezier points
      auto lPosition  = eCell.startParameters->position();
      auto lDirection = sMomentum.normalized();

      for (auto& es : eCell.extrapolationSteps) {
        if (es.parameters) {
          // Take the step parameters
          const T& pars       = (*es.parameters);
          auto     tPosition  = pars.position();
          auto     tDirection = pars.momentum().normalized();
          // Don't write the start position another time
          if (tPosition == sPosition) continue;
          // Write the start parameters because
          // at least one other space point is here
          if (vCounter == fCounter) {
            // increase the vertex counter
            ++vCounter;
            // write the actual point
            os << "v " << m_cfg.outputScalor * sPosition.x() << ", "
               << m_cfg.outputScalor * sPosition.y() << ", "
               << m_cfg.outputScalor * sPosition.z() << " # initial point "
               << '\n';
          }
          if (m_cfg.outputBezierSegment > 0.) {
            // construct P1 and P2
            // we take the nominal distance divided by segments
            double nDist = (tPosition - lPosition).norm();
            // calculate the number of segments
            size_t segments = size_t(nDist / m_cfg.outputBezierSegment);
            if (segments > 1) {
              // rescale
              nDist /= (double)segments;
              Acts::Vector3D p1
                  = lPosition + m_cfg.outputBezierSegment * lDirection;
              Acts::Vector3D p2
                  = tPosition - m_cfg.outputBezierSegment * tDirection;
              // loop over the bezier segments
              for (size_t ib = 1; ib <= size_t(segments - 1); ib++) {
                double t = ib / (double)segments;
                auto   bPoint
                    = calculateBezierPoint(t, lPosition, p1, p2, tPosition);
                ++vCounter;
                // write the space point
                os << "v " << m_cfg.outputScalor * bPoint.x() << ", "
                   << m_cfg.outputScalor * bPoint.y() << ", "
                   << m_cfg.outputScalor * bPoint.z() << " # bezier point "
                   << '\n';
              }  // end of bezier segent writing
            }    // protection against only one segment
          }      // end of bezier condition
          // increase the counter, and write
          ++vCounter;
          os << "v " << m_cfg.outputScalor * tPosition.x() << ", "
             << m_cfg.outputScalor * tPosition.y() << ", "
             << m_cfg.outputScalor * tPosition.z() << " # excell point "
             << '\n';

          // set for the next bezier loop
          lPosition  = tPosition;
          lDirection = tDirection;
        }  // end of parameter check

      }  // end of extrapolation step loop

      // write out the line - only if we have at least two points created
      if ((vCounter - fCounter) > 2)
        for (size_t iv = fCounter + 1; iv < vCounter; ++iv)
          os << "l " << iv << " " << iv + 1 << '\n';

    }  // end of eCells loop
    os << '\n' << '\n';

    // return success
    return FW::ProcessCode::SUCCESS;
  }

}  // namespace Obj
}  // namespace FW
