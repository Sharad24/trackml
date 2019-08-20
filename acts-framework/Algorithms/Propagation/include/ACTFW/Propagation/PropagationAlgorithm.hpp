// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <limits>
#include <memory>
#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Random/RandomNumberDistributions.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"

/// @brief this test algorithm performs test propagation
/// within the Acts::Propagator
///
/// If the propagator is equipped appropriately, it can
/// also be used to test the Extrapolator within the geomtetry
///
/// @tparam propagator_t Type of the Propagator to be tested
namespace FW {

template <typename propagator_t>
class PropagationAlgorithm : public BareAlgorithm
{
public:
  struct Config
  {
    // create a config object with the propagator
    Config(propagator_t prop) : propagator(std::move(prop)) {}

    /// the propagors to be tested
    propagator_t propagator;

    /// how to set it up
    std::shared_ptr<RandomNumbersSvc> randomNumberSvc = nullptr;

    /// proapgation mode
    int mode = 0;
    /// debug output
    bool debugOutput = false;
    /// number of particles
    size_t ntests = 100;
    /// d0 gaussian sigma
    double d0Sigma = 15. * Acts::units::_um;
    /// z0 gaussian sigma
    double z0Sigma = 55. * Acts::units::_mm;
    /// phi range
    std::pair<double, double> phiRange = {-M_PI, M_PI};
    /// eta range
    std::pair<double, double> etaRange = {-4., 4.};
    /// pt range
    std::pair<double, double> ptRange
        = {100. * Acts::units::_MeV, 100. * Acts::units::_GeV};
    /// looper protection
    double ptLoopers = 300. * Acts::units::_MeV;

    /// Max step size steering
    double maxStepSize = 1. * Acts::units::_mm;

    /// the step collection to be stored
    std::string propagationStepCollection = "PropagationSteps";

    /// covariance transport
    bool covarianceTransport = true;

    /// the covariance values
    Acts::ActsVectorD<5> covariances = Acts::ActsVectorD<5>::Zero();

    /// the correlation terms
    Acts::ActsSymMatrixD<5> correlations = Acts::ActsSymMatrixD<5>::Identity();
  };

  /// Constructor
  /// @param [in] cnf is the configuration struct
  /// @param [in] loglevel is the loggin level
  PropagationAlgorithm(const Config& cnf, Acts::Logging::Level loglevel);

  /// Framework execute method
  /// @param [in] the algorithm context for event consistency
  /// @return is a process code indicating succes or not
  FW::ProcessCode
  execute(AlgorithmContext ctx) const final override;

private:
  Config m_cfg;  ///< the config class

  /// Private helper method to create a corrleated covariance matrix
  /// @param[in] rnd is the random engine
  /// @param[in] gauss is a gaussian distribution to draw from
  std::unique_ptr<Acts::ActsSymMatrixD<5>>
  generateCovariance(FW::RandomEngine& rnd, FW::GaussDist& gauss) const;

  /// Templated execute test method for
  /// charged and netural particles
  ///
  // @tparam parameters_t type of the parameters objects (charged/neutra;)
  ///
  /// @param [in] startParameters the start parameters
  ///
  /// @return collection of Propagation steps for further analysis
  template <typename parameters_t>
  std::vector<Acts::detail::Step>
  executeTest(const parameters_t& startParameters,
              double pathLength = std::numeric_limits<double>::max()) const;
};

#include "PropagationAlgorithm.ipp"

}  // namespace FW
