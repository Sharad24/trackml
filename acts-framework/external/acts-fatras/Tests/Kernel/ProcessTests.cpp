// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE AbortList Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/Utilities/Definitions.hpp"
#include "Fatras/Kernel/PhysicsList.hpp"
#include "Fatras/Kernel/Process.hpp"
#include "Fatras/Kernel/SelectorList.hpp"
#include "Particle.hpp"
#include <algorithm>
#include <random>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {

namespace Test {

/// the generator
typedef std::mt19937 Generator;

/// The detector
struct Detector {};

/// The selector
struct Selector {

  /// call operator
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &) const {
    return true;
  }
};

/// The scattering formula
struct EnergyDecreaser {

  // constant 10 percent of enery loss
  double cvalue = 0.90;

  template <typename generator_t, typename detector_t, typename particle_t>
  std::vector<particle_t> operator()(generator_t &, const detector_t &,
                                     particle_t &in) const {

    in.energyLoss((1. - cvalue) * in.E());
    return {};
  }
};

/// Test the scattering implementation
BOOST_DATA_TEST_CASE(
    Process_test_,
    bdata::random(
        (bdata::seed = 20,
         bdata::distribution = std::uniform_real_distribution<>(0., 1.))) ^
        bdata::random(
            (bdata::seed = 21,
             bdata::distribution = std::uniform_real_distribution<>(0., 1.))) ^
        bdata::random(
            (bdata::seed = 22,
             bdata::distribution = std::uniform_real_distribution<>(0., 1.))) ^
        bdata::random((
            bdata::seed = 23,
            bdata::distribution = std::uniform_real_distribution<>(1., 100.))) ^
        bdata::xrange(100),
    x, y, z, p, index) {
  // standard generator
  Generator generator;

  // Dummy detctor
  Detector detector;

  // create the particle and set the momentum
  /// position at 0.,0.,0
  Acts::Vector3D position{0., 0., 0.};
  // pT of 1 GeV
  Acts::Vector3D momentum =
      p * Acts::units::_GeV * Acts::Vector3D(x, y, z).normalized();
  // positively charged
  double q = 1.;
  double m = 105.658367 * Acts::units::_MeV; // muon mass

  // create the particle
  Particle particle(position, momentum, m, q, 13, 1);

  // outgoing particles (always none for scattering)
  std::vector<Particle> outgoing;

  // T"{he select all list
  typedef SelectorListAND<Selector> All;
  typedef Process<EnergyDecreaser, All, All, All> EnergyLoss;
  EnergyLoss cEnergyLoss;

  // energy loss is not allowed to throw abort command
  BOOST_CHECK(!cEnergyLoss(generator, detector, particle, outgoing));

  // check the the particle momentum magnitude is NOT identical
  BOOST_CHECK(momentum.norm() != particle.momentum().norm());

  // let's test this as part of a physics list
  PhysicsList<EnergyLoss> energyLossPhysics;

  // scattering is not allowed to throw abort command
  BOOST_CHECK(!energyLossPhysics(generator, detector, particle, outgoing));
}

} // namespace Test
} // namespace Fatras
