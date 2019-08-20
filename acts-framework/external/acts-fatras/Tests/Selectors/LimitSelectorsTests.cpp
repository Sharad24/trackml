// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE LimitSelector Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Fatras/Selectors/LimitSelectors.hpp"
#include "Particle.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {

namespace Test {

// some material
Acts::Material berilium = Acts::Material(352.8, 407., 9.012, 4., 1.848e-3);

double m = 134.9766 * Acts::units::_MeV;

// This tests the implementation of kinematic cast operators
BOOST_AUTO_TEST_CASE(Kinematic_cast_tests) {

  Acts::MaterialProperties detector(berilium, 1. * Acts::units::_mm);

  // a central pion
  Acts::Vector3D position(0., 0., 0.);

  Acts::Vector3D momentum(1500. * Acts::units::_MeV, 0., 0.);
  Particle pion(position, momentum, m, -1.);

  // the limit of the particle
  pion.setLimits(0.15, 0.45);
  // the path of the particle
  pion.update(position, momentum, 0.10, 0.34);

  X0Limit x0LimitSelector;
  L0Limit l0LimitSelector;
  // the limit is not yet reached
  BOOST_CHECK(!x0LimitSelector(detector, pion));
  BOOST_CHECK(!l0LimitSelector(detector, pion));

  detector = Acts::MaterialProperties(berilium, 150. * Acts::units::_mm);

  // the limit is now reached
  BOOST_CHECK(x0LimitSelector(detector, pion));
  BOOST_CHECK(l0LimitSelector(detector, pion));
}

} // end of Test namespace
} // end of Fatras namespace
