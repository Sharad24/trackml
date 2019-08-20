// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE KinemtaicCast Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/Utilities/Units.hpp"
#include "Fatras/Selectors/KinematicCasts.hpp"
#include "Particle.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {

namespace Test {

double m = 134.9766 * Acts::units::_MeV;

// This tests the implementation of kinematic cast operators
BOOST_AUTO_TEST_CASE(Kinematic_cast_tests) {

  // a central pion
  Acts::Vector3D position(0., 0., 0.);
  Acts::Vector3D momentumCentral(1500. * Acts::units::_MeV, 0., 0.);
  Particle pionCentral(position, momentumCentral, m, -1.);

  // a forward pion
  Acts::Vector3D positionFwd(0., 0., 100.);
  Acts::Vector3D momentumFwd(10. * Acts::units::_MeV, 10. * Acts::units::_MeV,
                             1500. * Acts::units::_MeV);
  Particle pionFwd(positionFwd, momentumFwd, m, -1.);

  // the list of possible casts
  casts::eta eta_c;
  casts::absEta absEta_c;
  casts::pT pT_c;
  casts::p p_c;
  casts::E E_c;
  casts::vR vR_c;
  casts::vZ vZ_c;

  // test the central
  BOOST_TEST(eta_c(pionCentral) == 0., tt::tolerance(1e-10));
  BOOST_TEST(absEta_c(pionCentral) == 0., tt::tolerance(1e-10));
  BOOST_TEST(pT_c(pionCentral), 1500. * Acts::units::_MeV);
  BOOST_TEST(p_c(pionCentral), 1500. * Acts::units::_MeV);
  BOOST_CHECK(E_c(pionCentral) > p_c(pionCentral));

  BOOST_CHECK_CLOSE(vR_c(pionCentral), 0., 10e-5);
  BOOST_CHECK_CLOSE(vZ_c(pionCentral), 0., 10e-5);

  // test the forward
  BOOST_CHECK(eta_c(pionFwd) > eta_c(pionCentral));
  BOOST_TEST(vZ_c(pionFwd), 100. * Acts::units::_MeV);
}

} // namespace Test
} // namespace Fatras
