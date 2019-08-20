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
#include "Fatras/Selectors/SelectorHelpers.hpp"
#include "Particle.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {

namespace Test {

struct Detector {};

double m_pion = 134.9766 * Acts::units::_MeV; // pi0 rest mass

// This tests the implementation of kinematic cast operators
BOOST_AUTO_TEST_CASE(SelectorHelper_tests) {

  Detector detector;

  Acts::Vector3D position(0., 0., 0.);
  Acts::Vector3D momentumCast(1500. * Acts::units::_MeV, 0., 0.);

  // e central electron
  Particle pionCast(position, momentumCast, m_pion, -1.);
  Acts::Vector3D positionForward(0., 0., 100. * Acts::units::_mm);
  Acts::Vector3D momentumForward(10. * Acts::units::_MeV,
                                 10. * Acts::units::_MeV,
                                 1500. * Acts::units::_MeV);

  Particle pionForward(positionForward, momentumForward, m_pion, -1.);

  Acts::Vector3D positionBackward(0., 0., 0.);
  Acts::Vector3D momentumBackward(10. * Acts::units::_MeV,
                                  10. * Acts::units::_MeV,
                                  -1500. * Acts::units::_MeV);

  Particle pionBackward(positionBackward, momentumBackward, m_pion, -1.);

  // the list of possible casts
  casts::eta etaCast;
  casts::absEta etaAbsCast;

  // A minimum of 0.5 Eta is required
  Min<casts::eta> minEta05;
  minEta05.valMin = 0.5;

  Min<casts::absEta> minEtaAbs05;
  minEtaAbs05.valMin = 0.5;

  // the central will fail both
  BOOST_CHECK(!minEta05(detector, pionCast));
  BOOST_CHECK(!minEtaAbs05(detector, pionCast));

  // the forward will satisfy both
  BOOST_CHECK(minEta05(detector, pionForward));
  BOOST_CHECK(minEtaAbs05(detector, pionForward));

  // A maximum of 0.5 Eta is required
  Max<casts::eta> maxEta05;
  maxEta05.valMax = 0.5;

  // the central will satisfy both
  BOOST_CHECK(maxEta05(detector, pionCast));
  BOOST_CHECK(maxEta05(detector, pionCast));

  // the forward will fail both
  BOOST_CHECK(!maxEta05(detector, pionForward));
  BOOST_CHECK(!maxEta05(detector, pionForward));

  // a range test
  Range<casts::eta> rangeEtaM0;
  rangeEtaM0.valMin = -6.;
  rangeEtaM0.valMax = -0.5;

  Range<casts::absEta> rangeEtaM1;
  rangeEtaM1.valMin = -6.;
  rangeEtaM1.valMax = -0.5;

  // the central will fail both
  BOOST_CHECK(!rangeEtaM0(detector, pionCast));
  BOOST_CHECK(!rangeEtaM1(detector, pionCast));

  // the forward will fail both
  BOOST_CHECK(!rangeEtaM0(detector, pionForward));
  BOOST_CHECK(!rangeEtaM1(detector, pionForward));

  // the backeard will satsify the eta cast, not the abs(eta)
  BOOST_CHECK(rangeEtaM0(detector, pionBackward));
  BOOST_CHECK(!rangeEtaM1(detector, pionBackward));
}

} // namespace Test
} // namespace Fatras
