// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Ellipse Bounds Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>

#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit test for creating compliant/non-compliant EllipseBounds object
  BOOST_AUTO_TEST_CASE(EllipseBoundsConstruction)
  {
    double minRad1(10.), minRad2(15.), maxRad1(15.), maxRad2(20.),
        averagePhi(0.), phiSector(M_PI / 2.);
    // test default construction
    // EllipseBounds defaultConstructedEllipseBounds;  //deleted
    //
    /// Test construction with dimensions
    BOOST_CHECK_EQUAL(
        EllipseBounds(minRad1, minRad2, maxRad1, maxRad2, averagePhi, phiSector)
            .type(),
        SurfaceBounds::Ellipse);
    //
    /// Copy constructor
    EllipseBounds original(
        minRad1, minRad2, maxRad1, maxRad2, averagePhi, phiSector);
    EllipseBounds copied(original);
    BOOST_CHECK_EQUAL(copied.type(), SurfaceBounds::Ellipse);
  }

  /// Unit tests for EllipseBounds properties
  BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(CylinderBoundsProperties, 1)
  BOOST_AUTO_TEST_CASE(EllipseBoundsProperties)
  {
    double minRad1(10.), minRad2(15.), maxRad1(15.), maxRad2(20.),
        averagePhi(0.), phiSector(M_PI / 2.);
    /// Test clone
    EllipseBounds ellipseBoundsObject(
        minRad1, minRad2, maxRad1, maxRad2, averagePhi, phiSector);
    auto pClonedEllipseBoundsObject = ellipseBoundsObject.clone();
    BOOST_CHECK_NE(pClonedEllipseBoundsObject, nullptr);
    delete pClonedEllipseBoundsObject;
    //
    /// Test type() (redundant; already used in constructor confirmation)
    BOOST_CHECK_EQUAL(ellipseBoundsObject.type(), SurfaceBounds::Ellipse);
    //
    // clone already tested
    //
    /// Test distanceToBoundary
    Vector2D origin(0., 0.);
    Vector2D outsideBy10(0., 30.);
    Vector2D inRectangle(17., 11.);
    CHECK_CLOSE_REL(ellipseBoundsObject.distanceToBoundary(origin),
                    10.,
                    1e-6);  // makes sense
    CHECK_CLOSE_REL(ellipseBoundsObject.distanceToBoundary(outsideBy10),
                    10.,
                    1e-6);  // fails, not clear why
    //
    /// Test rMinX
    BOOST_CHECK_EQUAL(ellipseBoundsObject.rMinX(), minRad1);
    //
    /// Test rMinY
    BOOST_CHECK_EQUAL(ellipseBoundsObject.rMinY(), minRad2);
    //
    /// Test rMaxX
    BOOST_CHECK_EQUAL(ellipseBoundsObject.rMaxX(), maxRad1);
    //
    /// Test rMaxY
    BOOST_CHECK_EQUAL(ellipseBoundsObject.rMaxY(), maxRad2);
    //
    /// Test averagePhi
    BOOST_CHECK_EQUAL(ellipseBoundsObject.averagePhi(), averagePhi);
    //
    /// Test vertices
    std::vector<Vector2D> expectedVertices{
        {15, 0}, {0, 20}, {-15, 0}, {0, -20}};
    const auto& actualVertices = ellipseBoundsObject.vertices();
    BOOST_CHECK_EQUAL_COLLECTIONS(actualVertices.cbegin(),
                                  actualVertices.cend(),
                                  expectedVertices.cbegin(),
                                  expectedVertices.cend());
    //
    /// Test boundingBox
    BOOST_CHECK_EQUAL(ellipseBoundsObject.boundingBox(),
                      RectangleBounds(15., 20.));
    //
    /// Test halfPhiSector
    BOOST_CHECK_EQUAL(ellipseBoundsObject.halfPhiSector(), M_PI / 2.);
    //
    /// Test dump
    boost::test_tools::output_test_stream dumpOuput;
    ellipseBoundsObject.dump(dumpOuput);
    BOOST_CHECK(dumpOuput.is_equal(
        "Acts::EllipseBounds:  (innerRadiusX, innerRadiusY, outerRadiusX, "
        "outerRadiusY, hPhiSector) = (10.0000000, 15.0000000, 15.0000000, "
        "20.0000000, 0.0000000, 1.5707963)"));
    //
    /// Test inside
    BOOST_CHECK(!ellipseBoundsObject.inside(inRectangle, BoundaryCheck(true)));
    // dont understand why this is so:
    BOOST_CHECK(!ellipseBoundsObject.inside(outsideBy10, BoundaryCheck(true)));
  }
  /// Unit test for testing EllipseBounds assignment
  BOOST_AUTO_TEST_CASE(EllipseBoundsAssignment)
  {
    double minRad1(10.), minRad2(15.), maxRad1(15.), maxRad2(20.),
        averagePhi(0.), phiSector(M_PI / 2.);
    EllipseBounds ellipseBoundsObject(
        minRad1, minRad2, maxRad1, maxRad2, averagePhi, phiSector);
    EllipseBounds similarlyConstructeEllipseBoundsObject(
        minRad1, minRad2, maxRad1, maxRad2, averagePhi, phiSector);
    /// Test operator ==
    BOOST_CHECK_EQUAL(ellipseBoundsObject,
                      similarlyConstructeEllipseBoundsObject);
    //
    /// Test assignment
    EllipseBounds assignedEllipseBoundsObject(11., 12., 17., 18., 1.);
    // object, in some sense
    assignedEllipseBoundsObject = ellipseBoundsObject;
    BOOST_CHECK_EQUAL(assignedEllipseBoundsObject, ellipseBoundsObject);
  }

  BOOST_AUTO_TEST_CASE(EllipseBounds_toVariantData)
  {
    double minRad1(10.), minRad2(15.), maxRad1(15.), maxRad2(20.),
        averagePhi(0.), phiSector(M_PI / 2.);
    EllipseBounds ell(
        minRad1, minRad2, maxRad1, maxRad2, averagePhi, phiSector);
    variant_data var_data = ell.toVariantData();

    std::cout << var_data << std::endl;

    variant_map var_map = boost::get<variant_map>(var_data);
    BOOST_CHECK_EQUAL(var_map.get<std::string>("type"), "EllipseBounds");
    variant_map pl = var_map.get<variant_map>("payload");
    BOOST_CHECK_EQUAL(pl.get<double>("rMinX"), minRad1);
    BOOST_CHECK_EQUAL(pl.get<double>("rMinY"), minRad2);
    BOOST_CHECK_EQUAL(pl.get<double>("rMaxX"), maxRad1);
    BOOST_CHECK_EQUAL(pl.get<double>("rMaxY"), maxRad2);
    BOOST_CHECK_EQUAL(pl.get<double>("avgPhi"), averagePhi);
    BOOST_CHECK_EQUAL(pl.get<double>("halfPhi"), phiSector);

    EllipseBounds ell2(var_data);
    BOOST_CHECK_EQUAL(ell.rMinX(), ell2.rMinX());
    BOOST_CHECK_EQUAL(ell.rMinY(), ell2.rMinY());
    BOOST_CHECK_EQUAL(ell.rMaxX(), ell2.rMaxX());
    BOOST_CHECK_EQUAL(ell.rMaxY(), ell2.rMaxY());
    BOOST_CHECK_EQUAL(ell.averagePhi(), ell2.averagePhi());
    BOOST_CHECK_EQUAL(ell.halfPhiSector(), ell2.halfPhiSector());
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
