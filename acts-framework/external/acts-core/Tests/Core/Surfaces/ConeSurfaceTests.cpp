// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE ConeSurface Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>

#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"

namespace tt = boost::test_tools;
using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(ConeSurfaces)
  /// Unit test for creating compliant/non-compliant ConeSurface object
  BOOST_AUTO_TEST_CASE(ConeSurfaceConstruction)
  {
    // ConeSurface default constructor is deleted
    //
    /// Constructor with transform pointer, null or valid, alpha and symmetry
    /// indicator
    double alpha{M_PI / 8.}, halfPhiSector{M_PI / 16.}, zMin{1.0}, zMax{10.};
    bool   symmetric(false);
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    auto          pNullTransform = std::make_shared<const Transform3D>();
    BOOST_CHECK_EQUAL(
        Surface::makeShared<ConeSurface>(pNullTransform, alpha, symmetric)
            ->type(),
        Surface::Cone);
    BOOST_CHECK_EQUAL(
        Surface::makeShared<ConeSurface>(pTransform, alpha, symmetric)->type(),
        Surface::Cone);
    //
    /// Constructor with transform pointer, alpha,z min and max, halfPhiSector
    BOOST_CHECK_EQUAL(Surface::makeShared<ConeSurface>(
                          pTransform, alpha, zMin, zMax, halfPhiSector)
                          ->type(),
                      Surface::Cone);
    //

    /// Constructor with transform and ConeBounds pointer
    // ConeBounds (double alpha, double zmin, double zmax, double halfphi=M_PI,
    // double avphi=0.)
    auto pConeBounds = std::make_shared<const ConeBounds>(
        alpha, zMin, zMax, halfPhiSector, 0.);
    BOOST_CHECK_EQUAL(
        Surface::makeShared<ConeSurface>(pTransform, pConeBounds)->type(),
        Surface::Cone);
    //
    //
    /// Copy constructor
    auto coneSurfaceObject
        = Surface::makeShared<ConeSurface>(pTransform, alpha, symmetric);
    auto copiedConeSurface
        = Surface::makeShared<ConeSurface>(*coneSurfaceObject);
    BOOST_CHECK_EQUAL(copiedConeSurface->type(), Surface::Cone);
    BOOST_CHECK_EQUAL(*copiedConeSurface, *coneSurfaceObject);
    //
    /// Copied and transformed
    auto copiedTransformedConeSurface
        = Surface::makeShared<ConeSurface>(*coneSurfaceObject, *pTransform);
    BOOST_CHECK_EQUAL(copiedTransformedConeSurface->type(), Surface::Cone);

    /// Construct with nullptr bounds
    BOOST_CHECK_THROW(auto nullBounds
                      = Surface::makeShared<ConeSurface>(nullptr, nullptr),
                      AssertionFailureException);
  }
  //
  /// Unit test for testing ConeSurface properties
  BOOST_AUTO_TEST_CASE(ConeSurfaceProperties)
  {
    /// Test clone method
    double alpha{M_PI / 8.} /*,halfPhiSector{M_PI/16.}, zMin{1.0}, zMax{10.}*/;
    bool   symmetric(false);
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    // auto pNullTransform = std::make_shared<const Transform3D>();
    auto coneSurfaceObject
        = Surface::makeShared<ConeSurface>(pTransform, alpha, symmetric);
    //
    auto pClonedConeSurface = coneSurfaceObject->clone();
    BOOST_CHECK_EQUAL(pClonedConeSurface->type(), Surface::Cone);
    //
    /// Test type (redundant)
    BOOST_CHECK_EQUAL(coneSurfaceObject->type(), Surface::Cone);
    //
    /// Test binningPosition
    Vector3D binningPosition{0., 1., 2.};
    CHECK_CLOSE_ABS(coneSurfaceObject->binningPosition(BinningValue::binPhi),
                    binningPosition,
                    1e-6);
    //
    /// Test referenceFrame
    Vector3D         globalPosition{2.0, 2.0, 2.0};
    Vector3D         momentum{1.e6, 1.e6, 1.e6};
    double           rootHalf = std::sqrt(0.5);
    RotationMatrix3D expectedFrame;
    expectedFrame << -rootHalf, 0., rootHalf, rootHalf, 0., rootHalf, 0., 1.,
        0.;
    CHECK_CLOSE_OR_SMALL(
        coneSurfaceObject->referenceFrame(globalPosition, momentum),
        expectedFrame,
        1e-6,
        1e-9);
    //
    /// Test normal, given 3D position
    Vector3D origin{0., 0., 0.};
    Vector3D normal3D = {0., -1., 0.};
    CHECK_CLOSE_ABS(coneSurfaceObject->normal(origin), normal3D, 1e-6);
    //
    /// Test normal given 2D rphi position
    Vector2D positionPiBy2(1.0, M_PI / 2.);
    Vector3D normalAtPiBy2{0.0312768, 0.92335, -0.382683};

    CHECK_CLOSE_OR_SMALL(
        coneSurfaceObject->normal(positionPiBy2), normalAtPiBy2, 1e-2, 1e-9);
    //
    /// Test rotational symmetry axis
    Vector3D symmetryAxis{0., 0., 1.};
    CHECK_CLOSE_ABS(coneSurfaceObject->rotSymmetryAxis(), symmetryAxis, 1e-6);
    //
    /// Test bounds
    BOOST_CHECK_EQUAL(coneSurfaceObject->bounds().type(), SurfaceBounds::Cone);
    //
    /// Test localToGlobal
    Vector2D localPosition{1.0, M_PI / 2.0};
    coneSurfaceObject->localToGlobal(localPosition, momentum, globalPosition);
    // std::cout<<globalPosition<<std::endl;
    Vector3D expectedPosition{0.0220268, 1.65027, 3.5708};

    CHECK_CLOSE_REL(globalPosition, expectedPosition, 1e-2);
    //
    /// Testing globalToLocal
    coneSurfaceObject->globalToLocal(globalPosition, momentum, localPosition);
    // std::cout<<localPosition<<std::endl;
    Vector2D expectedLocalPosition{1.0, M_PI / 2.0};

    CHECK_CLOSE_REL(localPosition, expectedLocalPosition, 1e-6);
    //
    /// Test isOnSurface
    Vector3D offSurface{100, 1, 2};
    BOOST_CHECK(coneSurfaceObject->isOnSurface(globalPosition, momentum, true));
    BOOST_CHECK(!coneSurfaceObject->isOnSurface(offSurface, momentum, true));
    //
    /// intersectionEstimate
    Vector3D direction{-1., 0, 0};
    auto     intersect = coneSurfaceObject->intersectionEstimate(
        offSurface, direction, forward, false);
    Intersection expectedIntersect{Vector3D{0, 1, 2}, 100., true, 0};
    BOOST_CHECK(intersect.valid);
    CHECK_CLOSE_ABS(intersect.position, expectedIntersect.position, 1e-6);
    CHECK_CLOSE_ABS(intersect.pathLength, expectedIntersect.pathLength, 1e-6);
    CHECK_CLOSE_ABS(intersect.distance, expectedIntersect.distance, 1e-6);
    //
    /// Test pathCorrection
    CHECK_CLOSE_REL(coneSurfaceObject->pathCorrection(offSurface, momentum),
                    0.40218866453252877,
                    0.01);
    //
    /// Test name
    BOOST_CHECK_EQUAL(coneSurfaceObject->name(),
                      std::string("Acts::ConeSurface"));
    //
    /// Test dump
    // TODO 2017-04-12 msmk: check how to correctly check output
    //    boost::test_tools::output_test_stream dumpOuput;
    //    coneSurfaceObject.dump(dumpOuput);
    //    BOOST_CHECK(dumpOuput.is_equal(
    //      "Acts::ConeSurface\n"
    //      "    Center position  (x, y, z) = (0.0000, 1.0000, 2.0000)\n"
    //      "    Rotation:             colX = (1.000000, 0.000000, 0.000000)\n"
    //      "                          colY = (0.000000, 1.000000, 0.000000)\n"
    //      "                          colZ = (0.000000, 0.000000, 1.000000)\n"
    //      "    Bounds  : Acts::ConeBounds: (tanAlpha, minZ, maxZ, averagePhi,
    //      halfPhiSector) = (0.4142136, 0.0000000, inf, 0.0000000,
    //      3.1415927)"));
  }

  BOOST_AUTO_TEST_CASE(EqualityOperators)
  {
    double alpha{M_PI / 8.} /*, halfPhiSector{M_PI/16.}, zMin{1.0}, zMax{10.}*/;
    bool   symmetric(false);
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    auto          coneSurfaceObject
        = Surface::makeShared<ConeSurface>(pTransform, alpha, symmetric);
    //
    auto coneSurfaceObject2
        = Surface::makeShared<ConeSurface>(pTransform, alpha, symmetric);
    //
    /// Test equality operator
    BOOST_CHECK_EQUAL(*coneSurfaceObject, *coneSurfaceObject2);
    //
    BOOST_TEST_CHECKPOINT(
        "Create and then assign a ConeSurface object to the existing one");
    /// Test assignment
    auto assignedConeSurface
        = Surface::makeShared<ConeSurface>(nullptr, 0.1, true);
    *assignedConeSurface = *coneSurfaceObject;
    /// Test equality of assigned to original
    BOOST_CHECK_EQUAL(*assignedConeSurface, *coneSurfaceObject);
  }

  BOOST_AUTO_TEST_CASE(ConeSurface_toVariantData)
  {
    double        alpha = M_PI / 2., zMin = 1, zMax = 5, halfPhi = M_PI;
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    auto          cone       = Surface::makeShared<ConeSurface>(
        pTransform, alpha, zMin, zMax, halfPhi);

    variant_data var_cone = cone->toVariantData();
    std::cout << var_cone << std::endl;

    const variant_map& pl
        = boost::get<variant_map>(var_cone).get<variant_map>("payload");
    const variant_map& bounds_pl
        = pl.get<variant_map>("bounds").get<variant_map>("payload");
    BOOST_CHECK_EQUAL(bounds_pl.get<double>("alpha"), alpha);
    BOOST_CHECK_EQUAL(bounds_pl.get<double>("zMin"), zMin);
    BOOST_CHECK_EQUAL(bounds_pl.get<double>("zMax"), zMax);
    BOOST_CHECK_EQUAL(bounds_pl.get<double>("halfPhi"), halfPhi);

    auto cone2      = Surface::makeShared<ConeSurface>(var_cone);
    auto conebounds = dynamic_cast<const ConeBounds*>(&cone2->bounds());
    BOOST_CHECK_EQUAL(conebounds->alpha(), alpha);
    BOOST_CHECK_EQUAL(conebounds->halfPhiSector(), halfPhi);
    BOOST_CHECK_EQUAL(conebounds->minZ(), zMin);
    BOOST_CHECK_EQUAL(conebounds->maxZ(), zMax);
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
