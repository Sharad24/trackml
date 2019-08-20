// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE SurfaceArrayCreator
#include <boost/test/included/unit_test.hpp>
#include <boost/format.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include <fstream>
#include <random>

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tools/SurfaceArrayCreator.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"

using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::perp;

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

#define CHECK_ROTATION_ANGLE(t, a, tolerance)                                  \
  {                                                                            \
    Vector3D v = (*t) * Vector3D(1, 0, 0);                                     \
    CHECK_CLOSE_ABS(phi(v), (a), tolerance);                                   \
  }

  using SrfVec = std::vector<std::shared_ptr<const Surface>>;

  struct SurfaceArrayCreatorFixture
  {
    SurfaceArrayCreator                         m_SAC;
    std::vector<std::shared_ptr<const Surface>> m_surfaces;

    SurfaceArrayCreatorFixture()
      : m_SAC(SurfaceArrayCreator::Config(),
              Acts::getDefaultLogger("SurfaceArrayCreator",
                                     Acts::Logging::VERBOSE))
    {
      BOOST_TEST_MESSAGE("setup fixture");
    }
    ~SurfaceArrayCreatorFixture() { BOOST_TEST_MESSAGE("teardown fixture"); }

    template <typename... Args>
    SurfaceArrayCreator::ProtoAxis
    createEquidistantAxis(Args&&... args)
    {
      return m_SAC.createEquidistantAxis(std::forward<Args>(args)...);
    }

    template <typename... Args>
    SurfaceArrayCreator::ProtoAxis
    createVariableAxis(Args&&... args)
    {
      return m_SAC.createVariableAxis(std::forward<Args>(args)...);
    }

    template <detail::AxisBoundaryType bdtA,
              detail::AxisBoundaryType bdtB,
              typename... Args>
    std::unique_ptr<SurfaceArray::ISurfaceGridLookup>
    makeSurfaceGridLookup2D(Args&&... args)
    {
      return m_SAC.makeSurfaceGridLookup2D<bdtA, bdtB>(
          std::forward<Args>(args)...);
    }

    SrfVec
    fullPhiTestSurfacesEC(size_t n     = 10,
                          double shift = 0,
                          double zbase = 0,
                          double r     = 10,
                          double w     = 2,
                          double h     = 1)
    {

      SrfVec res;
      // TODO: The test is extremely numerically unstable in the face of upward
      //       rounding in this multiplication and division. Find out why.
      double phiStep = 2 * M_PI / n;
      for (size_t i = 0; i < n; ++i) {
        double z   = zbase + ((i % 2 == 0) ? 1 : -1) * 0.2;
        double phi = std::fma(i, phiStep, shift);

        Transform3D trans;
        trans.setIdentity();
        trans.rotate(Eigen::AngleAxisd(phi, Vector3D(0, 0, 1)));
        trans.translate(Vector3D(r, 0, z));

        auto bounds = std::make_shared<const RectangleBounds>(w, h);

        auto transptr = std::make_shared<const Transform3D>(trans);
        std::shared_ptr<Surface> srf
            = Surface::makeShared<PlaneSurface>(transptr, bounds);

        res.push_back(srf);
        m_surfaces.push_back(
            std::move(srf));  // keep shared, will get destroyed at the end
      }

      return res;
    }

    SrfVec
    fullPhiTestSurfacesBRL(size_t n     = 10,
                           double shift = 0,
                           double zbase = 0,
                           double incl  = M_PI / 9.,
                           double w     = 2,
                           double h     = 1.5)
    {

      SrfVec res;
      // TODO: The test is extremely numerically unstable in the face of upward
      //       rounding in this multiplication and division. Find out why.
      double phiStep = 2 * M_PI / n;
      for (size_t i = 0; i < n; ++i) {
        double z   = zbase;
        double phi = std::fma(i, phiStep, shift);

        Transform3D trans;
        trans.setIdentity();
        trans.rotate(Eigen::AngleAxisd(phi, Vector3D(0, 0, 1)));
        trans.translate(Vector3D(10, 0, z));
        trans.rotate(Eigen::AngleAxisd(incl, Vector3D(0, 0, 1)));
        trans.rotate(Eigen::AngleAxisd(M_PI / 2., Vector3D(0, 1, 0)));

        auto bounds = std::make_shared<const RectangleBounds>(w, h);

        auto transptr = std::make_shared<const Transform3D>(trans);
        std::shared_ptr<Surface> srf
            = Surface::makeShared<PlaneSurface>(transptr, bounds);

        res.push_back(srf);
        m_surfaces.push_back(
            std::move(srf));  // keep shared, will get destroyed at the end
      }

      return res;
    }

    SrfVec
    straightLineSurfaces(size_t             n        = 10.,
                         double             step     = 3,
                         const Vector3D&    origin   = {0, 0, 1.5},
                         const Transform3D& pretrans = Transform3D::Identity(),
                         const Vector3D&    dir      = {0, 0, 1})
    {
      SrfVec res;
      for (size_t i = 0; i < n; ++i) {
        Transform3D trans;
        trans.setIdentity();
        trans.translate(origin + dir * step * i);
        // trans.rotate(AngleAxis3D(M_PI/9., Vector3D(0, 0, 1)));
        trans.rotate(AngleAxis3D(M_PI / 2., Vector3D(1, 0, 0)));
        trans = trans * pretrans;

        auto bounds = std::make_shared<const RectangleBounds>(2, 1.5);

        auto transptr = std::make_shared<const Transform3D>(trans);
        std::shared_ptr<Surface> srf
            = Surface::makeShared<PlaneSurface>(transptr, bounds);

        res.push_back(srf);
        m_surfaces.push_back(
            std::move(srf));  // keep shared, will get destroyed at the end
      }

      return res;
    }

    SrfVec
    makeBarrel(int nPhi, int nZ, double w, double h)
    {
      double z0 = -(nZ - 1) * w;
      SrfVec res;

      for (int i = 0; i < nZ; i++) {
        double z = i * w * 2 + z0;
        // std::cout << "z=" << z << std::endl;
        SrfVec ring = fullPhiTestSurfacesBRL(nPhi, 0, z, M_PI / 9., w, h);
        res.insert(res.end(), ring.begin(), ring.end());
      }

      return res;
    }

    std::pair<SrfVec, std::vector<std::pair<const Surface*, const Surface*>>>
    makeBarrelStagger(int    nPhi,
                      int    nZ,
                      double shift = 0,
                      double incl  = M_PI / 9.,
                      double w     = 2,
                      double h     = 1.5)
    {
      double z0 = -(nZ - 1) * w;
      SrfVec res;
      std::vector<std::pair<const Surface*, const Surface*>> pairs;
      // TODO: The test is extremely numerically unstable in the face of upward
      //       rounding in this multiplication and division. Find out why.
      double phiStep = 2 * M_PI / nPhi;
      for (int i = 0; i < nZ; i++) {
        double z = i * w * 2 + z0;
        for (int j = 0; j < nPhi; ++j) {
          double      phi = std::fma(j, phiStep, shift);
          Transform3D trans;
          trans.setIdentity();
          trans.rotate(Eigen::AngleAxisd(phi, Vector3D(0, 0, 1)));
          trans.translate(Vector3D(10, 0, z));
          trans.rotate(Eigen::AngleAxisd(incl, Vector3D(0, 0, 1)));
          trans.rotate(Eigen::AngleAxisd(M_PI / 2., Vector3D(0, 1, 0)));

          auto bounds = std::make_shared<const RectangleBounds>(w, h);

          auto transAptr = std::make_shared<const Transform3D>(trans);

          std::shared_ptr<Surface> srfA
              = Surface::makeShared<PlaneSurface>(transAptr, bounds);

          Vector3D    nrm    = srfA->normal();
          Transform3D transB = trans;
          transB.pretranslate(nrm * 0.1);
          auto transBptr = std::make_shared<const Transform3D>(transB);
          std::shared_ptr<Surface> srfB
              = Surface::makeShared<PlaneSurface>(transBptr, bounds);

          pairs.push_back(std::make_pair(srfA.get(), srfB.get()));

          res.push_back(srfA);
          res.push_back(srfB);
          m_surfaces.push_back(std::move(srfA));
          m_surfaces.push_back(std::move(srfB));
        }
      }

      return std::make_pair(res, pairs);
    }
  };

  void
  draw_surfaces(SrfVec surfaces, const std::string& fname)
  {

    std::ofstream os;
    os.open(fname);

    os << std::fixed << std::setprecision(4);

    size_t nVtx = 0;
    for (const auto& srfx : surfaces) {
      std::shared_ptr<const PlaneSurface> srf
          = std::dynamic_pointer_cast<const PlaneSurface>(srfx);
      const PlanarBounds* bounds
          = dynamic_cast<const PlanarBounds*>(&srf->bounds());

      for (const auto& vtxloc : bounds->vertices()) {
        Vector3D vtx = srf->transform() * Vector3D(vtxloc.x(), vtxloc.y(), 0);
        os << "v " << vtx.x() << " " << vtx.y() << " " << vtx.z() << "\n";
      }

      // connect them
      os << "f";
      for (size_t i = 1; i <= bounds->vertices().size(); ++i) {
        os << " " << nVtx + i;
      }
      os << "\n";

      nVtx += bounds->vertices().size();
    }

    os.close();
  }

  BOOST_AUTO_TEST_SUITE(Tools)

  BOOST_FIXTURE_TEST_CASE(SurfaceArrayCreator_createEquidistantAxis_Phi,
                          SurfaceArrayCreatorFixture)
  {
    // fail on empty srf vector
    std::vector<const Surface*> emptyRaw;
    ProtoLayer                  pl(emptyRaw);
    auto                        tr = Transform3D::Identity();
    BOOST_CHECK_THROW(
        createEquidistantAxis(emptyRaw, BinningValue::binPhi, pl, tr),
        std::logic_error);

    std::vector<float> bdExp = {
        -3.14159, -2.93215, -2.72271, -2.51327, -2.30383,  -2.0944,   -1.88496,
        -1.67552, -1.46608, -1.25664, -1.0472,  -0.837758, -0.628319, -0.418879,
        -0.20944, 0,        0.20944,  0.418879, 0.628319,  0.837758,  1.0472,
        1.25664,  1.46608,  1.67552,  1.88496,  2.09439,   2.30383,   2.51327,
        2.72271,  2.93215,  3.14159};

    double step = 2 * M_PI / 30.;

    // endcap style modules

    for (int i = -1; i <= 2; i += 2) {
      double z = 10 * i;
      // case 1: one module sits at pi / -pi
      double angleShift = step / 2.;
      auto   surfaces   = fullPhiTestSurfacesEC(30, angleShift, z);
      std::vector<const Surface*> surfacesRaw = unpack_shared_vector(surfaces);
      pl                                      = ProtoLayer(surfacesRaw);
      tr                                      = Transform3D::Identity();
      auto axis
          = createEquidistantAxis(surfacesRaw, BinningValue::binPhi, pl, tr);

      BOOST_CHECK_EQUAL(axis.nBins, 30);
      CHECK_CLOSE_REL(axis.max, M_PI, 1e-6);
      CHECK_CLOSE_REL(axis.min, -M_PI, 1e-6);
      BOOST_CHECK_EQUAL(axis.bType, equidistant);
      CHECK_SMALL(phi(tr * Vector3D::UnitX()), 1e-6);

      // case 2: two modules sit symmetrically around pi / -pi
      angleShift  = 0.;
      surfaces    = fullPhiTestSurfacesEC(30, angleShift, z);
      surfacesRaw = unpack_shared_vector(surfaces);
      pl          = ProtoLayer(surfacesRaw);
      tr          = Transform3D::Identity();
      axis = createEquidistantAxis(surfacesRaw, BinningValue::binPhi, pl, tr);
      draw_surfaces(surfaces,
                    "SurfaceArrayCreator_createEquidistantAxis_EC_2.obj");
      BOOST_CHECK_EQUAL(axis.nBins, 30);
      CHECK_CLOSE_REL(axis.max, M_PI, 1e-6);
      CHECK_CLOSE_REL(axis.min, -M_PI, 1e-6);
      BOOST_CHECK_EQUAL(axis.bType, equidistant);
      // CHECK_CLOSE_REL(bdExp, axis.binEdges, 0.001);
      CHECK_CLOSE_REL(phi(tr * Vector3D::UnitX()), -0.5 * step, 1e-3);
      // case 3: two modules sit asymmetrically around pi / -pi shifted up
      angleShift  = step / -4.;
      surfaces    = fullPhiTestSurfacesEC(30, angleShift, z);
      surfacesRaw = unpack_shared_vector(surfaces);
      pl          = ProtoLayer(surfacesRaw);
      tr          = Transform3D::Identity();
      axis = createEquidistantAxis(surfacesRaw, BinningValue::binPhi, pl, tr);
      draw_surfaces(surfaces,
                    "SurfaceArrayCreator_createEquidistantAxis_EC_3.obj");
      BOOST_CHECK_EQUAL(axis.nBins, 30);
      CHECK_CLOSE_REL(axis.max, M_PI, 1e-6);
      CHECK_CLOSE_REL(axis.min, -M_PI, 1e-6);
      BOOST_CHECK_EQUAL(axis.bType, equidistant);
      CHECK_CLOSE_REL(phi(tr * Vector3D::UnitX()), step / -4., 1e-3);

      // case 4: two modules sit asymmetrically around pi / -pi shifted down
      angleShift  = step / 4.;
      surfaces    = fullPhiTestSurfacesEC(30, angleShift, z);
      surfacesRaw = unpack_shared_vector(surfaces);
      pl          = ProtoLayer(surfaces);
      surfacesRaw = unpack_shared_vector(surfaces);
      tr          = Transform3D::Identity();
      axis = createEquidistantAxis(surfacesRaw, BinningValue::binPhi, pl, tr);
      surfacesRaw = unpack_shared_vector(surfaces);
      draw_surfaces(surfaces,
                    "SurfaceArrayCreator_createEquidistantAxis_EC_4.obj");
      BOOST_CHECK_EQUAL(axis.nBins, 30);
      CHECK_CLOSE_REL(axis.max, M_PI, 1e-6);
      CHECK_CLOSE_REL(axis.min, -M_PI, 1e-6);
      BOOST_CHECK_EQUAL(axis.bType, equidistant);
      CHECK_CLOSE_REL(phi(tr * Vector3D::UnitX()), step / 4., 1e-3);
    }

    for (int i = -1; i <= 2; i += 2) {
      double z = 10 * i;
      // case 1: one module sits at pi / -pi
      double angleShift  = step / 2.;
      auto   surfaces    = fullPhiTestSurfacesBRL(30, angleShift, z);
      auto   surfacesRaw = unpack_shared_vector(surfaces);
      pl                 = ProtoLayer(surfacesRaw);
      tr                 = Transform3D::Identity();
      auto axis
          = createEquidistantAxis(surfacesRaw, BinningValue::binPhi, pl, tr);
      draw_surfaces(surfaces,
                    "SurfaceArrayCreator_createEquidistantAxis_BRL_1.obj");
      BOOST_CHECK_EQUAL(axis.nBins, 30);
      CHECK_CLOSE_REL(axis.max, M_PI, 1e-6);
      CHECK_CLOSE_REL(axis.min, -M_PI, 1e-6);
      BOOST_CHECK_EQUAL(axis.bType, equidistant);
      CHECK_SMALL(phi(tr * Vector3D::UnitX()), 1e-6);

      // case 2: two modules sit symmetrically around pi / -pi
      angleShift  = 0.;
      surfaces    = fullPhiTestSurfacesBRL(30, angleShift, z);
      surfacesRaw = unpack_shared_vector(surfaces);
      pl          = ProtoLayer(surfacesRaw);
      tr          = Transform3D::Identity();
      axis = createEquidistantAxis(surfacesRaw, BinningValue::binPhi, pl, tr);
      draw_surfaces(surfaces,
                    "SurfaceArrayCreator_createEquidistantAxis_BRL_2.obj");
      BOOST_CHECK_EQUAL(axis.nBins, 30);
      CHECK_CLOSE_REL(axis.max, M_PI, 1e-6);
      CHECK_CLOSE_REL(axis.min, -M_PI, 1e-6);
      BOOST_CHECK_EQUAL(axis.bType, equidistant);
      // CHECK_CLOSE_REL(bdExp, axis.binEdges, 0.001);
      CHECK_CLOSE_REL(phi(tr * Vector3D::UnitX()), -0.5 * step, 1e-3);

      // case 3: two modules sit asymmetrically around pi / -pi shifted up
      angleShift  = step / -4.;
      surfaces    = fullPhiTestSurfacesBRL(30, angleShift, z);
      surfacesRaw = unpack_shared_vector(surfaces);
      pl          = ProtoLayer(surfacesRaw);
      tr          = Transform3D::Identity();
      axis = createEquidistantAxis(surfacesRaw, BinningValue::binPhi, pl, tr);
      draw_surfaces(surfaces,
                    "SurfaceArrayCreator_createEquidistantAxis_BRL_3.obj");
      BOOST_CHECK_EQUAL(axis.nBins, 30);
      CHECK_CLOSE_REL(axis.max, M_PI, 1e-6);
      CHECK_CLOSE_REL(axis.min, -M_PI, 1e-6);
      BOOST_CHECK_EQUAL(axis.bType, equidistant);
      // CHECK_CLOSE_REL(bdExp, axis.binEdges, 0.001);
      CHECK_CLOSE_REL(phi(tr * Vector3D::UnitX()), step / -4., 1e-3);

      // case 4: two modules sit asymmetrically around pi / -pi shifted down
      angleShift  = step / 4.;
      surfaces    = fullPhiTestSurfacesBRL(30, angleShift, z);
      surfacesRaw = unpack_shared_vector(surfaces);
      pl          = ProtoLayer(surfacesRaw);
      tr          = Transform3D::Identity();
      axis = createEquidistantAxis(surfacesRaw, BinningValue::binPhi, pl, tr);
      draw_surfaces(surfaces,
                    "SurfaceArrayCreator_createEquidistantAxis_BRL_4.obj");
      BOOST_CHECK_EQUAL(axis.nBins, 30);
      CHECK_CLOSE_REL(axis.max, M_PI, 1e-6);
      CHECK_CLOSE_REL(axis.min, -M_PI, 1e-6);
      BOOST_CHECK_EQUAL(axis.bType, equidistant);
      // CHECK_CLOSE_REL(bdExp, axis.binEdges, 0.001);
      CHECK_CLOSE_REL(phi(tr * Vector3D::UnitX()), step / 4., 1e-3);
    }

    SrfVec surfaces;

    // single element in phi
    surfaces         = fullPhiTestSurfacesEC(1);
    auto surfacesRaw = unpack_shared_vector(surfaces);
    draw_surfaces(surfaces,
                  "SurfaceArrayCreator_createEquidistantAxis_EC_Single.obj");

    pl = ProtoLayer(surfacesRaw);
    tr = Transform3D::Identity();
    auto axis
        = createEquidistantAxis(surfacesRaw, BinningValue::binPhi, pl, tr);
    BOOST_CHECK_EQUAL(axis.nBins, 1);

    CHECK_CLOSE_ABS(axis.max, phi(Vector3D(8, 1, 0)), 1e-3);
    CHECK_CLOSE_ABS(axis.min, phi(Vector3D(8, -1, 0)), 1e-3);
    BOOST_CHECK_EQUAL(axis.bType, equidistant);
  }

  BOOST_FIXTURE_TEST_CASE(SurfaceArrayCreator_createEquidistantAxis_Z,
                          SurfaceArrayCreatorFixture)
  {

    // single element in z
    auto       surfaces    = straightLineSurfaces(1);
    auto       surfacesRaw = unpack_shared_vector(surfaces);
    ProtoLayer pl          = ProtoLayer(surfacesRaw);
    auto       trf         = Transform3D::Identity();
    auto axis = createEquidistantAxis(surfacesRaw, BinningValue::binZ, pl, trf);
    draw_surfaces(surfaces,
                  "SurfaceArrayCreator_createEquidistantAxis_Z_1.obj");
    BOOST_CHECK_EQUAL(axis.nBins, 1);
    CHECK_CLOSE_ABS(axis.max, 3, 1e-6);
    CHECK_CLOSE_ABS(axis.min, 0, 1e-6);
    BOOST_CHECK_EQUAL(axis.bType, equidistant);

    // z rows with varying starting point
    for (size_t i = 0; i <= 20; i++) {
      double z0   = -10 + 1. * i;
      surfaces    = straightLineSurfaces(10, 3, Vector3D(0, 0, z0 + 1.5));
      surfacesRaw = unpack_shared_vector(surfaces);
      pl          = ProtoLayer(surfacesRaw);
      trf         = Transform3D::Identity();
      axis = createEquidistantAxis(surfacesRaw, BinningValue::binZ, pl, trf);
      draw_surfaces(
          surfaces,
          (boost::format(
               "SurfaceArrayCreator_createEquidistantAxis_Z_2_%1%.obj")
           % i)
              .str());
      BOOST_CHECK_EQUAL(axis.nBins, 10);
      CHECK_CLOSE_ABS(axis.max, 30 + z0, 1e-6);
      CHECK_CLOSE_ABS(axis.min, z0, 1e-6);
      BOOST_CHECK_EQUAL(axis.bType, equidistant);
    }

    // z row where elements are rotated around y
    Transform3D tr = Transform3D::Identity();
    tr.rotate(AngleAxis3D(M_PI / 4., Vector3D(0, 0, 1)));
    surfaces    = straightLineSurfaces(10, 3, Vector3D(0, 0, 0 + 1.5), tr);
    surfacesRaw = unpack_shared_vector(surfaces);
    pl          = ProtoLayer(surfacesRaw);
    trf         = Transform3D::Identity();
    axis = createEquidistantAxis(surfacesRaw, BinningValue::binZ, pl, trf);
    draw_surfaces(surfaces,
                  "SurfaceArrayCreator_createEquidistantAxis_Z_3.obj");
    BOOST_CHECK_EQUAL(axis.nBins, 10);
    CHECK_CLOSE_ABS(axis.max, 30.9749, 1e-3);
    CHECK_CLOSE_ABS(axis.min, -0.974873, 1e-3);
    BOOST_CHECK_EQUAL(axis.bType, equidistant);
  }

  BOOST_FIXTURE_TEST_CASE(SurfaceArrayCreator_createEquidistantAxis_R,
                          SurfaceArrayCreatorFixture)
  {

    // single element in r
    auto surfaces    = fullPhiTestSurfacesEC(1, 0, 0, 15);
    auto surfacesRaw = unpack_shared_vector(surfaces);
    draw_surfaces(surfaces,
                  "SurfaceArrayCreator_createEquidistantAxis_R_1.obj");
    auto       trf = Transform3D::Identity();
    ProtoLayer pl  = ProtoLayer(surfacesRaw);
    auto axis = createEquidistantAxis(surfacesRaw, BinningValue::binR, pl, trf);
    BOOST_CHECK_EQUAL(axis.nBins, 1);
    CHECK_CLOSE_ABS(axis.max, perp(Vector3D(17, 1, 0)), 1e-3);
    CHECK_CLOSE_ABS(axis.min, 13, 1e-3);
    BOOST_CHECK_EQUAL(axis.bType, equidistant);

    // multiple rings
    surfaces.resize(0);
    auto ringa = fullPhiTestSurfacesEC(30, 0, 0, 10);
    surfaces.insert(surfaces.end(), ringa.begin(), ringa.end());
    auto ringb = fullPhiTestSurfacesEC(30, 0, 0, 15);
    surfaces.insert(surfaces.end(), ringb.begin(), ringb.end());
    auto ringc = fullPhiTestSurfacesEC(30, 0, 0, 20);
    surfaces.insert(surfaces.end(), ringc.begin(), ringc.end());
    draw_surfaces(surfaces,
                  "SurfaceArrayCreator_createEquidistantAxis_R_2.obj");

    surfacesRaw = unpack_shared_vector(surfaces);
    pl          = ProtoLayer(surfacesRaw);
    trf         = Transform3D::Identity();
    axis = createEquidistantAxis(surfacesRaw, BinningValue::binR, pl, trf);

    BOOST_CHECK_EQUAL(axis.nBins, 3);
    CHECK_CLOSE_REL(axis.max, perp(Vector3D(20 + 2, 1, 0)), 1e-3);
    CHECK_CLOSE_ABS(axis.min, 8, 1e-3);
    BOOST_CHECK_EQUAL(axis.bType, equidistant);
  }

  // if there are concentring disc or barrel modules, the bin count might be off
  // we want to create _as few bins_ as possible, meaning the r-ring with
  // the lowest number of surfaces should be used for the bin count or
  // as basis for the variable edge procedure
  // double filling will make sure no surfaces are dropped
  BOOST_FIXTURE_TEST_CASE(SurfaceArrayCreator_dependentBinCounts,
                          SurfaceArrayCreatorFixture)
  {
    auto ringA = fullPhiTestSurfacesEC(10, 0, 0, 10, 2, 3);
    auto ringB = fullPhiTestSurfacesEC(15, 0, 0, 15, 2, 3.5);
    auto ringC = fullPhiTestSurfacesEC(20, 0, 0, 20, 2, 3.8);

    std::vector<std::shared_ptr<const Surface>> surfaces;
    std::copy(ringA.begin(), ringA.end(), std::back_inserter(surfaces));
    std::copy(ringB.begin(), ringB.end(), std::back_inserter(surfaces));
    std::copy(ringC.begin(), ringC.end(), std::back_inserter(surfaces));
    draw_surfaces(surfaces, "SurfaceArrayCreator_dependentBinCounts.obj");

    std::unique_ptr<SurfaceArray> sArray
        = m_SAC.surfaceArrayOnDisc(surfaces, equidistant, equidistant);
    std::cout << (*sArray) << std::endl;
    auto axes = sArray->getAxes();
    BOOST_CHECK_EQUAL(axes.at(0)->getNBins(), 3);
    BOOST_CHECK_EQUAL(axes.at(1)->getNBins(), 10);
  }

  BOOST_FIXTURE_TEST_CASE(SurfaceArrayCreator_completeBinning,
                          SurfaceArrayCreatorFixture)
  {
    SrfVec                      brl    = makeBarrel(30, 7, 2, 1);
    std::vector<const Surface*> brlRaw = unpack_shared_vector(brl);
    draw_surfaces(brl, "SurfaceArrayCreator_completeBinning_BRL.obj");

    detail::Axis<detail::AxisType::Equidistant,
                 detail::AxisBoundaryType::Closed>
        phiAxis(-M_PI, M_PI, 30u);
    detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Bound>
        zAxis(-14, 14, 7u);

    double R           = 10.;
    auto globalToLocal = [](const Vector3D& pos) {
      return Vector2D(phi(pos) + 2 * M_PI / 30 / 2, pos.z());
    };
    auto localToGlobal = [R](const Vector2D& loc) {
      double phi = loc[0] - 2 * M_PI / 30 / 2;
      return Vector3D(R * std::cos(phi), R * std::sin(phi), loc[1]);
    };

    auto sl
        = std::make_unique<SurfaceArray::SurfaceGridLookup<decltype(phiAxis),
                                                           decltype(zAxis)>>(
            globalToLocal,
            localToGlobal,
            std::make_tuple(std::move(phiAxis), std::move(zAxis)));
    sl->fill(brlRaw);
    SurfaceArray sa(std::move(sl), brl);

    // actually filled SA
    for (const auto& srf : brl) {
      Vector3D ctr        = srf->binningPosition(binR);
      auto     binContent = sa.at(ctr);

      BOOST_CHECK_EQUAL(binContent.size(), 1);
      BOOST_CHECK_EQUAL(srf.get(), binContent.at(0));
    }
  }

  BOOST_FIXTURE_TEST_CASE(SurfaceArrayCreator_barrelStagger,
                          SurfaceArrayCreatorFixture)
  {

    auto                        barrel = makeBarrelStagger(30, 7, 0, M_PI / 9.);
    auto                        brl    = barrel.first;
    std::vector<const Surface*> brlRaw = unpack_shared_vector(brl);
    draw_surfaces(brl, "SurfaceArrayCreator_barrelStagger.obj");

    ProtoLayer pl(brl);

    // EQUIDISTANT
    Transform3D tr = Transform3D::Identity();

    auto pAxisPhi = createEquidistantAxis(brlRaw, BinningValue::binPhi, pl, tr);
    auto pAxisZ   = createEquidistantAxis(brlRaw, BinningValue::binZ, pl, tr);

    double      R   = 10.;
    Transform3D itr = tr.inverse();

    auto globalToLocal = [tr](const Vector3D& pos) {
      Vector3D rot = tr * pos;
      return Vector2D(phi(rot), rot.z());
    };
    auto localToGlobal = [R, itr](const Vector2D& loc) {
      return itr * Vector3D(R * std::cos(loc[0]), R * std::sin(loc[0]), loc[1]);
    };

    auto sl = makeSurfaceGridLookup2D<detail::AxisBoundaryType::Closed,
                                      detail::AxisBoundaryType::Bound>(
        globalToLocal, localToGlobal, pAxisPhi, pAxisZ);

    sl->fill(brlRaw);
    SurfaceArray sa(std::move(sl), brl);
    auto         axes = sa.getAxes();
    BOOST_CHECK_EQUAL(axes.at(0)->getNBins(), 30);
    BOOST_CHECK_EQUAL(axes.at(1)->getNBins(), 7);
    std::cout << sa << std::endl;

    for (const auto& pr : barrel.second) {
      auto A = pr.first;
      auto B = pr.second;

      Vector3D ctr        = A->binningPosition(binR);
      auto     binContent = sa.at(ctr);
      BOOST_CHECK_EQUAL(binContent.size(), 2);
      std::set<const Surface*> act;
      act.insert(binContent[0]);
      act.insert(binContent[1]);

      std::set<const Surface*> exp;
      exp.insert(A);
      exp.insert(B);

      BOOST_CHECK(act == exp);
    }

    // VARIABLE
    BOOST_TEST_CONTEXT("Barrel Stagger Variable binning")
    {
      tr = Transform3D::Identity();

      auto pAxisPhiVar
          = createVariableAxis(brlRaw, BinningValue::binPhi, pl, tr);
      auto pAxisZVar = createVariableAxis(brlRaw, BinningValue::binZ, pl, tr);

      itr = tr.inverse();

      auto globalToLocalVar = [tr](const Vector3D& pos) {
        Vector3D rot = tr * pos;
        return Vector2D(phi(rot), rot.z());
      };
      auto localToGlobalVar = [R, itr](const Vector2D& loc) {
        return itr
            * Vector3D(R * std::cos(loc[0]), R * std::sin(loc[0]), loc[1]);
      };

      auto sl2 = makeSurfaceGridLookup2D<detail::AxisBoundaryType::Closed,
                                         detail::AxisBoundaryType::Bound>(
          globalToLocalVar, localToGlobalVar, pAxisPhiVar, pAxisZVar);

      sl2->fill(brlRaw);
      SurfaceArray sa2(std::move(sl2), brl);
      axes = sa2.getAxes();
      BOOST_CHECK_EQUAL(axes.at(0)->getNBins(), 30);
      BOOST_CHECK_EQUAL(axes.at(1)->getNBins(), 7);

      // check bin edges
      std::vector<double> phiEdgesExp
          = {-3.14159,  -2.93215,  -2.72271, -2.51327,    -2.30383, -2.0944,
             -1.88496,  -1.67552,  -1.46608, -1.25664,    -1.0472,  -0.837758,
             -0.628319, -0.418879, -0.20944, 4.44089e-16, 0.20944,  0.418879,
             0.628319,  0.837758,  1.0472,   1.25664,     1.46608,  1.67552,
             1.88496,   2.0944,    2.30383,  2.51327,     2.72271,  3.00831,
             3.14159};
      std::vector<double> zEdgesExp = {-14, -10, -6, -2, 2, 6, 10, 14};
      size_t              i         = 0;
      for (const auto& edge : axes.at(0)->getBinEdges()) {
        BOOST_TEST_INFO("phi edge index " << i);
        auto phiEdge = phiEdgesExp.at(i);
        CHECK_CLOSE_ABS(edge, phiEdge, 1e-5 * M_PI);
        i++;
      }
      i = 0;
      for (const auto& edge : axes.at(1)->getBinEdges()) {
        BOOST_TEST_INFO("z edge index " << i);
        CHECK_CLOSE_ABS(edge, zEdgesExp.at(i), 1e-6);
        i++;
      }

      std::cout << sa2 << std::endl;

      for (const auto& pr : barrel.second) {
        auto A = pr.first;
        auto B = pr.second;

        Vector3D ctr        = A->binningPosition(binR);
        auto     binContent = sa2.at(ctr);
        BOOST_CHECK_EQUAL(binContent.size(), 2);
        std::set<const Surface*> act;
        act.insert(binContent[0]);
        act.insert(binContent[1]);

        std::set<const Surface*> exp;
        exp.insert(A);
        exp.insert(B);

        BOOST_CHECK(act == exp);
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test

}  // namespace Acts
