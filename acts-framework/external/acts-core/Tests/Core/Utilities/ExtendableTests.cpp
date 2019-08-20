// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Extendable Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/Utilities/detail/Extendable.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // This tests the implementation of the ActionList
  // and the standard aborters
  BOOST_AUTO_TEST_CASE(Extendable_)
  {
    struct TypeA
    {
      double vaA = 0.;
    };

    struct TypeB
    {
      int vaB = 0;
    };

    struct TypeC
    {
      char vaC = '0';
    };

    // Test the empty list
    detail::Extendable<> nullist;
    (void)nullist;
    BOOST_CHECK_EQUAL(std::tuple_size<std::tuple<>>::value, 0);

    detail::Extendable<TypeA> alist;
    auto&                     a0_object = alist.get<TypeA>();
    a0_object.vaA                       = 1.;
    BOOST_CHECK_EQUAL(alist.get<TypeA>().vaA, 1.);

    detail::Extendable<TypeA, TypeB> ablist;
    auto& a1_object = ablist.get<TypeA>();
    a1_object.vaA   = 2.;
    auto& b1_object = ablist.get<TypeB>();
    b1_object.vaB   = 3;
    BOOST_CHECK_EQUAL(ablist.get<TypeA>().vaA, 2.);
    BOOST_CHECK_EQUAL(ablist.get<TypeB>().vaB, 3);

    TypeC c;
    c.vaC = '4';
    detail::Extendable<TypeA, TypeB, TypeC> abcList = ablist.append<TypeC>(c);
    BOOST_CHECK_EQUAL(abcList.get<TypeA>().vaA, 2.);
    BOOST_CHECK_EQUAL(abcList.get<TypeB>().vaB, 3);
    BOOST_CHECK_EQUAL(abcList.get<TypeC>().vaC, '4');
  }

}  // namespace Test
}  // namespace Acts
