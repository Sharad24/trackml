// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <ostream>
#include <string>
#include <vector>

using read_series  = std::vector<int>;
using read_range   = std::vector<double>;
using read_strings = std::vector<std::string>;

namespace std {
namespace detail {
  namespace {

    /// Helper function to print multiple elements in a container
    template <typename Iterator>
    inline std::ostream&
    printVector(Iterator      begin,
                Iterator      end,
                const char*   separator,
                std::ostream& os)
    {
      for (auto it = begin; it != end; ++it) {
        if (it != begin) {
          os << separator;
        }
        os << *it;
      }
      return os;
    }

  }  // namespace
}  // namespace detail

inline std::ostream&
operator<<(std::ostream& os, const read_series& vec)
{
  return detail::printVector(vec.begin(), vec.end(), " ", os);
}

inline std::ostream&
operator<<(std::ostream& os, const read_range& vec)
{
  return detail::printVector(vec.begin(), vec.end(), " ", os);
}

inline std::ostream&
operator<<(std::ostream& os, const read_strings& vec)
{
  return detail::printVector(vec.begin(), vec.end(), " ", os);
}

}  // namespace std
