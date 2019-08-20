// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DataClassTwo_h
#define DataClassTwo_h

#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace FWE {

class DataClassTwo
{
public:
  DataClassTwo(const std::string& stringData, double eventData)
    : m_dataString(stringData), m_dataDouble(eventData)
  {
  }

  /// the contained data : string
  const std::string
  data() const;

private:
  std::string m_dataString;  /// data member string
  double      m_dataDouble;  /// data member size_t
};

inline const std::string
DataClassTwo::data() const
{
  std::ostringstream oss;
  oss << "Data : " << m_dataString << " | " << m_dataDouble;
  return oss.str();
}

typedef std::vector<DataClassTwo> DataClassTwoCollection;

}  // namespace FWE

#endif
