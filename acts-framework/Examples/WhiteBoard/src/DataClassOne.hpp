// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DataClassOne_h
#define DataClassOne_h

#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace FWE {

class DataClassOne
{
public:
  DataClassOne(const std::string& stringData, size_t eventData)
    : m_dataString(stringData), m_dataSizeT(eventData)
  {
  }

  /// the contained data : string
  const std::string
  data() const;

private:
  std::string m_dataString;  ///< data member string
  size_t      m_dataSizeT;   ///< data member size_t
};

inline const std::string
DataClassOne::data() const
{
  std::ostringstream oss;
  oss << "Data : " << m_dataString << " | " << m_dataSizeT;
  return oss.str();
}

typedef std::vector<DataClassOne> DataClassOneCollection;

}  // namespace FWE

#endif
