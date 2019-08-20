// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GeometryStatics.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Utilities/Definitions.hpp"

/// Define statics for Geometry in Tracking
///
namespace Acts {

// transformations

static const Transform3D s_idTransform
    = Transform3D::Identity();  //!< idendity transformation
static const Rotation3D s_idRotation
    = Rotation3D::Identity();  //!< idendity rotation

// axis system
static const Vector3D s_xAxis(1, 0, 0);  //!< global x Axis;
static const Vector3D s_yAxis(0, 1, 0);  //!< global y Axis;
static const Vector3D s_zAxis(0, 0, 1);  //!< global z Axis;

// unit vectors
static const Vector2D s_origin2D(0., 0.);

// origin

static const Vector3D s_origin(0, 0, 0);  //!< origin position

static const double helper[9] = {0., 1., 0., 1., 0., 0., 0., 0., -1.};

static const RotationMatrix3D s_idRotationZinverse(helper);
}