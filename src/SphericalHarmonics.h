//=============================================================================
// Copyright 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/Types.h>

using namespace pmp;

//=============================================================================

// spherical harmonic function at a point p for a band l and its range m
// see: http://silviojemma.com/public/papers/lighting/spherical-harmonic-lighting.pdf
double spherical_harmonic(Point p, int l, int m);

//=============================================================================
