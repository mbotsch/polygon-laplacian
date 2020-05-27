//=============================================================================
// Copyright 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>

//=============================================================================

using namespace pmp;

//=============================================================================

class Parameterization
{
public:

    // give a mesh in the constructor
    Parameterization(SurfaceMesh& mesh) : mesh_(mesh) {}

    // compute discrete harmonic parameterization
    bool harmonic();
    
    // compute discrete harmonic parameterization using free boundaries following
    // "Spectral Conformal Parametrization" [Mullen et al., 2008]
    bool harmonic_free_boundary();

private:
    // map boundary to unit circle
    bool setup_boundary_constraints();

private:
    SurfaceMesh& mesh_;
};

//=============================================================================
