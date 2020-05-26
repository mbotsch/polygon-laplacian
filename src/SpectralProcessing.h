//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>

//=============================================================================

using namespace pmp;

//=============================================================================

class SpectralProcessing {
public:
    SpectralProcessing(SurfaceMesh &mesh) : mesh(mesh) {}

    /// renormalisation constant for SH function
    double scale(int l, int m);

    // spherical harmonic function at a point p for a band l and its range m
    // see: http://silviojemma.com/public/papers/lighting/spherical-harmonic-lighting.pdf
    double sphericalHarmonic(Point p, int l, int m);

    // evaluate an Associated Legendre Polynomial P(l,m,x) at x
    double P(int l, int m, double x);

    void analyze_sphericalHarmonics(bool lumped);

private:
    SurfaceMesh &mesh;
};

//=============================================================================
