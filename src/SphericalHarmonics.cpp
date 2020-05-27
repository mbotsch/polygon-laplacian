//=============================================================================
// Copyright 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "SphericalHarmonics.h"
#include "PolyDiffGeo.h"

//=============================================================================

using namespace pmp;

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

double factorial(int n)
{
    if (n == 0)
        return 1.0;
    return (double)n * factorial(n - 1);
}

//----------------------------------------------------------------------------

double scale(int l, int m)
{
    double temp =
        ((2.0 * l + 1.0) * factorial(l - m)) / (4.0 * M_PI * factorial(l + m));
    return sqrt(temp);
}

//----------------------------------------------------------------------------

double P(int l, int m, double x)
{
    // evaluate an Associated Legendre Polynomial P(l,m,x) at x
    double pmm = 1.0;
    if (m > 0)
    {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int i = 1; i <= m; i++)
        {
            pmm *= (-fact) * somx2;
            fact += 2.0;
        }
    }
    if (l == m)
        return pmm;
    double pmmp1 = x * (2.0 * m + 1.0) * pmm;
    if (l == m + 1)
        return pmmp1;
    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ++ll)
    {
        pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

//----------------------------------------------------------------------------

double spherical_harmonic(Point p, int l, int m)
{
    // l is the band, range [0..n]
    // m in the range [-l..l]
    // transform cartesian to spherical coordinates, assuming r = 1

    double phi = atan2(p[0], p[2]) + M_PI;
    double theta = acos(p[1] / norm(p));

    // return a point sample of a Spherical Harmonic basis function
    const double sqrt2 = sqrt(2.0);
    if (m == 0)
        return scale(l, 0) * P(l, m, cos(theta));
    else if (m > 0)
        return sqrt2 * scale(l, m) * cos(m * phi) * P(l, m, cos(theta));
    else
        return sqrt2 * scale(l, -m) * sin(-m * phi) * P(l, -m, cos(theta));
}

//=============================================================================
