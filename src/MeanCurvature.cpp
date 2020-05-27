//=============================================================================
// Copyright 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "MeanCurvature.h"
#include "PolyDiffGeo.h"

//=============================================================================

using namespace pmp;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

Curvature::Curvature(SurfaceMesh &mesh) : mesh_(mesh) {}

//-----------------------------------------------------------------------------

Curvature::~Curvature()
{
    if (curvatures_)
        mesh_.remove_vertex_property(curvatures_);
}

//-----------------------------------------------------------------------------

void Curvature::compute()
{
    if (!mesh_.n_vertices())
        return;

    // properties
    auto points = mesh_.vertex_property<Point>("v:point");
    curvatures_ = mesh_.add_vertex_property<Scalar>("v:curv");

    const unsigned int nv = mesh_.n_vertices();

    Eigen::SparseMatrix<double> M, S;
    Eigen::MatrixXd B(nv, 3), B_(mesh_.n_edges(), 3);
    Eigen::VectorXd H(nv), test(3);

    for (auto v : mesh_.vertices())
    {
        B.row(v.idx()) = static_cast<Eigen::Vector3d>(points[v]);
    }

    setup_stiffness_matrix(mesh_, S);
    setup_mass_matrix(mesh_, M);

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);
    Eigen::MatrixXd X = solver.solve(S * B);

    unsigned i = 0;
    for (auto v : mesh_.vertices())
    {
        curvatures_[v] = 0.5 * X.row(i).norm();
        i++;
    }
}

//-----------------------------------------------------------------------------

void Curvature::curvature_to_texture_coordinates() const
{
    if (!curvatures_)
        return;

    // sort curvature values
    std::vector<Scalar> values;
    values.reserve(mesh_.n_vertices());
    for (auto v : mesh_.vertices())
    {
        values.push_back(curvatures_[v]);
    }
    std::sort(values.begin(), values.end());
    unsigned int n = values.size() - 1;

    // clamp upper/lower 5%
    unsigned int i = n / 20;
    Scalar kmin = values[i];
    Scalar kmax = values[n - 1 - i];

    // generate 1D texture coordiantes
    auto tex = mesh_.vertex_property<TexCoord>("v:tex");
    if (kmin < 0.0) // signed
    {
        kmax = std::max(fabs(kmin), fabs(kmax));
        for (auto v : mesh_.vertices())
        {
            tex[v] = TexCoord((0.5f * curvatures_[v] / kmax) + 0.5f, 0.0);
        }
    }
    else // unsigned
    {
        for (auto v : mesh_.vertices())
        {
            tex[v] = TexCoord((curvatures_[v] - kmin) / (kmax - kmin), 0.0);
        }
    }

    // remove per-halfedge texture coordinates
    auto htex = mesh_.get_halfedge_property<TexCoord>("h:tex");
    if (htex)
        mesh_.remove_halfedge_property(htex);
}

//=============================================================================
