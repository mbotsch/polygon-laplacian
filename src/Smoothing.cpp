//=============================================================================
// Copyright 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "Smoothing.h"
#include "PolyDiffGeo.h"

//=============================================================================

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

void Smoothing::implicit_smoothing(Scalar timestep)
{
    if (!mesh_.n_vertices())
        return;

    // properties
    auto points = mesh_.vertex_property<Point>("v:point");

    // update stiffness matrix (if required)
    update_stiffness_matrix();

    // store center and area
    Point center_before = area_weighted_centroid(mesh_);
    Scalar area_before = mesh_area(mesh_);

    for (auto v : mesh_.vertices())
    {
        mesh_.position(v) -= center_before;
    }
    for (auto v : mesh_.vertices())
    {
        mesh_.position(v) *= sqrt(1.0 / area_before);
    }

    const unsigned int nv = mesh_.n_vertices();

    unsigned k = 0;
    SparseMatrix L, I(nv, nv), M;
    Eigen::MatrixXd B(nv, 3);

    I.setIdentity();

    for (auto v : mesh_.vertices())
    {
        B(k, 0) = points[v][0];
        B(k, 1) = points[v][1];
        B(k, 2) = points[v][2];
        k++;
    }

    setup_mass_matrix(mesh_, M);
    L = M - timestep * S_;

    Eigen::MatrixXd M_B = M * B;
    Eigen::MatrixXd X;
    static Eigen::SimplicialLLT<SparseMatrix> solver;
    solver.compute(L);
    X = solver.solve(M_B);

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "SurfaceSmoothing: Could not solve linear system\n";
    }
    else
    {
        // copy solution
        k = 0;
        for (unsigned int i = 0; i < nv; ++i)
        {
            Vertex v(i);
            if (!mesh_.is_boundary(v))
            {
                points[v][0] = X(k, 0);
                points[v][1] = X(k, 1);
                points[v][2] = X(k, 2);
                k++;
            }
        }
    }

    // restore original surface area
    Scalar area_after = mesh_area(mesh_);
    Scalar scale = sqrt(1 / area_after);
    Point center_after = area_weighted_centroid(mesh_);
    for (auto v : mesh_.vertices())
    {
        mesh_.position(v) -= center_after;
        mesh_.position(v) *= scale;
    }
}

//-----------------------------------------------------------------------------

void Smoothing::update_stiffness_matrix()
{
    if (mesh_.n_faces() != faces_ || mesh_.n_vertices() != vertices_ ||
        clamp_ != clamp_cotan_)
    {
        vertices_ = mesh_.n_vertices();
        faces_ = mesh_.n_faces();
        clamp_ = clamp_cotan_;

        std::cout << "Stiffness matrix has been updated" << std::endl;
        setup_stiffness_matrix(mesh_, S_);
    }
}

//=============================================================================
