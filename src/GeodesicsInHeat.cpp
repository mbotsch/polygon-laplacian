//=============================================================================
// Copyright 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "GeodesicsInHeat.h"
#include "PolyDiffGeo.h"
#include <cfloat>

//=============================================================================

GeodesicsInHeat::GeodesicsInHeat(SurfaceMesh &mesh) : mesh_(mesh)
{
    mesh_.add_face_property<Point>("f:point");
    mesh_.add_face_property<Eigen::VectorXd>("f:weights");

    distance_ = mesh_.add_vertex_property<Scalar>("v:dist");

    setup_face_point_properties(mesh_);
}

//-----------------------------------------------------------------------------

GeodesicsInHeat::~GeodesicsInHeat()
{
    auto area_points = mesh_.face_property<Point>("f:point");
    if (area_points)
        mesh_.remove_face_property(area_points);

    auto area_weights = mesh_.face_property<Eigen::VectorXd>("f:weights");
    if (area_weights)
        mesh_.remove_face_property(area_weights);

    mesh_.remove_vertex_property(distance_);
}

//-----------------------------------------------------------------------------

double GeodesicsInHeat::avg_edge_length() const
{
    double avgLen = 0.;
    for (auto e : mesh_.edges())
        avgLen += mesh_.edge_length(e);
    return avgLen / (Scalar)mesh_.n_edges();
}

//-----------------------------------------------------------------------------

void GeodesicsInHeat::compute_distance_from(Vertex source)
{
    // diffuse heat
    const int N = mesh_.n_vertices();

    Eigen::SparseVector<double> b(N);
    Eigen::SparseMatrix<double> W;
    setup_gradient_mass_matrix(mesh_, W);
    b.coeffRef(source.idx()) = 1.;

    // compute gradients
    Eigen::VectorXd heat = cholA.solve(b);
    Eigen::VectorXd grad = gradient_ * heat;

    // normalize gradients
    for (int i = 0; i < grad.rows(); i += 3)
    {
        dvec3 &g = *reinterpret_cast<dvec3 *>(&grad[i]);
        double n = norm(g);
        if (n > DBL_MIN)
            g /= n;
    }

    Eigen::VectorXd dist = cholL.solve(divergence_ * (-grad));

    double mi = dist.minCoeff();
    for (int i = 0; i < dist.rows(); ++i)
        dist[i] -= mi;

    int k = 0;
    for (auto v : mesh_.vertices())
    {
        distance_[v] = dist[k];
        k++;
    }
}

//-----------------------------------------------------------------------------

void GeodesicsInHeat::precompute()
{
    Eigen::SparseMatrix<double> S, M, A, M_bar;

    // setup divergence and gradient
    setup_gradient_matrix(mesh_, gradient_);
    setup_divergence_matrix(mesh_, divergence_);

    // setup stiffness and mass matrix
    setup_stiffness_matrix(mesh_, S);
    setup_mass_matrix(mesh_, M);

    // setup heat flow matrix
    const double h = pow(avg_edge_length(), 2);
    A = M - h * S;

    // factorize stiffness matrix
    cholL.analyzePattern(S);
    cholL.factorize(S);

    // factorize heat flow matrix
    cholA.analyzePattern(A);
    cholA.factorize(A);
}

//-----------------------------------------------------------------------------

void GeodesicsInHeat::distance_to_texture_coordinates() const
{
    auto distances = mesh_.get_vertex_property<Scalar>("v:dist");
    assert(distances);

    // find maximum distance
    Scalar maxdist(0);
    for (auto v : mesh_.vertices())
    {
        if (distances[v] <= FLT_MAX)
        {
            maxdist = std::max(maxdist, distances[v]);
        }
    }

    auto tex = mesh_.vertex_property<TexCoord>("v:tex");
    for (auto v : mesh_.vertices())
    {
        if (distances[v] <= FLT_MAX)
        {
            tex[v] = TexCoord(distances[v] / maxdist, 0.0);
        }
        else
        {
            tex[v] = TexCoord(1.0, 0.0);
        }
    }

    // remove per-halfedge texture coordinates
    auto htex = mesh_.get_halfedge_property<TexCoord>("h:tex");
    if (htex)
        mesh_.remove_halfedge_property(htex);
}

//=============================================================================
