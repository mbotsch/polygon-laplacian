//=============================================================================
// Copyright 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include <Eigen/Sparse>

//=============================================================================

using namespace pmp;

//=============================================================================

class GeodesicsInHeat
{
public:
    // construct with mesh
    GeodesicsInHeat(SurfaceMesh &mesh);

    // destructor: cleans up allocated properties
    ~GeodesicsInHeat();

    // setup and factorize matrices
    void precompute();

    // compute geodesic distance of all vertices from source vertex
    void compute_distance_from(Vertex source);

    // get previously computed distance of vertex v
    Scalar operator()(Vertex v) const
    {
        assert(distance_);
        return distance_[v];
    }

    // convert distances to texture coordinates for visualization
    void distance_to_texture_coordinates() const;

private:
    double avg_edge_length() const;

private:
    SurfaceMesh &mesh_;
    VertexProperty<Scalar> distance_;

    Eigen::SparseMatrix<double> divergence_, gradient_;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholL, cholA;
};
