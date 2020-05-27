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

class Curvature
{
public:
    // construct with mesh
    Curvature(SurfaceMesh &mesh);

    // destructor (deletes vertex property "v:curv")
    ~Curvature();

    // compute per-vertex mean curvatures, store them in vertex
    // property "v:curv"
    void compute();

    // return (absolute) mean curvature of vertex v
    Scalar operator()(Vertex v) const
    {
        assert(curvatures_);
        return curvatures_[v];
    }

    // convert curvature values ("v:curv") to 1D texture coordinates
    void curvature_to_texture_coordinates() const;

private:
    SurfaceMesh &mesh_;
    VertexProperty<Scalar> curvatures_;
};

//=============================================================================
