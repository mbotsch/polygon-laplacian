//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>

//=============================================================================

using namespace pmp;

//=============================================================================

class Curvature {
public:
    Curvature(SurfaceMesh &mesh, bool compare)
            : mesh_(mesh), compare_to_sphere(compare) {
    }

    //! Visualizes the mean curvature of our mesh.
    void analyze_curvature(bool lumped = true);


private:
    SurfaceMesh &mesh_;
    bool compare_to_sphere;

    //! convert curvature values ("v:curv") to 1D texture coordinates
    void curvature_to_texture_coordinates() const;
};

//=============================================================================
