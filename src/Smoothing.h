//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include <Eigen/Sparse>

//=============================================================================

using namespace pmp;

//=============================================================================

class Smoothing 
{
public:

    Smoothing(SurfaceMesh &mesh)
        : mesh_(mesh), vertices_(0), faces_(0), clamp_(false) {}

    //! Perform implicit Laplacian smoothing with \c timestep.
    void implicit_smoothing(Scalar timestep);

private:
    void update_stiffness_matrix();

private:
    SurfaceMesh &mesh_;
    Eigen::SparseMatrix<double> S_;
    unsigned int vertices_, faces_;
    bool clamp_;
};

//=============================================================================
