//=============================================================================
#pragma once
//=============================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>

//=============================================================================

using namespace pmp;

//=============================================================================

extern bool clamp_cotan_;

//=============================================================================

void setup_stiffness_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &S);

void setup_mass_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &M, bool lumped = true);

void localStiffnessMatrix(const Eigen::MatrixXd &poly, const Eigen::Vector3d &min, Eigen::VectorXd &w,
                          Eigen::MatrixXd &L);

void localMassMatrix(const Eigen::MatrixXd &poly, const Eigen::Vector3d &min, Eigen::VectorXd &w,
                     Eigen::MatrixXd &M);

//!compute the transformation matrix for an arbitrary mesh that inserts a chosen point per face .
void setup_prolongation_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &A);

//! Computes a diagonal matrix where the entries are the sums of the rows.
void lump_matrix(Eigen::SparseMatrix<double> &D);

//----------------------------------Centroid and area computations------------------------------------------------------

//! barycenter/my_centroid of mesh, computed as area-weighted mean of vertices.
Point area_weighted_centroid(const SurfaceMesh &mesh);

//! computes the area of a face.
double face_area(const SurfaceMesh &mesh, Face f);

//! surface area of the polygon mesh.
Scalar polygon_surface_area(const SurfaceMesh &mesh);

//---------------------Gradient/ Divergence --------------------------------------------------------------

//! Computes the gradient on a triangle formed by the points i, j and k.
Eigen::Vector3d gradient_hat_function(Point i, Point j, Point k);

//! Computes the Gradient matrix for any given mesh.
void setup_Gradient_Matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &G);

//! Computes the Divergence matrix for any given mesh.
void setup_Divergence_Matrix(SurfaceMesh &mesh,
                             Eigen::SparseMatrix<double> &Gt);

//! Computes the diagonal mass matrix W needed for L = DWG.
void setup_Gradient_Mass_Matrix(SurfaceMesh &mesh,
                                Eigen::SparseMatrix<double> &M);

//! Computes the squared triangle area minimizing points and its affine combination weights
//! for each face and stores it in a prior defined property.
void setup_face_point_properties(SurfaceMesh &mesh);

//------------------------Weight computation -----------------------------------------------------------------


//! Computes the affine weights for the polygon vertices to form the implicit vertex.
void find_polygon_weights(const Eigen::MatrixXd &poly,
                          Eigen::VectorXd &weights);


//=============================================================================
