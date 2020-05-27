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
    GeodesicsInHeat(pmp::SurfaceMesh &mesh, bool geodist, bool euklid);

    ~GeodesicsInHeat();

    void getDistance(const int vertex, Eigen::VectorXd &dist,
                     Eigen::VectorXd &orthodist);

    void distance_to_texture_coordinates() const;

    void compute_geodesics(bool lumped = true);

private:
    SurfaceMesh &mesh_;

    Eigen::MatrixXd pos;

    Eigen::SparseMatrix<double> divOperator, gradOperator;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholL, cholA;

    double averageEdgeLength(const pmp::SurfaceMesh &mesh);

    void geodist_intrinsic_delauney(const int vertex, Eigen::VectorXd &dist);

    void buildDivOperator();

    void buildGradientOperator();

    double great_circle_distance(Vertex v, Vertex vv, double r = 1.0);

    double haversine_distance(Vertex v, Vertex vv, double r = 1.0);

    double vincenty_distance(Vertex v, Vertex vv, double r = 1.0);

    bool geodist_sphere_, geodist_cube_;
};
