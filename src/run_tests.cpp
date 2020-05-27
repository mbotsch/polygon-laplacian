//=============================================================================
// Copyright 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include <pmp/algorithms/DifferentialGeometry.h>
#include <pmp/algorithms/SurfaceSubdivision.h>
#include <pmp/Timer.h>

#include "MeanCurvature.h"
#include "SpectralProcessing.h"
#include "GeodesicsInHeat.h"
#include "PolyDiffGeo.h"

//=============================================================================

using namespace pmp;

//=============================================================================

enum MassMatrix
{
    lumped = true,
    unlumped = false
};

enum ObjectShape
{
    Sphere = 1,
    Plane = 2
};

//=============================================================================

void dualize(SurfaceMesh &mesh)
{
    SurfaceMesh dual;

    auto fvertex = mesh.add_face_property<Vertex>("f:vertex");
    for (auto f : mesh.faces())
    {
        fvertex[f] = dual.add_vertex(centroid(mesh, f));
    }

    for (auto v : mesh.vertices())
    {
        if (!mesh.is_boundary(v))
        {
            std::vector<Vertex> vertices;
            for (auto f : mesh.faces(v))
                vertices.push_back(fvertex[f]);
            dual.add_face(vertices);
        }
    }

    mesh = dual;
}

//----------------------------------------------------------------------------

void normalize(SurfaceMesh &mesh)
{
    for (auto v : mesh.vertices())
    {
        mesh.position(v) = normalize(mesh.position(v));
    }
}

//----------------------------------------------------------------------------

void test_curvatures(SurfaceMesh &mesh)
{
    Curvature analyzer(mesh, true);

    std::cout << "lumped mass matrix: " << std::endl;
    analyzer.analyze_curvature(lumped);
    std::cout << "un-lumped mass matrix: " << std::endl;
    analyzer.analyze_curvature(unlumped);
}

//----------------------------------------------------------------------------

void test_geodesics(SurfaceMesh &mesh, ObjectShape shape)
{
    Eigen::VectorXd dist, geodist;

    bool sphere = (shape == Sphere);
    bool plane = (shape == Plane);

    std::cout << "lumped mass matrix: " << std::endl;
    {
        GeodesicsInHeat heat(mesh, sphere, plane);
        heat.compute_geodesics(lumped);
        heat.getDistance(0, dist, geodist);
    }

    std::cout << "un-lumped mass matrix: " << std::endl;
    {
        GeodesicsInHeat heat(mesh, sphere, plane);
        heat.compute_geodesics(unlumped);
        heat.getDistance(0, dist, geodist);
    }
}

//----------------------------------------------------------------------------

void test_sphericalH(SurfaceMesh &mesh)
{
    std::cout << "lumped mass matrix: " << std::endl;
    {
        SpectralProcessing analyzer_(mesh);
        analyzer_.analyze_sphericalHarmonics(lumped);
    }

    std::cout << "un-lumped mass matrix: " << std::endl;
    {
        SpectralProcessing analyzer_(mesh);
        analyzer_.analyze_sphericalHarmonics(unlumped);
    }
}

//----------------------------------------------------------------------------

void timing_construction(SurfaceMesh mesh)
{
    Eigen::SparseMatrix<double> S;

    pmp::Timer t;
    t.start();

    for (int i = 0; i < 20; i++)
        setup_stiffness_matrix(mesh, S);

    t.stop();
    std::cout << "Time for stiffness matrix construction : "
              << t.elapsed() / 20.0 << "ms" << std::endl;
}

//----------------------------------------------------------------------------

void timing_solution(SurfaceMesh mesh)
{
    auto points = mesh.vertex_property<Point>("v:point");
    Eigen::SparseMatrix<double> S, M, A;

    setup_stiffness_matrix(mesh, S);
    setup_mass_matrix(mesh, M);

    Eigen::MatrixXd B(mesh.n_vertices(), 3);

    for (auto v : mesh.vertices())
    {
        B(v.idx(), 0) = points[v][0];
        B(v.idx(), 1) = points[v][1];
        B(v.idx(), 2) = points[v][2];
    }

    A = M - 0.1 * S;
    Eigen::MatrixXd M_B = M * B;

    static Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;

    pmp::Timer t;

    t.start();
    for (int i = 0; i < 10; i++)
    {
        solver.compute(A);
    }
    t.stop();
    std::cout << "Time to factorize system matrix A : " << t.elapsed() / 10.0
              << "ms" << std::endl;

    t.start();
    for (int i = 0; i < 10; i++)
        Eigen::MatrixXd X = solver.solve(M_B);
    t.stop();
    std::cout << "Time to solve: " << t.elapsed() / 10.0 << "ms" << std::endl;
}

//----------------------------------------------------------------------------

int main(int argc, char **argv)
{
    SurfaceMesh mesh;

    // which test to run
    int mytest = (argc > 1) ? atoi(argv[1]) : 0;

    // counter for tests
    int test = 1;

    // compute mean curvature
    if (!mytest || mytest == test)
    {
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "Curvature: Hex sphere\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/unit-sphere.off");
            dualize(mesh);
            normalize(mesh);
            test_curvatures(mesh);
        }
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "Curvature: Fine sphere\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/unit-sphere.off");
            SurfaceSubdivision(mesh).catmull_clark();
            normalize(mesh);
            test_curvatures(mesh);
        }
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "Curvature: Regular sphere\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/quad-sphere.off");
            normalize(mesh);
            test_curvatures(mesh);
        }
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "Curvature : Noisy sphere\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/noisy-sphere.off");
            normalize(mesh);
            test_curvatures(mesh);
        }
    }
    ++test;

    // compute spherical harmonics
    if (!mytest || mytest == test)
    {
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "SH: Hex sphere\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/unit-sphere.off");
            dualize(mesh);
            normalize(mesh);
            test_sphericalH(mesh);
        }
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "SH: Fine sphere\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/unit-sphere.off");
            SurfaceSubdivision(mesh).catmull_clark();
            normalize(mesh);
            test_sphericalH(mesh);
        }
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "SH: Regular sphere\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/quad-sphere.off");
            normalize(mesh);
            test_sphericalH(mesh);
        }
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "SH: Noisy sphere\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/noisy-sphere.off");
            normalize(mesh);
            test_sphericalH(mesh);
        }
    }
    ++test;

    // compute geodesic distance
    if (!mytest || mytest == test)
    {
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "Geodist: Subdivided quad plane \n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/quad-plane.obj");
            SurfaceSubdivision(mesh).catmull_clark();
            test_geodesics(mesh, Plane);
        }
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "Geodist: Quad plane\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/quad-plane.obj");
            test_geodesics(mesh, Plane);
        }
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "Geodist: L-plane\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/L-plane.obj");
            test_geodesics(mesh, Plane);
        }
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "Geodist: Tetris-plane (non-starshaped)\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/tetris.obj");
            test_geodesics(mesh, Plane);
        }
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "Geodist:  Tetris-plane (non-convex)\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/tetris_2.obj");
            test_geodesics(mesh, Plane);
        }
    }

    // measure timings
    if (!mytest || mytest == test)
    {
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "Timings: Hexasphere\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/Spheres/hexaSphere.off");
            timing_construction(mesh);
            timing_solution(mesh);
        }
        {
            std::cout << "----------------------------------------------------"
                      << std::endl;
            std::cout << "Timings: Fine Sphere\n";
            std::cout << "----------------------------------------------------"
                      << std::endl;

            mesh.read("../data/Spheres/fineSphere.off");
            timing_construction(mesh);
            timing_solution(mesh);
        }
    }
}

//=============================================================================
