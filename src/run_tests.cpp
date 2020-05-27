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

// to lump or not to lump the mass matrix?
bool LUMP_MASS_MATRIX = true;

//=============================================================================

using namespace pmp;

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

// normalize all vertex position to norm=1.
// used to turn an almost spherical mesh into a spherical mesh.
void normalize(SurfaceMesh &mesh)
{
    for (auto v : mesh.vertices())
    {
        mesh.position(v) = normalize(mesh.position(v));
    }
}

//----------------------------------------------------------------------------

// compute curvature on a spherical mesh, compare to ground truth (=1)
void test_curvatures(SurfaceMesh &mesh)
{
    Curvature curvature(mesh);
    curvature.compute();

    double rms = 0.0;
    for (auto v : mesh.vertices())
    {
        double c = curvature(v);
        rms += (c - 1.0) * (c - 1.0);
    }
    rms /= (double)mesh.n_vertices();
    rms = sqrt(rms);

    std::cout << "rms error = " << rms << std::endl;
}

//----------------------------------------------------------------------------

void test_geodesics(SurfaceMesh &mesh)
{
    Vertex source(0);
    const Point p = mesh.position(source);

    GeodesicsInHeat heat(mesh);
    heat.precompute();
    heat.compute_distance_from(source);

    double rms = 0.0;
    for (auto v : mesh.vertices())
    {
        double geodesic_dist = heat(v);
        double euclidean_dist = distance(p, mesh.position(v));
        rms += pow(geodesic_dist - euclidean_dist, 2.0);
    }
    rms /= (double)mesh.n_vertices();
    rms = sqrt(rms);

    std::cout << "rms error = " << rms << std::endl;
}

//----------------------------------------------------------------------------

void test_sphericalH(SurfaceMesh &mesh)
{
    SpectralProcessing analyzer_(mesh);
    analyzer_.analyze_sphericalHarmonics();
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

    // whether or not we lump the mass matrix
    lump_mass_matrix_ = LUMP_MASS_MATRIX;
    std::cout << "\nlump mass matrix = " << lump_mass_matrix_ << "\n\n";

    // compute mean curvature
    if (!mytest || mytest == test)
    {
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "Curvature: Hex sphere\n";
            mesh.read("../data/unit-sphere.off");
            dualize(mesh);
            normalize(mesh);
            test_curvatures(mesh);
        }
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "Curvature: Fine sphere\n";
            mesh.read("../data/unit-sphere.off");
            SurfaceSubdivision(mesh).catmull_clark();
            normalize(mesh);
            test_curvatures(mesh);
        }
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "Curvature: Regular sphere\n";
            mesh.read("../data/quad-sphere.off");
            normalize(mesh);
            test_curvatures(mesh);
        }
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "Curvature : Noisy sphere\n";
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
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "SH: Hex sphere\n";
            mesh.read("../data/unit-sphere.off");
            dualize(mesh);
            normalize(mesh);
            test_sphericalH(mesh);
        }
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "SH: Fine sphere\n";
            mesh.read("../data/unit-sphere.off");
            SurfaceSubdivision(mesh).catmull_clark();
            normalize(mesh);
            test_sphericalH(mesh);
        }
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "SH: Regular sphere\n";
            mesh.read("../data/quad-sphere.off");
            normalize(mesh);
            test_sphericalH(mesh);
        }
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "SH: Noisy sphere\n";
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
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "Geodist: Subdivided quad plane \n";
            mesh.read("../data/quad-plane.obj");
            SurfaceSubdivision(mesh).catmull_clark();
            test_geodesics(mesh);
        }
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "Geodist: Quad plane\n";
            mesh.read("../data/quad-plane.obj");
            test_geodesics(mesh);
        }
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "Geodist: L-plane\n";
            mesh.read("../data/L-plane.obj");
            test_geodesics(mesh);
        }
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "Geodist: Tetris-plane (non-starshaped)\n";
            mesh.read("../data/tetris.obj");
            test_geodesics(mesh);
        }
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "Geodist:  Tetris-plane (non-convex)\n";
            mesh.read("../data/tetris_2.obj");
            test_geodesics(mesh);
        }
    }
    ++test;

    // measure timings
    if (!mytest || mytest == test)
    {
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "Timings: Hexasphere\n";
            mesh.read("../data/Spheres/hexaSphere.off");
            timing_construction(mesh);
            timing_solution(mesh);
        }
        {
            std::cout
                << "----------------------------------------------------\n";
            std::cout << "Timings: Fine Sphere\n";
            mesh.read("../data/Spheres/fineSphere.off");
            timing_construction(mesh);
            timing_solution(mesh);
        }
    }
}

//=============================================================================
