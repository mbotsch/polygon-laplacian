//=============================================================================
// Copyright 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "PolyDiffGeo.h"
#include <pmp/algorithms/DifferentialGeometry.h>
#include <pmp/algorithms/SurfaceTriangulation.h>

//=============================================================================

using namespace std;
using namespace pmp;

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

const double eps = 1e-10;
bool clamp_cotan_ = false;
bool lump_mass_matrix_ = true;

//=============================================================================

void setup_stiffness_matrix(const SurfaceMesh &mesh, Eigen::SparseMatrix<double> &S)
{
    const int nv = mesh.n_vertices();

    std::vector<Vertex> vertices; // polygon vertices
    Eigen::MatrixXd polygon;      // positions of polygon vertices
    Eigen::VectorXd weights;      // affine weights of virtual vertex
    Eigen::MatrixXd Si;           // local stiffness matrix

    std::vector<Eigen::Triplet<double>> trip;

    for (Face f : mesh.faces())
    {
        // collect polygon vertices
        vertices.clear();
        for (Vertex v : mesh.vertices(f))
        {
            vertices.push_back(v);
        }
        const int n=vertices.size();

        // collect their positions
        polygon.resize(n, 3);
        for (int i=0; i<n; ++i)
        {
            polygon.row(i) = (Eigen::Vector3d) mesh.position(vertices[i]);
        }

        // compute virtual vertex, setup local stiffness matrix
        compute_virtual_vertex(polygon, weights);
        setup_polygon_stiffness_matrix(polygon, weights, Si);

        // assemble to global stiffness matrix
        for (int j=0; j<n; ++j)
        {
            for (int k=0; k<n; ++k)
            {
                trip.emplace_back(vertices[k].idx(), vertices[j].idx(), -Si(k, j));
            }
        }
    }

    // build sparse matrix from triplets
    S.resize(nv, nv);
    S.setFromTriplets(trip.begin(), trip.end());
}

//----------------------------------------------------------------------------------

void setup_polygon_stiffness_matrix(const Eigen::MatrixXd &polygon,
                                    const Eigen::VectorXd &vweights, 
                                    Eigen::MatrixXd &S)
{
    const int n = (int)polygon.rows();
    S.resize(n, n);
    S.setZero();

    // compute position of virtual vertex
    Eigen::Vector3d vvertex = polygon.transpose() * vweights;

    Eigen::VectorXd ln(n + 1);
    ln.setZero();

    double l[3], l2[3];

    for (int i = 0; i < n; ++i)
    {
        const int i1 = (i + 1) % n;

        l2[2] = (polygon.row(i) - polygon.row(i1)).squaredNorm();
        l2[0] = (polygon.row(i1) - vvertex.transpose()).squaredNorm();
        l2[1] = (polygon.row(i) - vvertex.transpose()).squaredNorm();

        l[0] = sqrt(l2[0]);
        l[1] = sqrt(l2[1]);
        l[2] = sqrt(l2[2]);

        const double arg = (l[0] + (l[1] + l[2])) * (l[2] - (l[0] - l[1])) *
                           (l[2] + (l[0] - l[1])) * (l[0] + (l[1] - l[2]));
        const double area = 0.5 * sqrt(arg);
        if (area > 1e-7)
        {
            l[0] = 0.25 * (l2[1] + l2[2] - l2[0]) / area;
            l[1] = 0.25 * (l2[2] + l2[0] - l2[1]) / area;
            l[2] = 0.25 * (l2[0] + l2[1] - l2[2]) / area;

            S(i1, i1) += l[0];
            S(i, i) += l[1];
            S(i1, i) -= l[2];
            S(i, i1) -= l[2];
            S(i, i) += l[2];
            S(i1, i1) += l[2];

            ln(i1) -= l[0];
            ln(i) -= l[1];
            ln(n) += l[0] + l[1];
        }
    }

    // sandwiching with (local) restriction and prolongation matrix
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            S(i, j) += vweights(i) * ln(j) + vweights(j) * ln(i) + vweights(i) * vweights(j) * ln(n);
}

//----------------------------------------------------------------------------------

void setup_mass_matrix(const SurfaceMesh &mesh, Eigen::SparseMatrix<double> &M)
{
    const int nv = mesh.n_vertices();

    std::vector<Vertex> vertices; // polygon vertices
    Eigen::MatrixXd polygon;      // positions of polygon vertices
    Eigen::VectorXd weights;      // affine weights of virtual vertex
    Eigen::MatrixXd Mi;           // local mass matrix

    std::vector<Eigen::Triplet<double>> trip;

    for (Face f : mesh.faces())
    {
        // collect polygon vertices
        vertices.clear();
        for (Vertex v : mesh.vertices(f))
        {
            vertices.push_back(v);
        }
        const int n=vertices.size();

        // collect their positions
        polygon.resize(n, 3);
        for (int i=0; i<n; ++i)
        {
            polygon.row(i) = (Eigen::Vector3d) mesh.position(vertices[i]);
        }

        // compute virtual vertex, setup local mass matrix
        compute_virtual_vertex(polygon, weights);
        setup_polygon_mass_matrix(polygon, weights, Mi);

        // assemble into global mass matrix
        for (int j=0; j<n; ++j)
        {
            for (int k=0; k<n; ++k)
            {
                trip.emplace_back(vertices[k].idx(), vertices[j].idx(), Mi(k, j));
            }
        }
    }

    // build sparse matrix from triplets
    M.resize(nv, nv);
    M.setFromTriplets(trip.begin(), trip.end());

    // optional: lump mass matrix
    if (lump_mass_matrix_)
    {
        lump_matrix(M);
    }
}

//----------------------------------------------------------------------------------

void setup_polygon_mass_matrix(const Eigen::MatrixXd &polygon,
                               const Eigen::VectorXd &vweights,
                               Eigen::MatrixXd &M)
{
    const int n = (int)polygon.rows();
    M.resize(n, n);
    M.setZero();

    // compute position of virtual vertex
    Eigen::Vector3d vvertex = polygon.transpose() * vweights;

    Eigen::VectorXd ln(n + 1);
    ln.setZero();

    double l[3], l2[3];

    for (int i = 0; i < n; ++i)
    {
        const int i1 = (i + 1) % n;

        l2[2] = (polygon.row(i) - polygon.row(i1)).squaredNorm();
        l2[0] = (polygon.row(i1) - vvertex.transpose()).squaredNorm();
        l2[1] = (polygon.row(i) - vvertex.transpose()).squaredNorm();

        l[0] = sqrt(l2[0]);
        l[1] = sqrt(l2[1]);
        l[2] = sqrt(l2[2]);

        const double arg = (l[0] + (l[1] + l[2])) * (l[2] - (l[0] - l[1])) *
                           (l[2] + (l[0] - l[1])) * (l[0] + (l[1] - l[2]));
        const double area = 0.25 * sqrt(arg);

        l[0] = 1.0 / 6.0 * area;
        l[1] = 1.0 / 12.0 * area;

        M(i1, i1) += 1.0 / 6.0 * area;
        M(i, i) += 1.0 / 6.0 * area;
        M(i1, i) += 1.0 / 12.0 * area;
        M(i, i1) += 1.0 / 12.0 * area;

        ln(i1) += l[1];
        ln(i) += l[1];
        ln(n) += l[0];
    }

    // sandwiching with (local) restriction and prolongation matrix
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            M(i, j) += vweights(i) * ln(j) + vweights(j) * ln(i) + vweights(i) * vweights(j) * ln(n);
}

//-----------------------------------------------------------------------------

void lump_matrix(SparseMatrix &D)
{
    std::vector<Triplet> triplets;
    triplets.reserve(D.rows() * 6);

    for (int k = 0; k < D.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(D, k); it; ++it)
        {
            triplets.emplace_back(it.row(), it.row(), it.value());
        }
    }

    D.setFromTriplets(triplets.begin(), triplets.end());
}

//----------------------------------------------------------------------------------

void setup_prolongation_matrix(const SurfaceMesh &mesh, SparseMatrix &A)
{
    auto area_weights = mesh.get_face_property<Eigen::VectorXd>("f:weights");

    const unsigned int nv = mesh.n_vertices();
    const unsigned int nf = mesh.n_faces();

    std::vector<Triplet> tripletsA;

    for (auto v : mesh.vertices())
    {
        tripletsA.emplace_back(v.idx(), v.idx(), 1.0);
    }

    unsigned int j = 0;
    for (auto f : mesh.faces())
    {
        const Eigen::VectorXd& w = area_weights[f];

        unsigned int i = 0;
        for (auto v : mesh.vertices(f))
        {
            tripletsA.emplace_back(nv + j, v.idx(), w(i));
            i++;
        }
        j++;
    }

    // build sparse matrix from triplets
    A.resize(nv + nf, nv);
    A.setFromTriplets(tripletsA.begin(), tripletsA.end());
}

//-----------------------------------------------------------------------------

Scalar mesh_area(const SurfaceMesh &mesh)
{
    Scalar area(0);
    for (auto f : mesh.faces())
    {
        area += face_area(mesh, f);
    }
    return area;
}

//-----------------------------------------------------------------------------

double face_area(const SurfaceMesh &mesh, Face f)
{
#if 0
    // vector area of a polygon
    Normal n(0, 0, 0);
    for (auto h : mesh.halfedges(f))
        n += cross(mesh.position(mesh.from_vertex(h)), 
                   mesh.position(mesh.to_vertex(h)));
    return 0.5 * norm(n);
#else
    double a = 0.0;
    Point C = centroid(mesh, f);
    Point Q, R;
    for (auto h : mesh.halfedges(f)) {
        Q = mesh.position(mesh.from_vertex(h));
        R = mesh.position(mesh.to_vertex(h));
        a += pmp::triangle_area(C, Q, R);
    }
    return a;
#endif
}

//-----------------------------------------------------------------------------

Point area_weighted_centroid(const SurfaceMesh &mesh)
{
    Point center(0, 0, 0), c;
    Scalar area(0), a;
    for (auto f : mesh.faces())
    {
        int count = 0;
        c = Point(0, 0, 0);
        for (auto v : mesh.vertices(f))
        {
            c += mesh.position(v);
            count++;
        }
        c /= (Scalar)count;
        a = (Scalar)face_area(mesh, f);
        area += a;
        center += a * c;
    }
    return center /= area;
}

//-----------------------------------------------------------------------------

Eigen::Vector3d gradient_hat_function(const Point& i, const Point& j, const Point& k)
{
    Point base, site, grad;
    Eigen::Vector3d gradient;
    double area;
    area = triangle_area(i, j, k);
    site = i - j;
    base = k - j;
    grad = site - (dot(site, base) / norm(base)) * base / norm(base);
    if (area < eps)
    {
        gradient = Eigen::Vector3d(0, 0, 0);
    }
    else
    {
        grad = norm(base) * grad / norm(grad);
        gradient = Eigen::Vector3d(grad[0], grad[1], grad[2]) / (2.0 * area);
    }

    return gradient;
}

//-----------------------------------------------------------------------------

void setup_gradient_matrix(const SurfaceMesh &mesh, SparseMatrix &G)
{
    SparseMatrix A;
    setup_prolongation_matrix(mesh, A);

    const unsigned int nv = mesh.n_vertices();
    const unsigned int nf = mesh.n_faces();
    Point p, p0, p1;
    Vertex v0, v1;
    int nr_triangles = 0;
    int k = 0;
    auto area_points = mesh.get_face_property<Point>("f:point");
    Eigen::Vector3d gradient_p, gradient_p0, gradient_p1;
    // nonzero elements of G as triplets: (row, column, value)
    std::vector<Triplet> triplets;

    for (Face f : mesh.faces())
    {
        nr_triangles += mesh.valence(f);
        p = area_points[f];
        for (auto h : mesh.halfedges(f))
        {
            v0 = mesh.from_vertex(h);
            v1 = mesh.to_vertex(h);

            p0 = mesh.position(v0);
            p1 = mesh.position(v1);

            gradient_p = gradient_hat_function(p, p0, p1);
            gradient_p0 = gradient_hat_function(p0, p1, p);
            gradient_p1 = gradient_hat_function(p1, p, p0);

            for (int j = 0; j < 3; j++)
            {
                triplets.emplace_back(3 * k + j, nv + f.idx(), gradient_p(j));
                triplets.emplace_back(3 * k + j, v0.idx(), gradient_p0(j));
                triplets.emplace_back(3 * k + j, v1.idx(), gradient_p1(j));
            }
            k++;
        }
    }

    G.resize(3 * nr_triangles, nv + nf);
    G.setFromTriplets(triplets.begin(), triplets.end());
    G = G * A;
}

//-----------------------------------------------------------------------------

void setup_divergence_matrix(const SurfaceMesh &mesh, SparseMatrix &Gt)
{
    SparseMatrix G, M;
    setup_gradient_matrix(mesh, G);
    setup_gradient_mass_matrix(mesh, M);
    Gt = -G.transpose() * M;
}

//-----------------------------------------------------------------------------

void setup_gradient_mass_matrix(const SurfaceMesh &mesh,
                                Eigen::SparseMatrix<double> &M)
{
    auto area_points = mesh.get_face_property<Point>("f:point");
    double area;
    std::vector<Eigen::Triplet<double>> triplets;
    int valence, idx, c = 0;
    for (auto f : mesh.faces())
    {
        valence = mesh.valence(f);
        int i = 0;
        for (auto h : mesh.halfedges(f))
        {
            Point p0 = mesh.position(mesh.from_vertex(h));
            Point p1 = mesh.position(mesh.to_vertex(h));
            area = triangle_area(p0, p1, area_points[f]);
            for (int j = 0; j < 3; j++)
            {
                idx = c + 3 * i + j;
                triplets.emplace_back(idx, idx, area);
            }
            i++;
        }
        c += valence * 3;
    }
    M.resize(c, c);

    M.setFromTriplets(triplets.begin(), triplets.end());
}

//-----------------------------------------------------------------------------

void setup_virtual_vertices(SurfaceMesh &mesh)
{
    auto area_points = mesh.get_face_property<Point>("f:point");
    auto area_weights = mesh.get_face_property<Eigen::VectorXd>("f:weights");

    Eigen::VectorXd w;
    Eigen::MatrixXd poly;

    for (Face f : mesh.faces())
    {
        const int n = mesh.valence(f);
        poly.resize(n, 3);
        int i = 0;
        for (Vertex v : mesh.vertices(f))
        {
            poly.row(i) = (Eigen::Vector3d) mesh.position(v);
            i++;
        }

        compute_virtual_vertex(poly, w);

        area_points[f]  = poly.transpose() * w;
        area_weights[f] = w;
    }
}

//-----------------------------------------------------------------------------

void compute_virtual_vertex(const Eigen::MatrixXd &poly, Eigen::VectorXd &weights)
{
    int val = poly.rows();
    Eigen::MatrixXd J(val, val);
    Eigen::VectorXd b(val);
    weights.resize(val);

    for (int i = 0; i < val; i++)
    {
        Eigen::Vector3d pk = poly.row(i);

        double Bk1_d2 = 0.0;
        double Bk1_d1 = 0.0;

        double Bk2_d0 = 0.0;
        double Bk2_d2 = 0.0;

        double Bk3_d0 = 0.0;
        double Bk3_d1 = 0.0;

        double CBk = 0.0;
        Eigen::Vector3d d = Eigen::MatrixXd::Zero(3, 1);

        for (int j = 0; j < val; j++)
        {
            Eigen::Vector3d pi = poly.row(j);
            Eigen::Vector3d pj = poly.row((j + 1) % val);
            d = pi - pj;

            double Bik1 = d(1) * pk(2) - d(2) * pk(1);
            double Bik2 = d(2) * pk(0) - d(0) * pk(2);
            double Bik3 = d(0) * pk(1) - d(1) * pk(0);

            double Ci1 = d(1) * pi(2) - d(2) * pi(1);
            double Ci2 = d(2) * pi(0) - d(0) * pi(2);
            double Ci3 = d(0) * pi(1) - d(1) * pi(0);

            Bk1_d1 += d(1) * Bik1;
            Bk1_d2 += d(2) * Bik1;

            Bk2_d0 += d(0) * Bik2;
            Bk2_d2 += d(2) * Bik2;

            Bk3_d0 += d(0) * Bik3;
            Bk3_d1 += d(1) * Bik3;

            CBk += Ci1 * Bik1 + Ci2 * Bik2 + Ci3 * Bik3;
        }
        for (int k = 0; k < val; k++)
        {
            Eigen::Vector3d xj = poly.row(k);
            J(i, k) = 0.5 * (xj(2) * Bk1_d1 - xj(1) * Bk1_d2 + xj(0) * Bk2_d2 -
                             xj(2) * Bk2_d0 + xj(1) * Bk3_d0 - xj(0) * Bk3_d1);
        }
        b(i) = 0.5 * CBk;
    }

    Eigen::MatrixXd M(val + 1, val);
    M.block(0, 0, val, val) = 4 * J;
    M.block(val, 0, 1, val).setOnes();

    Eigen::VectorXd b_(val + 1);
    b_.block(0, 0, val, 1) = 4 * b;

    b_(val) = 1.;

    weights = M.completeOrthogonalDecomposition().solve(b_).topRows(val);
}

//=============================================================================
