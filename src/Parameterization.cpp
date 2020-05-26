//=============================================================================

#include "Parameterization.h"
#include "PolyDiffGeo.h"
#include <slice.h>
#include <slice_into.h>

//=============================================================================

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

bool Parameterization::harmonic_free_boundary() 
{
    // get properties
    auto tex = mesh_.vertex_property<TexCoord>("v:tex");

    Vertex vh;
    Halfedge hh;

    // Initialize all texture coordinates to the origin.
    for (auto v : mesh_.vertices())
        tex[v] = TexCoord(0.5, 0.5);
    
    // construct area operator from boundary loop
    const int nv = mesh_.n_vertices();
    std::vector<Eigen::Triplet<double>> tripA;
    std::vector<Eigen::Triplet<double>> tripQ;
    std::vector<char> flag(mesh_.n_halfedges(), 0);
    
    SurfaceMesh::HalfedgeIterator hit;
        
    while(1)
    {
        for(hit = mesh_.halfedges_begin(); hit != mesh_.halfedges_end(); ++hit)
        {
            if(!flag[(*hit).idx()] && mesh_.is_boundary(*hit))
            {
                break;
            }
        }
        
        if(hit ==  mesh_.halfedges_end()) break;

        hh = *hit;
        auto start = hh;
        
        do {
            const int i = mesh_.to_vertex(hh).idx();
            const int j = mesh_.to_vertex(mesh_.next_halfedge(hh)).idx();
                    
            tripA.emplace_back(i, nv + j, -0.5);
            tripA.emplace_back(j, nv + i, 0.5);
            
            tripA.emplace_back(nv + j, i, -0.5);
            tripA.emplace_back(nv + i, j, 0.5);
            
            tripQ.emplace_back(i, i, 1.);
            tripQ.emplace_back(nv + i, nv + i, 1.);
            
            flag[hh.idx()] = 1;
            
            hh = mesh_.next_halfedge(hh);
        } while (hh != start);
    }
    
    // no boundary found ?
    if (tripQ.empty()) {
        std::cerr << "Mesh has no boundary." << std::endl;
        return false;
    }
    
    Eigen::SparseMatrix<double> A(2 * nv, 2 * nv);
    A.setFromTriplets(tripA.begin(), tripA.end());
    
    SparseMatrix L, L2(2 * nv, 2 * nv);
    setup_stiffness_matrix(mesh_, L);
    
    // build L2 = [L 0; 0 L]
    Eigen::VectorXi range(nv);
    
    for(int i = 0; i < nv; ++i) range(i) = i;
    igl::slice_into(L, range, range, L2);

    for(int i = 0; i < nv; ++i) range(i) = nv + i;
    igl::slice_into(L, range, range, L2);

    // build full matrix
    Eigen::SparseMatrix<double> B = -(L2 + A);
    
    // select partial matrices
    int nb = 0;
    for (auto v : mesh_.vertices()) {
        if (mesh_.is_boundary(v)) {
            nb++;
        }
    }
    
    Eigen::VectorXi boundary(2 * nb), inner(2 * (nv - nb));
    int cnt0 = 0;
    int cnt1 = 0;
    
    for (auto v : mesh_.vertices()) {
        if (mesh_.is_boundary(v)) {
            boundary(cnt0++) = v.idx();
            boundary(cnt0++) = nv + v.idx();
        
        } else
        {
            inner(cnt1++) = v.idx();
            inner(cnt1++) = nv + v.idx();
        }
    }
    
    Eigen::SparseMatrix<double> BII, BBB, BIB;
    igl::slice(B, inner, inner, BII);
    igl::slice(B, inner, boundary, BIB);
    igl::slice(B, boundary, boundary, BBB);
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol;
    chol.analyzePattern(BII);
    chol.factorize(BII);
    
    Eigen::MatrixXd B2 = BBB.toDense() - BIB.transpose() * chol.solve(BIB);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigs(B2);
   
    // find first eigenvalue above threshold (1e-8).
    // Eigen::SelfAdjointEigenSolver seems to produce sorted eigenvalues. Not sure if this is guaranteed though! Eigen::EigenSolver does not follow this convention!
    int evid = 0;
    
    for(;evid < B2.cols(); ++evid)
        if(std::abs(eigs.eigenvalues()(evid)) > 1e-8) break;
    
    Eigen::VectorXd xb = eigs.eigenvectors().col(evid);
    Eigen::VectorXd xi = chol.solve(- BIB * xb );
    Eigen::VectorXd x(2 * nv);
    
    for(int i = 0; i < inner.size(); ++i)
    {
        x(inner(i)) = xi(i);
    }
    
    for(int i = 0; i < boundary.size(); ++i)
    {
        x(boundary(i)) = xb(i);
    }
    
    x.array() -= x.minCoeff();
    x.array() /= x.maxCoeff();
        
    for (auto v : mesh_.vertices())
    {
        tex[v] = TexCoord(x(v.idx()), x(nv + v.idx()));
    }
    
    return true;
}

//-----------------------------------------------------------------------------

bool Parameterization::harmonic() 
{
    // map boundary to circle
    if (!setup_boundary_constraints()) {
        std::cerr << "Could not perform setup of boundary constraints.\n";
        return false;
    }

    // get property for 2D texture coordinates
    auto tex = mesh_.vertex_property<TexCoord>("v:tex");
    // Compute the chosen implicit points and its convex combination weights for each face

    int nb = 0;
    for (auto v : mesh_.vertices()) {
        if (mesh_.is_boundary(v)) {
            nb++;
        }
    }
    const unsigned int nv = mesh_.n_vertices();

    unsigned k = 0;
    unsigned ins = 0;
    unsigned out = 0;

    Eigen::MatrixXd B(nv, 2);
    std::vector<Vertex> vertices;
    Vertex v;
    vertices.reserve(mesh_.n_vertices());
    Eigen::VectorXi in(nv - nb), b(nb), x(2);
    x << 0, 1;
    for (auto v : mesh_.vertices()) {
        vertices.push_back(v);
        if (!mesh_.is_boundary(v)) {
            in(ins) = k;
            B(k, 0) = 0.0;
            B(k, 1) = 0.0;
            ins++;
        } else {
            b(out) = k;
            B(k, 0) = tex[v][0];
            B(k, 1) = tex[v][1];
            out++;
        }
        k++;
    }

    SparseMatrix L;
    Eigen::SparseMatrix<double> L_in_in, L_in_b;
    Eigen::MatrixXd b_in, b_out, X;

    setup_stiffness_matrix(mesh_, L);

    igl::slice(L, in, in, L_in_in);
    igl::slice(L, in, b, L_in_b);
    igl::slice(B, in, x, b_in);
    igl::slice(B, b, x, b_out);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L_in_in);
    X = solver.solve(b_in - L_in_b * b_out);

    if (solver.info() != Eigen::Success) {
        std::cerr << "harmonic(): Could not solve linear system\n";
    } else {
        // copy solution
        k = 0;
        for (unsigned int i = 0; i < nv; ++i) {
            v = vertices[i];
            if (!mesh_.is_boundary(v)) {
                tex[v][0] = X(k, 0);
                tex[v][1] = X(k, 1);
                k++;
            }
        }
    }
    return true;
}

//-----------------------------------------------------------------------------

bool Parameterization::setup_boundary_constraints() 
{
    // get properties
    auto points = mesh_.vertex_property<Point>("v:point");
    auto tex = mesh_.vertex_property<TexCoord>("v:tex");

    SurfaceMesh::VertexIterator vit, vend = mesh_.vertices_end();
    Vertex vh;
    Halfedge hh;
    std::vector<Vertex> loop;

    // Initialize all texture coordinates to the origin.
    for (auto v : mesh_.vertices())
        tex[v] = TexCoord(0.5, 0.5);

    // find 1st boundary vertex
    for (vit = mesh_.vertices_begin(); vit != vend; ++vit)
        if (mesh_.is_boundary(*vit))
            break;

    // no boundary found ?
    if (vit == vend) {
        std::cerr << "Mesh has no boundary." << std::endl;
        return false;
    }

    // collect boundary loop
    vh = *vit;
    hh = mesh_.halfedge(vh);
    do {
        loop.push_back(mesh_.to_vertex(hh));
        hh = mesh_.next_halfedge(hh);
    } while (hh != mesh_.halfedge(vh));

    // map boundary loop to unit circle in texture domain
    unsigned int i, n = loop.size();
    Scalar angle, l, length;
    TexCoord t;

    // compute length of boundary loop
    for (i = 0, length = 0.0; i < n; ++i)
        length += distance(points[loop[i]], points[loop[(i + 1) % n]]);

    // map length intervalls to unit circle intervals
    for (i = 0, l = 0.0; i < n;) {
        // go from 2pi to 0 to preserve orientation
        angle = 2.0 * M_PI * (1.0 - l / length);

        t[0] = 0.5 + 0.5 * cosf(angle);
        t[1] = 0.5 + 0.5 * sinf(angle);

        tex[loop[i]] = t;

        ++i;
        if (i < n) {
            l += distance(points[loop[i]], points[loop[(i + 1) % n]]);
        }
    }

    return true;
}

//=============================================================================
