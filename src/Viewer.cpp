//=============================================================================
// Copyright 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "Viewer.h"
#include "Parameterization.h"
#include "Smoothing.h"
#include "PolyDiffGeo.h"
#include "MeanCurvature.h"
#include "GeodesicsInHeat.h"

#include <pmp/algorithms/DifferentialGeometry.h>
#include <pmp/algorithms/SurfaceSubdivision.h>

#include <imgui.h>
#include <random>

//=============================================================================

using namespace pmp;

//=============================================================================

bool Viewer::load_mesh(const char *filename)
{
    bool success = MeshViewer::load_mesh(filename);
    set_draw_mode("Hidden Line");
    return success;
}

//----------------------------------------------------------------------------

void Viewer::process_imgui()
{
    // add standard mesh info stuff
    pmp::MeshViewer::process_imgui();

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Settings", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Clamp cotan", &clamp_cotan_);
    }

    ImGui::Spacing();
    ImGui::Spacing();

    // turn mesh into non-triangles
    if (ImGui::CollapsingHeader("Polygons!", ImGuiTreeNodeFlags_DefaultOpen))
    {
        // Catmull-Clark subdivision
        if (ImGui::Button("Catmull-Clark"))
        {
            SurfaceSubdivision(mesh_).catmull_clark();
            update_mesh();
        }

        // dualize the mesh
        if (ImGui::Button("Dualize mesh"))
        {
            dualize();
        }

        if (ImGui::Button("Kugelize"))
        {
            for (auto v : mesh_.vertices())
                mesh_.position(v) = normalize(mesh_.position(v));
            update_mesh();
        }
    }

    ImGui::Spacing();
    ImGui::Spacing();

    // discrete harmonic parameterization
    if (ImGui::CollapsingHeader("Parametrization",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Discrete Harmonic"))
        {
            Parameterization(mesh_).harmonic_free_boundary();
            mesh_.use_checkerboard_texture();
            set_draw_mode("Texture");

            update_mesh();
        }
    }

    ImGui::Spacing();
    ImGui::Spacing();

    // implicit smoothing
    if (ImGui::CollapsingHeader("Smoothing", ImGuiTreeNodeFlags_DefaultOpen))
    {
        static float timestep = 0.1;
        float lb = 0.001;
        float ub = 1.0;
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("TimeStep", &timestep, lb, ub);
        ImGui::PopItemWidth();

        if (ImGui::Button("Implicit Smoothing"))
        {
            close_holes();
            Scalar dt = timestep;
            smooth_.implicit_smoothing(dt);
            update_mesh();
            BoundingBox bb = mesh_.bounds();
            set_scene((vec3)bb.center(), 0.5 * bb.size());
            open_holes();
        }
    }

    ImGui::Spacing();
    ImGui::Spacing();

    // curvature visualization
    if (ImGui::CollapsingHeader("Curvature", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Mean Curvature"))
        {
            Curvature curv(mesh_);
            curv.compute();
            curv.curvature_to_texture_coordinates();
            mesh_.use_cold_warm_texture();
            update_mesh();
            set_draw_mode("Texture");
        }
    }
    if (ImGui::CollapsingHeader("Geodesics in Heat",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Compute Distances Vertex 0"))
        {
            GeodesicsInHeat heat(mesh_);
            heat.precompute();
            heat.compute_distance_from(Vertex(0));
            heat.distance_to_texture_coordinates();
            mesh_.use_checkerboard_texture();
            update_mesh();
            set_draw_mode("Texture");
        }
    }
}

//----------------------------------------------------------------------------

void Viewer::draw(const std::string &draw_mode)
{
    mesh_.draw(projection_matrix_, modelview_matrix_, draw_mode);
}

//----------------------------------------------------------------------------

void Viewer::dualize()
{
    SurfaceMeshGL dual;

    auto fvertex = mesh_.add_face_property<Vertex>("f:vertex");
    for (auto f : mesh_.faces())
    {
        fvertex[f] = dual.add_vertex(centroid(mesh_, f));
    }

    for (auto v : mesh_.vertices())
    {
        if (!mesh_.is_boundary(v))
        {
            std::vector<Vertex> vertices;
            for (auto f : mesh_.faces(v))
                vertices.push_back(fvertex[f]);
            dual.add_face(vertices);
        }
    }

    mesh_ = dual;
    update_mesh();
}

//----------------------------------------------------------------------------

void Viewer::update_mesh()
{
    // re-compute face and vertex normals
    mesh_.update_opengl_buffers();
}

//----------------------------------------------------------------------------

void Viewer::close_holes()
{
    bool finished = false;
    std::vector<Face> holes;
    while (!finished)
    {
        finished = true;

        // loop through all vertices
        for (auto v : mesh_.vertices())
        {
            // if we find a boundary vertex...
            if (mesh_.is_boundary(v))
            {
                // trace boundary loop
                std::vector<Vertex> vertices;
                vertices.push_back(v);
                for (Halfedge h = mesh_.halfedge(v); mesh_.to_vertex(h) != v;
                     h = mesh_.next_halfedge(h))
                {
                    vertices.push_back(mesh_.to_vertex(h));
                }

                // add boudary loop as polygonal face
                Face f = mesh_.add_face(vertices);
                holes.push_back(f);
                // start over
                finished = false;
                break;
            }
        }
    }
    holes_ = holes;
    update_mesh();
}

//----------------------------------------------------------------------------

void Viewer::open_holes()
{
    for (Face f : holes_)
        mesh_.delete_face(f);
    mesh_.garbage_collection();
    update_mesh();
}

//----------------------------------------------------------------------------

void Viewer::mouse(int button, int action, int mods)
{
    if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_MIDDLE &&
        mods == GLFW_MOD_SHIFT)
    {
        double x, y;
        cursor_pos(x, y);
        Vertex v = pick_vertex(x, y);
        if (mesh_.is_valid(v))
        {
            GeodesicsInHeat heat(mesh_);
            heat.precompute();
            heat.compute_distance_from(v);
            heat.distance_to_texture_coordinates();
            update_mesh();
            mesh_.use_checkerboard_texture();
            set_draw_mode("Texture");
        }
    }
    else
    {
        MeshViewer::mouse(button, action, mods);
    }
}

//=============================================================================
