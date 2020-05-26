//=============================================================================
#pragma once
//=============================================================================

#include <pmp/visualization/MeshViewer.h>
#include "Smoothing.h"
#include "SpectralProcessing.h"

//=============================================================================

class Viewer : public pmp::MeshViewer 
{
public:

    Viewer(const char *title, int width, int height)
            : MeshViewer(title, width, height),
              smooth_(mesh_),
              analyzer_(mesh_),
              compare_sphere(false),
              compare_cube(false) {
        set_draw_mode("Hidden Line");
    }

    virtual bool load_mesh(const char *filename) override;

protected:

    virtual void keyboard(int key, int code, int action, int mod) override;
    virtual void process_imgui() override;
    virtual void update_mesh() override;
    virtual void draw(const std::string &_draw_mode) override;
    virtual void mouse(int button, int action, int mods) override;

    void dualize();
    void close_holes();
    void open_holes();

private:

    Smoothing smooth_;
    SpectralProcessing analyzer_;

    bool compare_sphere;
    bool compare_cube;

    std::vector<Face> holes_;
};

//=============================================================================
