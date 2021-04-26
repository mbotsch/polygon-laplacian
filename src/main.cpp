//=============================================================================
// Copyright 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "Viewer.h"

//=============================================================================

int main(int argc, char **argv)
{
    Viewer window("Polygon Laplacian", 800, 600);

#ifndef __EMSCRIPTEN__
    if (argc == 2)
        window.load_mesh(argv[1]);
#else
    window.load_mesh(argc == 2 ? argv[1] : "input.off");
#endif

    return window.run();
}

//=============================================================================
