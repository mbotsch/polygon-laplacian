# Polygon Laplacian Made Simple

This repository contains the implementation of our paper [Bunge, Herholz, Kazhdan, Botsch, *Polygon Laplacian Made Simple*, Eurographics 2020](https://diglib.eg.org/handle/10.1111/cgf13931).

The code is based on the [pmp-library](http://www.pmp-library.org/) and uses it as a git-submodule. You therefore have to checkout the repository **recursively**:

    git clone --recursive https://github.com/mbotsch/polygon-laplacian.git

Configure and build:

    cd polygon-laplacian && mkdir build && cd build && cmake .. && make

This will automatically build the pmp-library, its dependecies, and the polygon-laplacian demo and tests. You can start the interactive app with a polygon mesh:

    ./polylaplace ../data/boar.obj

Have fun!
