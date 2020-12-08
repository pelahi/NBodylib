================================================================================================

# NBodylib

================================================================================================

## developed by:

Main Developers
 - Pascal Jahan Elahi (continuously)

Additional contributors:
 - James Willis
 - Rodrigo Canas
 - Rodrigo Tobar

================================================================================================

## Content

This is brief description of the package. This library contains several classes
that are useful in analysing discrete element simulations. It contains a Particle class
for storing information, a KDTree binary tree class for analysing a system of particles
and several other associated math routines. It is principally developed to
provide the classes used by the [VELOCIraptor](https://github.com/pelahi/VELOCIraptor-STF)
halo finder used in analysing astrophysical simulations.

    doc/        contains Doxygen generated latex and html file of code
    src/        contains source code for the libraries
    lib/        library containing libMath, libNBody, libKD, libAnalysis, libCosmology,
                libInitCond
    include/     include files

================================================================================================

## Compiling (see documentation for more information)

NBodylib uses CMake as its build tool. cmake is used to perform system-level checks,
like looking for libraries and setting up the rules for the build, and then generates the
actual build scripts in one of the supported build systems. Among other things, cmake supports
out-of-tree builds (useful to keep more than one build with different settings, and to avoid
cluttering the original source code directories) and several build system, like make and
ninja files.

The simplest way of building is, standing on the root your repository, run cmake to produce
Makefiles and then compile with these steps:

    mkdir build
    cd build
    cmake .. # By default will generate Makefiles
    make all

There are a variety of options that can be invoked and these can be viewed using

    cmake -LH ./
    ccmake ./

(though this only works after having run cmake at least once)

Although documentation is present on the readthedocs site, extra html and latex
documentation can be produced with Doxygen.

    make doc
    doc/html/index.html
    doc/latex/refman.tex

Note that this repo and all variants do not support non-Unix environments. (Mac OS X is fine; Windows is not).

### CI State

| Platform | Branch          | Status |
| -------- | --------------  | ------ |
| UNIX     | ``master``      | [![master](https://travis-ci.com/pelahi/NBodylib.png?branch=master)](https://travis-ci.com/pelahi/NBodylib/branches) |
| UNIX     | ``development`` | [![development](https://travis-ci.com/pelahi/NBodylib.png?branch=development)](https://travis-ci.com/pelahi/NBodylib/branches) |

================================================================================================

## Using Library
