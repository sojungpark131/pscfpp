
# PSCF - Polymer Self-Consistent Field theory (C++/CUDA version)

PSCF is a package of software for solving the Edwards-Helfand self-consistent 
field theory for polymer liquids. The version in the directory containing
this file is written primarily in C++, with GPU accelerated code in CUDA.
This C++/CUDA version of PSCF is still under development, but is intended 
to eventually supersede the existing PSCF Fortran program. The older 
Fortran program is maintained in a separate github repostory dmorse/pscf.

## Overview

Differences between this C++/CUDA version of PSCF from the older Fortran 
version and expected advantages of the new code include:

   - PSCF (C++/CUDA) is an extensible package of several different programs 
     designed for use with different geometries and boundary conditions, 
     different algorithms or different hardware, designed around a 
     common software framework. 

   - PSCF (C++/CUDA) allows simulations of mixtures containing arbitrary 
     acyclic branched copolymers, in addition to the linear block copolymers 
     and homopolymers allowed by Fortran PSCF code.

   - An object oriented design allows creation of independent objects to 
     represent different phases, which will facilitate analysis of phase
     coexistence in mixtures.

   - Adoption of C/C++ as a base language has simplified implementation 
     of variant that that uses graphical process units (GPUs).

The C++/CUDA version of PSCF is free, open source software. It is distributed 
under the terms of the GNU General Public License as (GPL) published by the 
Free Software Foundation, either version 3 of the License or (at your option) 
any later version.  PSCF is distributed without any warranty, without even 
the implied warranty of merchantability or fitness for a particular purpose. 
See the LICENSE file or the 
<a href=http://www.gnu.org/licenses/> gnu web page </a> for details.

## Programs

Currently, the package contains the following SCFT solvers:

   - A one-dimensional finite difference solver for problems that involve 
     variation in a single spatial coordinate in Cartesian, cylindrical 
     or spherical coordinates.

   - A CPU-based pseudo-spectral solver periodic microstructures that are
     periodic in 1, 2, or 3 coordinates.

   - A GPU accelerated pseudo-spectral solver for period microstructures. 

The one-dimensional finite different solver is useful for treating problems
involving flat or curved interfaces, and cylindrical or spherical micelles.
The executable for the this program is named pscf_fd.

The CPU-based pseudo-spectral solver for periodic microstructures is
similar in most respects to the existing PSCF Fortran program, and 
provides similar level of performance. Like the Fortran program, it 
allows the user to search for a solution with any specified crystal
system and space group symmetry, and provides efficient algorithms
to relax the unit cell parameters so as to minimize the free energy.
The new code can read and write the same file formats for representing 
a field in terms of symmetry-adapted basis functions as those used by 
the PSCF Fortran program.  Currently, the most important features of
the Fortran code that have not yet been reimplemented in the new codes
are the "sweep" continuation feature and specialized code to represent
point-particle solvents. Separate executables are used to solve 1, 2 
and 3 dimensionally periodic structures, which are named pscf_pc1d, 
pscf_pc2d and pscf_pc3d, respectively. Here, "pc" stands for 
"periodic CPU".

The GPU-accelerated pseudo-spectral solver for periodic structures 
is based on algorithms similar to those used in the CPU pseudo-spectral 
solver, but is somewhat less mature. Like the corresponding CPU code,
the GPU-accelerated code allows the use of any unit cell type, 
including no-orthogonal unit cells, and automatic relaxation of unit 
cell parameters. The most important difference in features is that the 
GPU-accelerated code does yet allow the user to use symmetry-adapted
basis functions to constrain the space group symmetry of the solution.
The GPU accelerated programs for solving 1, 2 and 3 dimensionally 
periodic structures are named pscf_pg1d, pscf_pg2d and pscf_pg3d, 
respectively, where "pg" stands for "periodic GPU".

## Getting the source code

The PSCF C++/CUDA source code is maintained in the github repository

   <https://github.com/dmorse/pscfpp>.

It may be obtained by using a git version control system client to
clone the repository. To do so, enter the command:
``` 
git clone --recursive https://github.com/dmorse/pscfpp.git
```
The use of the --recursive option to the git clone command:
This is necessary to clone some git submodules that are maintained
in separate repositories. This command will create a directory 
called pscfpp/ that contains all of the source code and associated
documentation.

## Documentation

PSCF is distributed with source files for an html web manual.
After cloning the source code, you can use the doxygen documentation
generator to generate a local copy of this documentation. To do this,
doxygen must be installed on your computer, and the directory 
containing the doxygen executable must be in your command search
PATH. To generate documentation:

   - Change directory (cd) to the pscfpp/ root directory

   - Enter "make html"

This should create many html files in the pscfpp/doc/html directory.
To begin reading the documentation, point a browser at the file
pscfpp/doc/html/index.html, which is the main page of the manual.

## Dependencies

The PSCF source code is written in a combination of C++ and (for
the GPU accelerated program) CUDA, and must be compiled from source.
The package was developed on linux and and Mac OS X operating systems 
using standard unix utilities, and is designed to run on these 
systems. To compile linux-like software on a Mac, you must first 
install the XCode Mac development environment and the unix command 
line tools.  

The CPU-based programs within the PSCF package depend on the 
following external libraries:

  - Gnu scientific library (GSL)

  - FFTW fast Fourier transform library

The one-dimensional finite difference program pscf_fd requires 
only GSL, and not FFTW. The CPU-based programs for spatially
periodic structures require both GSL and FFTW libraries.

The GPU-accelerated programs can only run on a computer with an
appropriate NVIDIA graphics card. To compile these programs, the
system must also have an NVIDIA CUDA development environment 
that provides the CUFFT fast Fourier transform library. 

## Compiling

Complete directions for compiling and installing PSCF are
provided in section 2 of the html documentation. Short instructions
for compiling, after installing all of the required dependencies,
are given below:

   - Add the pscfpp/bin directory to your linux command search PATH
     environment variable.
   
   - Add the pscfpp/scripts/python directory to your PYTHONPATH
     environment variable.
   
   - cd to the pscfpp/ root directory
   
   - Enter "./setup" from this root directory to run a setup script
     (you only need to do this once, before compiling the first time).

   - Change directory (cd) to the subdirectory pscfpp/bld/.
   
   - To compile and install all CPU-based programs in the package 
     (excluding GPU-accelerated programs), enter "make all-cpu"
   
   - To compile the GPU-accelerated programs on a machine with an
     NVIDIA GPU, a CUDA compiler and the CUFFT library, enter
     "make pspg". 

The setup script installs a file pscfpp/bld/config.mk that contains
makefile variables that define compiler executable names, compiler options 
and paths to head and library files for external dependencies.  If the 
default options are not adequate, edit this file as needed.

## Command line syntax (invoking a program)

PSCF is a package containing several different SCFT programs designed 
for different geometries, different algorithms or different hardware. 
Executable names (given above) are:

   - pscf_fd : 1D finite-difference program 

   - pscf_pcNd : CPU based programs for N=1,2, or 3 dimensional periodic
     structures

   - pscf_pgNd : GPU based programs for N=1,2, or 3 dimensional periodic
     structures

In the names pscf_pcdN and pscf_pgdN, N denotes a dimension dimension
of space that can be N=1, 2, 3. The CPU program for three-dimensionally
periodic microstructures is thus pscf_pc3d.

Each of these programs reads a parameter file and a command file. The 
parameter file is fixed-format file that contains parameters required 
to initialize the program. The command file is a more flexible script 
containing a sequence of commands that are read and executed sequentially 
to specify a sequence of computational steps.  The command file for 
a standard SCFT calculation also specifies the name of a file that 
contain an initial guess for monomer chemical potential fields and 
names of files to which final chemical potential and monomer 
concentration fields should be written.

The command line syntax for invoking any pscfp++ program is:
```
program -p param -c command
```
where "program" denotes the name of the program, "param" denotes the path 
to a parameter file, and "command" denotes the path to a parameter file. 
For example, one might enter
```
pscf_pc3d -p param -c command
```
to run the pscf_pc3d CPU program for three dimensional periodic structures. 
This form of the command would write log output to the screen.  Output 
produced during a computation may also be redirected to a log file by 
using the unix ">" standard output redirection operator. In addition, 
and -e command line option may also be used to cause the program to echo 
the parameter file to standard out as this file is read. With echoing and 
redirection, the command syntax would be 
```
program -e -p param -c command > log
```
where "log" denotes the name of a log file to which output will be
written during the computation.

## Examples

Directory pscfpp/examples contains a set of examples of simple 
calculations, each of which contains a sample parameter and command file. 
Top level subdirectories of pscfpp/examples contain examples for different 
PSCF programs.

Subdirectory examples/fd1d subdirectory contains examples for the 1D 
finite-difference program pscf_fd. Top level subdirectories of 
examples/fd1d contain examples for planar, cylindrical and spherical 
geometries, as indicated by the subdirectory names. One or more example 
is given for each type of geometry.

Subdirectory examples/pspc contains examples contains examples the CPU 
based pscf_pcNd programs for N=1, 2 and 3. Top level subdirectories
contain solutions for a particular type of physical system, e.g.,
examples/pspc/diblock contains examples for a diblock copolymer melt. 
Subdirectories of examples/pspc/diblock contain examples for lamellar 
(N=1), hexagonal (N=2) and BCC (N=3) structures, each of which has a 
different number of spatially periodic dimensions.
