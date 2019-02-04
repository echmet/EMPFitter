ElmigParamsFitter
===

Introduction
---
ElmigParamsFitter library implements a simple algorithm based on nonlinear regression whose purpose is to calculate electrophoretic mobilities and pKa constants of a chemical compound. Input data for the algorithm is a set of effective electrophoretic mobilities of the compound in a series of background buffers of various pH. The algorithm is built on top of a model implemented in [ECHMETCoreLibs](https://github.com/echmet/ECHMETCoreLibs) and can account for Debye-Hückel and Onsager-Fuoss ionic interactions effects.

Building
---
ElmigParamsFitter depends of the following tools and libraries:

- C++14-aware compiler
- [Eigen matrix library](https://eigen.tuxfamily.org/)
- [CMake build system](http://cmake.org)
- [ECHMETCoreLibs](https://github.com/echmet/ECHMETCoreLibs)
- [ECHMETRegressor](https://github.com/echmet/ECHMETRegressor)

### Generate makefiles with CMake
- Linux
1. `cd` into the source directory
2. Run `mkdir build` and `cd build`
3. Run `cmake .. -DCMAKE_BUILD_TYPE=Release -DEIGEN_INCLUDE_DIR=<path_to_eigen_library> -DECHMET_CORE_LIBS_DIR=<path_to_ECHMETCoreLibs_installation> -DECHMET_REGRESSOR_DIR=<path_to_ECHMETRegressor_headers>`
4. Run `make` and `make install`
- Windows
1. You may use CMake GUI to set up the required configuration options as it is described in the Linux section of this README.
2. Generate project files. MSVC 14 is currently the preferred compiler on Windows.

The project is known to build correctly with GCC and Clang on Linux and MSVC 14 on Windows. While other compilers may by used as well it has not been tested by the author.

Licensing
---
The ElmigParamsFitter project is distributed under the terms of **The GNU General Public License v3** (GNU GPLv3). See the enclosed `LICENSE` file for details.

As permitted by section 7. *Additional Terms* of The GNU GPLv3 license, the authors require that any derivative work based on ElmigParamsFitter clearly refers to the origin of the software and its authors. Such reference must include the address of this source code repository (https://github.com/echmet/ElmigParamsFitter) and names of all authors and their affiliation stated in section [Authors](#Authors) of this README file.

<a name="Authors"></a>
Authors
---
Michal Malý    

Group of Electromigration and Chromatographic Methods (http://echmet.natur.cuni.cz)

Department of Physical and Macromolecular Chemistry  
Faculty of Science, Charles University, Czech Republic
