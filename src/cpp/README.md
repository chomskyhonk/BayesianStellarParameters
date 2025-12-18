
# C++ Components for Bayesian Stellar Parameters

This directory contains C++ source files used for various aspects of Bayesian stellar parameter analysis. The main functionality includes automating PDF calculations, performing likelihood and isogrid constraints, and calculating astrometric values.

## Purpose

The C++ code handles:
- Data processing: Reading and preprocessing stellar parameter data
- Statistical analysis: Implementing Bayesian inference algorithms
- Performance-critical computations: Likelihood calculations

## Source Files

### Key Source Files and Their Purpose

1. `automatepdfs.cpp`:
   - Automates the process of calculating PDFs for all stars.
   - Calls the likelihood (`gaialikelihood.cpp`) and isogrid constraint (`gaiaisogrid.cpp`) modules in tandem for each star.

2. `gaialikelihood.cpp`:
   - Computes the likelihood distributions for individual stars based on input data.

3. `gaiaisogrid.cpp`:
   - Handles isochrone grid constraints for individual stars.

4. `astrometry.cpp`:
   - Performs astrometric calculations for stellar parameter analysis across all stars with Bayesian parameters.

---

## Compilation and Build Instructions

To compile the C++ source files and create the desired executables, follow these steps:

### Prerequisites

1. Compiler & Tools:
   - `clang++` (C++17 or higher)
   - `make`

2. Libraries & Dependencies:
   - Python:
     - Include path: `/Users/billyharrison/anaconda3/include/python3.11`
     - Library path: `/Users/billyharrison/anaconda3/lib`
   - `libomp`: For OpenMP parallelization support.
   - Eigen3: Linear algebra library (`/usr/local/opt/eigen`).

3. Ensure all paths to external libraries are correctly defined in the `Makefile` for successful compilation.

---

### Build Commands

Use the provided `Makefile` to build and run individual programs.

1. Automate PDF Calculations Across All Stars:
   Compile the default target to automate PDF calculations:
   ```bash
   make TARGET=automatepdfs
   ./automatepdfs path/to/input_data.dat
   ```

2. Individual Operations:
   Update the `TARGET` in the `Makefile` as follows to compile and run individual components:
   - **Likelihood Calculations**:
     ```bash
     make TARGET=gaialikelihood
     ./gaialikelihood path/to/input_likelihood_data.dat
     ```

   - Isochrone Grid Constraints:
     ```bash
     make TARGET=gaiaisogrid
     ./gaiaisogrid path/to/input_isogrid_data.dat
     ```

   - Astrometry Calculations:
     ```bash
     make TARGET=astrometry
     ./astrometry path/to/input_astrometry_data.dat
     ```

3. Clean Build Artifacts:
   Use the `clean` target to remove compiled files:
   ```bash
   make clean
   ```

---

## Dependencies

The following dependencies must be installed to build and run the C++ components in this directory:

### Compiler and Tools
- C++ Compiler:
  - Required: A C++17-compatible compiler (GCC 8.1+, Clang 7+, or MSVC with `/std:c++17`).
- Build Tools:
  - `make` (for build automation).

### Libraries
1. OpenMP:
   - For parallelization using `<omp.h>`.
   - Installation:
     - macOS: `brew install libomp`
     - Linux: Ensure GCC is installed (`sudo apt install gcc`).
     - Windows: Check MSVC compiler options for OpenMP support.
2. Eigen:
   - For matrix operations using `<Eigen/Dense>`.
   - Installation:
     - macOS: `brew install eigen`
     - Linux: `sudo apt install libeigen3-dev`
     - Windows: Download from [Eigen's website](https://eigen.tuxfamily.org/).
3. pybind11:
   - For embedding Python and interacting with Python objects (`<pybind11/numpy.h>`, `<pybind11/stl.h>`).
   - Installation:
     ```bash
     pip install pybind11
     ```
4. matplotlib-cpp:
   - For creating plots from C++ code (`"matplotlibcpp.h"`).
   - Requirements:
     - Python: Install Matplotlib and NumPy:
       ```bash
       pip install matplotlib numpy
       ```
     - Install `pybind11` as described above.
     - Clone the `matplotlib-cpp` repository into your `src` folder.

### Project-Specific Dependencies
- Gaia-related modules (`gaiaisogrid.h`, `gaialikelihood.h`, `isochronetypes.h`) must be present in the source tree or properly linked using the `Makefile`.

---

### Verifying Dependencies
- After installing the required dependencies, build the project using:
  ```bash
  make all
  ```
- If any errors occur, ensure the paths to external libraries (e.g., Eigen, pybind11) are specified correctly in the `Makefile`.

## Example Input and Output

### Input Files
Each C++ program expects specific input data files, such as:
- Star parameter data (`.csv` files)

### Output Files
The executables generate the following output files:
- `automatepdfs`: Combined PDFs for all stars.
- `gaialikelihood`: Likelihood distributions.
- `gaiaisogrid`: Isochrone grid constraints.
- `astrometry`: Astrometry calculations.

Ensure you place input data files in the correct directory and pass the file paths as arguments when running the executables.
Replace all input and output filepaths for data and save location as my current directories are unchanged across .cpp files.

---

## Troubleshooting

If you encounter build or run-time errors, try the following:
- Verify that paths to the required libraries are correctly specified in the `Makefile`.
- Enable debugging flags during compilation:
  ```bash
  make CXXFLAGS="-g -fsanitize=address"
  ```
- Ensure your input data files are correctly formatted.

---
