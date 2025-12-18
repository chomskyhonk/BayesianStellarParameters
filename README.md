# Probabilistic Parameter Determination of Stars and Galaxies

### MSc Research Project – University College London (UCL)

This project explores the use of Bayesian inference to estimate the ages and metallicities of stars in the Solar neighbourhood, leveraging large-scale stellar datasets from Gaia. Accurate age determination is one of the most challenging problems in Galactic astrophysics due to observational uncertainties and degeneracies between stellar parameters and within stellar models like isochrones.  
By applying a probabilistic approach, this project demonstrates how Bayesian modelling and isochrone fitting can be used to generate strong posterior distributions for fundamental stellar parameters.

---

## Project Overview

Traditional stellar parameter estimation methods suffer from large uncertainties and model degeneracies. This project develops a Bayesian inference pipeline to overcome these challenges by:

- Integrating prior information from stellar evolution models and previous literature
- Constraining parameter inference with stellar model probability distribution functions (PDFs)
- Performing probabilistic isochrone fitting to observed Gaia data in the form of the likelihood function  
- Quantifying uncertainty through full posterior PDFs
- Testing the resulting age–velocity and metallicity–dispersion relations against the literature  

The core goal is to validate Bayesian methods as a reliable tool for stellar parameter estimation, and to understand how biases and selection effects influence inferences about the structure and evolution of the Galactic disk, particualrly in the Solar neighbourhood.

---

## Methodology

1. Data Acquisition  
   - Parallax-limited sample of stars from the Gaia database  
   - Queried using ADQL (SQL-like language for astronomical databases)  
   - Preprocessed and cleaned using Python (NumPy, Pandas)

2. Bayesian Modelling  
   - Constructed likelihoods for observed stellar parameters  
   - Combined with priors and constraining models based on stellar evolution and isochrones  
   - Implemented inference via combining these using Bayesian formalism in C++

3. Posterior Analysis 
   - Extracted PDFs for stellar age and metallicity using marginalisation  
   - Evaluated global relations such as:
     - Age–velocity dispersion relation (β = 0.29 ± 0.01)
     - Metallicity–dispersion trends
     - Age-metallicity relation
   - Assessed model bias, overprediction of old/metal-rich stars, and uncertainty propagation

4. Visualisation & Validation  
   - Produced distribution maps and statistical summaries with Matplotlib
   - Compared results with literature to validate physical consistency

---

## Key Results

- Derived an age–velocity dispersion relation (AVR) consistent with theoretical and empirical expectations  
- Identified negative metallicity–dispersion trends, revealing higher kinematic heating at intermediate metallicities  
- Validated Bayesian isochrone fitting as a credible framework for probabilistic stellar parameter estimation  
- Highlighted limitations such as overestimation of old, metal-rich stars, suggesting future refinement of priors and likelihood constraints  

---

## Tools & Technologies

| Category | Tools / Libraries |
|-----------|-------------------|
| Languages | C++, Python |
| Data Handling | Pandas, NumPy |
| Statistical Modelling | SciPy, custom C++ routines |
| Database Querying | ADQL (Astronomical SQL), API |
| Visualisation | Matplotlib, Seaborn |
| Workflow Automation | API-based data collection of isochrones and automated processing pipelines |

---

## Lessons Learned

- Developed strong expertise in probabilistic inference,*data pipeline design, and statistical model validation 
- Gained experience integrating multi-source datasets and ensuring data reproducibility  
- Improved understanding of how biases and selection effects influence model inference  
- Strengthened ability to communicate complex analytical results 

---

## Abstract

> Decoding the age and metallicity distributions of stars in the Solar neighbourhood is central to tracing the kinematic evolution of the Galactic disk, yet direct age determinations are notoriously uncertain. Traditional techniques are limited by large observational uncertainties and isochrone degeneracies. This project employs a **Bayesian framework** to fit isochrones to *Gaia* stars, integrating prior information and probabilistic error treatment. The resulting posterior probability distributions enable robust statistical inference of stellar ages and compositions. The derived age–velocity dispersion relation (β = 0.29 ± 0.01) aligns with literature, with higher vertical heating observed compared to in-plane motion. The model reproduces key Solar neighbourhood trends in disk heating and chemical evolution, while highlighting overprediction of old, metal-rich stars as a limitation for future improvement.  

---


## Repository Structure

```
BayesianStellarParameters/
├── docs/                  # Thesis materials
│   └── README.md         # Documentation guide
├── src/
│   ├── cpp/              # C++ code for data handling and statistical analysis
│   │   └── README.md    # C++ code documentation
│   └── python/           # Python code for visualization
│       └── README.md    # Python code documentation
├── data/                 # Data files (excluded from git)
│   └── README.md        # Data documentation
├── .gitignore           # Git ignore rules
└── README.md            # This file
```

## Getting Started

## Prerequisites

Ensure your development environment meets the following prerequisites:

### Compiler & Build Tools
- `clang++` (with support for C++17 or higher)
- `make` (for build automation)

### Python Requirements
- Python 3.11
- Required Python Libraries:
  - `matplotlib-cpp`
  - `pybind11`
  - `numpy`

### Additional Dependencies
- Eigen3: `/usr/local/opt/eigen` - did not end up using
- Anaconda (for Python library and runtime management):
  - Include path: `/Users/billyharrison/anaconda3/include/python3.11`
  - Library path: `/Users/billyharrison/anaconda3/lib`
- Project-specific Paths for External Libraries:
  - Automate PDFs (`automatepdfs`): `/Users/billyharrison/projects/vsc/iniplots/.../isochrone_datasets/automatepdfs/include`
  - Gaia Likelihood (`gaialikelihood`): `/Users/billyharrison/projects/vsc/iniplots/.../isochrone_datasets/gaialikelihood/include`
  - Gaia Isochrone Grid (`gaiaisogrid`): `/Users/billyharrison/projects/vsc/iniplots/.../isochrone_datasets/gaiaisogrid/include`

---

## Installation

Follow the instructions below to set up and build the project.

1. Clone the Repository:
   ```bash
   git clone https://github.com/chomskyhonk/BayesianStellarParameters.git
   cd BayesianStellarParameters
   ```

2. Customize `Makefile` (Optional):
   Depending on the specific task you want to perform, update the `TARGET` and `SRCS` variables in the `Makefile`:
   - For Astrometry Calculations:
     ```make
     TARGET = astrometry
     SRCS = astrometry.cpp
     ```
   - For Likelihood Calculations Only:
     ```make
     TARGET = gaialikelihood
     SRCS = gaialikelihood.cpp
     ```
   - For Isochrone Grid Constraints Only:
     ```make
     TARGET = gaiaisogrid
     SRCS = gaiaisogrid.cpp
     ```
   - For Full Automation Across Desired Stars (default):
     ```make
     TARGET = automatepdfs
     SRCS = automatepdfs.cpp
     ```

3. Build the Project:
   Run the `make` command to compile the project:
   ```bash
   make
   ```
   The specified `TARGET` will be created as an executable in the project root directory.

4. Clean Build Files (Optional):
   To remove all compiled files:
   ```bash
   make clean
   ```

---

## Data Processing

### Automated PDF Calculations for All Stars
1. Ensure your data is properly organized and paths are set in `automatepdfs.cpp`.
2. Build the project with `TARGET = automatepdfs` in the `Makefile`.
3. Run the automation program:
   ```bash
   make run
   ```
   This will:
   - Perform likelihood and isogrid constraints for each star, looping through the dataset, and calculate the posterior.
   - Save output PDFs and likelihood constraints in specified directories.

### Individual Operations
You can also run modular calculations for specific needs w.r.t dispersion, star kinematics:
- Astrometry Calculations:
  1. Update the `Makefile`: `TARGET = astrometry`, `SRCS = astrometry.cpp`.
  2. Compile and run:
     ```bash
     make
     ./astrometry path/to/input_file.dat
     ```

- Likelihood Calculations:
  1. Update the `Makefile`: `TARGET = gaialikelihood`, `SRCS = gaialikelihood.cpp`.
  2. Compile and run:
     ```bash
     make
     ./gaialikelihood path/to/input_file.dat
     ```

- Isochrone Grid Constraints:
  1. Update the `Makefile`: `TARGET = gaiaisogrid`, `SRCS = gaiaisogrid.cpp`.
  2. Compile and run:
     ```bash
     make
     ./gaiaisogrid path/to/input_file.dat
     ```

---

## Example Usage

Here’s how to run the programs for different use cases:

1. Automate PDFs for All Stars:
   ```bash
   make run
   ```

2. Individual Calculations:
   - For `astrometry`:
     ```bash
     ./astrometry path/to/star_data.dat
     ```
   - For `likelihood`:
     ```bash
     ./gaialikelihood path/to/star_data.dat
     ```
   - For `isogrid`:
     ```bash
     ./gaiaisogrid path/to/star_data.dat
     ```
   - ^^here add your path to the desired data/csv file taken from the data folder, or queried yourself
---

## Troubleshooting

1. Build Errors:
   - Verify that all paths in the `Makefile` are configured correctly.
   - Ensure required libraries are installed and accessible.

2. Debugging:
   - Rebuild with debugging flags enabled:
     ```bash
     make CXXFLAGS="-g -fsanitize=address"
     ```

3. Run-time Errors:
   - Ensure input data files match the required format.
   - Use `--help` flags (if supported) to explore further options.

---


## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{bayesian_stellar_parameters,
  author = {Billy Harrison},
  title = {Bayesian Stellar Parameters - Framework for Age and Metallicity PDFs},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/chomskyhonk/BayesianStellarParameters}
}
```

## Contact

For questions or collaboration inquiries:
- GitHub Issues: [Project Issues](https://github.com/chomskyhonk/BayesianStellarParameters/issues)
- (LinkedIn: https://www.linkedin.com/in/billy-harrison-74ab66251)


