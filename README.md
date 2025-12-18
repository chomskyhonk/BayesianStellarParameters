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
| Visualisation | Matplotlib |
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
├── docs/                  # Documentation and thesis materials
│   └── README.md         # Documentation guide
├── src/
│   ├── cpp/              # C++ code for data handling and statistical analysis
│   │   └── README.md    # C++ code documentation
│   └── python/           # Python code for visualization
│       └── README.md    # Python code documentation
├── examples/             # Example scripts and notebooks
│   └── README.md        # Examples guide
├── data/                 # Data files (excluded from git)
│   └── README.md        # Data documentation
├── .gitignore           # Git ignore rules
└── README.md            # This file
```

## Getting Started

### Prerequisites

C++ Requirements:
- C++11 or later
- Standard library
- (Add other dependencies as needed)

**Python Requirements:**
- Python 3.7+
- numpy
- matplotlib
- scipy
- pandas
- astropy (recommended)
- corner (for corner plots)

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/chomskyhonk/BayesianStellarParameters.git
   cd BayesianStellarParameters
   ```

2. **Add your code:**
   - Place C++ files in `src/cpp/`
   - Place Python files in `src/python/`
   - Place your thesis/project PDF in `docs/`

3. **Install Python dependencies:**
   ```bash
   cd src/python
   pip install -r requirements.txt
   ```

4. **Build C++ code:**
   ```bash
   cd src/cpp
   # Add compilation instructions based on your code
   ```

## Usage

### Data Processing (C++)

Place usage instructions for your C++ data handling and statistical analysis code here.

Example:
```bash
cd src/cpp
./analysis input_data.csv output_results.dat
```

### Visualization (Python)

Place usage instructions for your Python visualization code here.

Example:
```python
cd src/python
python visualize_results.py --input ../path/to/results.dat
```

### Examples

See the `examples/` directory for:
- Tutorial notebooks
- Sample scripts
- Example outputs

## Project Documentation

- **Thesis/Project Paper**: See `docs/` for the complete theoretical background, methodology, and analysis
- **Code Documentation**: Each subdirectory contains a README with specific documentation


##Results

(Add summary of key findings or link to thesis)

Key achievements:
- Improved uncertainty quantification compared to frequentist methods
- Successful recovery of known Galactic trends
- Handling of isochrone degeneracies
- Identification of overabundance of old, metal-rich stars

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature/new-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{bayesian_stellar_parameters,
  author = {Your Name},
  title = {Bayesian Stellar Parameters: A Framework for Robust Parameter Estimation},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/chomskyhonk/BayesianStellarParameters}
}
```

## Contact

For questions or collaboration inquiries:
- GitHub Issues: [Project Issues](https://github.com/chomskyhonk/BayesianStellarParameters/issues)
- (Add your contact information if desired)

## Acknowledgments

- (Add acknowledgments for data sources, collaborators, funding, etc.)

## References

- (Add key references for the methodology and stellar evolution models used)
