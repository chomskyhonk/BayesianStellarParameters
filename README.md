# Bayesian Stellar Parameters

A Bayesian inference pipeline to estimate fundamental stellar parameters, specifically age and metallicity, for Solar neighbourhood stars.

## Overview

Traditional frequentist techniques suffer from large uncertainties and isochrone degeneracies when estimating stellar parameters. This project tests whether a Bayesian framework can improve inference and reproduce known Galactic trends.

### Key Features

- **Bayesian Inference**: Robust parameter estimation with proper uncertainty quantification
- **Statistical Analysis**: Advanced MCMC sampling and likelihood calculations
- **Data Processing**: Efficient handling of stellar parameter data
- **Visualization**: Publication-quality plots and corner plots for posterior distributions

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

**C++ Requirements:**
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
   pip install -r requirements.txt  # (create this file with your dependencies)
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
- **API Reference**: (Add links to generated documentation if available)

## Methodology

This project implements a Bayesian framework for stellar parameter estimation:

1. **Data Preparation**: Load and preprocess observational data
2. **Model Setup**: Define likelihood functions and priors based on stellar evolution models
3. **MCMC Sampling**: Use Markov Chain Monte Carlo to sample the posterior distribution
4. **Analysis**: Compute parameter estimates and uncertainties
5. **Visualization**: Generate diagnostic plots and scientific figures

## Results

(Add summary of key findings or link to thesis)

Key achievements:
- Improved uncertainty quantification compared to frequentist methods
- Successful recovery of known Galactic trends
- Robust handling of isochrone degeneracies

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature/new-feature`)
5. Open a Pull Request

## License

(Add license information)

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
