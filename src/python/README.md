# Python Source Code

This directory contains Python code for data visualization and analysis.

## Purpose

The Python code handles:
- Visualization: Creating plots and graphs of results
- Figure generation: Producing publication-quality figures

## Setup

Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

After running the calculations, you can plot the results using the provided Python scripts. Each script is designed to visualize specific aspects of the data.

### Available Python Scripts

The following scripts are available in the `/src/python` directory:

1. Plot Priors:
   - Script: `priors.py`
   - Description: Handles the calculation and visualization of prior distributions.
   - Usage: Ensure parameters are passed within the script or modify for your data.

2. 3D Scatter Visualization:
   - Script: `3dscatter.py`
   - Description: Creates 3D scatter plots of age, metallicity bins, and velocities with dispersion data.
   - Usage: The script should reference your dataset and configurations pre-defined in the code.

3. Plot Hertzsprung-Russell Diagram (HRD):
   - Script: `HRDplot.py`
   - Description: Visualizes the HRD using stellar properties, including Bayesian ages and metallicities.
   - Usage: Replace data paths in the script and provide a `csv_file` during execution.

4. Plot Velocity Dispersion:
   - Script: `plotdispersion.py`
   - Description: Analyzes and plots velocity dispersion components based on temperature and luminosity bins.
   - Usage: 
     ```bash
     python3 src/python/plotdispersion.py
     ```

### Requirements for Plotting Scripts

Ensure you have the following Python libraries installed:
- `matplotlib`
- `numpy`
- `pandas` (for data manipulation, if applicable)

Install the required libraries using pip:
```bash
pip install matplotlib numpy pandas
```
