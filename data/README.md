## Data Requirements

The project requires specific astrophysical and isochrone data files for calculations. Due to the size of the data files, they are not included in the repository. Follow the instructions below to obtain and format the datasets.

### 1. Gaia Dataset Requirements

The main dataset is a CSV file containing astrophysical data queried from the Gaia database. You can obtain this data using the Gaia Archive or other suitable services (e.g., `astroquery`). 

#### Mandatory Columns

Ensure the dataset includes the following columns:

```
SOURCE_ID, parallax, parallax_error, distance_gspphot, distance_gspphot_lower, distance_gspphot_upper,
teff_gspphot, teff_gspphot_lower, teff_gspphot_upper, logg_gspphot, logg_gspphot_lower, logg_gspphot_upper,
bp_rp, mh_gspphot, mh_gspphot_lower, mh_gspphot_upper, ebpminrp_gspphot, phot_g_mean_mag, ra, ra_error,
dec, dec_error, pmra, pmra_error, pmdec, pmdec_error, radial_velocity, radial_velocity_error, source_id,
age_flame, age_flame_lower, age_flame_upper, mass_flame, mass_flame_lower, mass_flame_upper, lum_flame,
lum_flame_lower, lum_flame_upper, mh_gspspec, mh_gspspec_lower, mh_gspspec_upper, U, V, W, error_U, error_V, error_W
```

If your query does not include these columns, ensure to modify your query to include them.

#### Query Instructions

To query Gaia data with all the necessary fields:
- Use the [Gaia Archive Service](https://gea.esac.esa.int/archive/) or `astroquery` in Python.
- Example query (use ADQL for bulk querying):
    ```sql
    SELECT SOURCE_ID, parallax, parallax_error, distance_gspphot, distance_gspphot_lower, distance_gspphot_upper,
           teff_gspphot, teff_gspphot_lower, teff_gspphot_upper, logg_gspphot, logg_gspphot_lower, logg_gspphot_upper,
           bp_rp, mh_gspphot, mh_gspphot_lower, mh_gspphot_upper, ebpminrp_gspphot, phot_g_mean_mag, ra, ra_error,
           dec, dec_error, pmra, pmra_error, pmdec, pmdec_error, radial_velocity, radial_velocity_error, source_id,
           age_flame, age_flame_lower, age_flame_upper, mass_flame, mass_flame_lower, mass_flame_upper, lum_flame,
           lum_flame_lower, lum_flame_upper, mh_gspspec, mh_gspspec_lower, mh_gspspec_upper, U, V, W, error_U, error_V, error_W
    FROM gaiadr3.gaia_source
    ```
- Save the output as a CSV file and place it in the `data/` directory of the project.

#### Directory Placement
Place the Gaia data file in the following directory:
```
data/gaia_data.csv
```

---

### 2. Isochrone Files

The project relies on a set of isochrone data files. Isochrones must range from:
- Age: 100 Myr to 15,000 Myr (100 Myr increments)
- Metallicity: -1.5 dex to 0.42 dex (48 increments)

#### Generating Isochrone Data
If you do not have pre-existing isochrone data files, use the [MIST](http://waps.cfa.harvard.edu/MIST/) or other isochrone generators (I used BaSTI isochrones) with the following parameters:
- Age Range: 100 Myr to 14,000 Myr (in 100 Myr increments)
- [M/H] Range (dex): -1.5 to 0.42 (48 increments).

### Example Isochrone File

Below is an example of an isochrone file using data from the BaSTI-IAC database. Each file contains metadata about the isochrone (e.g., metallicity, initial/final mass, age) followed by columns of tabular data.

```plaintext
# Isochrone from BaSTI-IAC database
# Scaled solar models & transformations  -  GAIA-DR3
#========================================================================================================================================================================================================
#========================================================================================================================================================================================================
#  Np = 2100   [M/H] = -1.500   Z = 0.0004968   Y = 0.24756614   Age (Myr) = 100.000
#========================================================================================================================================================================================================
#    M/Mo(ini)     M/Mo(fin)    log(L/Lo)  logTe        G        G_BP      G_RP     G_RVS
#========================================================================================================================================================================================================
   0.1000000000   0.1000000000  -2.53793   3.54406   11.6460   12.6495   10.6629   10.2228
   0.1134406865   0.1134406865  -2.44052   3.55409   11.3646   12.3283   10.4038    9.9787
   0.1268813730   0.1268813730  -2.35297   3.56344   11.1103   12.0356   10.1710    9.7572
   0.1403220595   0.1403220595  -2.27922   3.56982   10.9007   11.7996    9.9764    9.5701
   0.1537627459   0.1537627459  -2.21084   3.57558   10.7070   11.5814    9.7963    9.3962
   0.1672034324   0.1672034324  -2.15090   3.57984   10.5402   11.3964    9.6398    9.2440
   0.1806441189   0.1806441189  -2.09416   3.58389   10.3823   11.2211    9.4916    9.0998
   0.1940848054   0.1940848054  -2.04320   3.58729   10.2411   11.0651    9.3588    8.9703
```

#### Explanation of the Columns
- **`M/Mo(ini)`**: Initial mass (solar mass).
- **`M/Mo(fin)`**: Final mass (solar mass).
- **`log(L/Lo)`**: Logarithm of luminosity relative to the Sun.
- **`logTe`**: Logarithm of effective temperature.
- **`G`, `G_BP`, `G_RP`, `G_RVS`**: Gaia magnitudes in different bands.

#### Metadata
The metadata at the top of the file contains parameters for this specific isochrone:
- **`Np`**: Number of points in the isochrone.
- **`[M/H]`**: Metallicity.
- **`Z`**: Heavy element abundance fraction.
- **`Y`**: Helium abundance.
- **`Age (Myr)`**: Age of the isochrone in millions of years.

---

Repeat this for other isochrone files or customize the example path as required.



Save each isochrone file as a separate `.dat` or `.isc_gaia-dr3` file, following the directory structure below.

Use meaningful filenames to reflect the age and metallicity of each isochrone file.
Throughout this repository, I use  functions that load the isochrone data when they are in the form:
`100z0004968y248P00O0D0E0.isc_gaia-dr3`


### Example Directory Structure

After preparing the data files, your directory structure should look like this:
```
project/
│
├── data/
│   ├── gaia_data.csv
│   └── isochrones/
│       ├── 100z0004968y248P00O0D0E0.isc_gaia-dr3
│       ├── 100z0005458y248P00O0D0E0.isc_gaia-dr3
│       ├── ...
│       ├── 14000z0363071y298P00O0D0E0.isc_gaia-dr3
│
├── src/
│   ├── cpp/
│   ├── python/
```
