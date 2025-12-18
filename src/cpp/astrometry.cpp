#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <filesystem>
#include <map>
#include <utility>
#include <algorithm>
#include <cmath>
#include <limits>
#include <execution>
#include "matplotlibcpp.h"
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <omp.h>
#include <iomanip>
#include <regex>
#include <chrono>
#include <random>
#include <filesystem>
#include <tuple>
#include <Eigen/Dense>

namespace fs = std::filesystem;
namespace plt = matplotlibcpp;
namespace py = pybind11;

struct DispersionStats {
    double sigma_U;
    double sigma_V;
    double sigma_W;
    double sigma_total;
};

class GAIA {
public:
    std::vector<double> parallax;
    std::vector<double> parallax_error;
    std::vector<double> distance_gspphot;
    std::vector<double> distance_gspphot_lower;
    std::vector<double> distance_gspphot_upper;
    std::vector<double> teff_gspphot;
    std::vector<double> teff_gspphot_lower;
    std::vector<double> teff_gspphot_upper;
    std::vector<double> logg_gspphot;
    std::vector<double> logg_gspphot_lower;
    std::vector<double> logg_gspphot_upper;
    std::vector<double> bp_rp;
    std::vector<double> mh_gspphot;
    std::vector<double> mh_gspphot_lower;
    std::vector<double> mh_gspphot_upper;
    std::vector<double> phot_g_mean_mag;
    std::vector<double> ra;
    std::vector<double> ra_error;
    std::vector<double> dec;
    std::vector<double> dec_error;
    std::vector<double> pmra;
    std::vector<double> pmra_error;
    std::vector<double> pmdec;
    std::vector<double> pmdec_error;
    std::vector<double> radial_velocity;
    std::vector<double> radial_velocity_error;
    std::vector<long long int> source_id;
    std::vector<double> age_flame;
    std::vector<double> age_flame_lower;
    std::vector<double> age_flame_upper;
    std::vector<double> mass_flame;
    std::vector<double> mass_flame_lower;
    std::vector<double> mass_flame_upper;
    std::vector<double> lum_flame;
    std::vector<double> lum_flame_lower;
    std::vector<double> lum_flame_upper;
    std::vector<double> mh_gspspec;
    std::vector<double> mh_gspspec_lower;
    std::vector<double> mh_gspspec_upper;
    std::vector<double> U;
    std::vector<double> V;
    std::vector<double> W;
    std::vector<double> age_mean; // Added for mean age
    std::vector<double> age_std; // Added for standard deviation of age
    std::vector<double> met_mean; // Added for mean metallicity
    std::vector<double> met_std; // Added for standard deviation of metallicity
    std::vector<double> error_U;
    std::vector<double> error_V;
    std::vector<double> error_W;

    GAIA() = default;
    GAIA (std::vector<double> parallax, std::vector<double> parallax_error, std::vector<double> distance_gspphot, std::vector<double> distance_gspphot_lower, std::vector<double> distance_gspphot_upper, std::vector<double> teff_gspphot, std::vector<double>  teff_gspphot_lower, std::vector<double> teff_gspphot_upper, std::vector<double> logg_gspphot, std::vector<double> logg_gspphot_lower, std::vector<double> logg_gspphot_upper, std::vector<double> bp_rp, std::vector<double> mh_gspphot, std::vector<double> mh_gspphot_lower, std::vector<double> mh_gspphot_upper, std::vector<double> phot_g_mean_mag, std::vector<double> ra, std::vector<double> ra_error, std::vector<double> dec, std::vector<double> dec_error, std::vector<double> pmra, std::vector<double> pmra_error, std::vector<double> pmdec, std::vector<double> pmdec_error, std::vector<double> radial_velocity, std::vector<double> radial_velocity_error, std::vector<long long int> source_id, std::vector<double> age_flame, std::vector<double> age_flame_lower, std::vector<double> age_flame_upper, std::vector<double> mass_flame, std::vector<double> mass_flame_lower, std::vector<double> mass_flame_upper, std::vector<double> lum_flame, std::vector<double> lum_flame_lower, std::vector<double> lum_flame_upper, std::vector<double> mh_gspspec, std::vector<double> mh_gspspec_lower, std::vector<double> mh_gspspec_upper, std::vector<double> U, std::vector<double> V, std::vector<double> W, std::vector<double> age_mean, std::vector<double> age_std, std::vector<double> met_mean, std::vector<double> met_std, std::vector<double> error_U, std::vector<double> error_V, std::vector<double> error_W)
        : parallax(parallax), parallax_error(parallax_error), distance_gspphot(distance_gspphot), distance_gspphot_lower(distance_gspphot_lower), distance_gspphot_upper(distance_gspphot_upper), teff_gspphot(teff_gspphot), teff_gspphot_lower(teff_gspphot_lower), teff_gspphot_upper(teff_gspphot_upper), logg_gspphot(logg_gspphot), logg_gspphot_lower(logg_gspphot_lower), logg_gspphot_upper(logg_gspphot_upper), bp_rp(bp_rp), mh_gspphot(mh_gspphot), mh_gspphot_lower(mh_gspphot_lower), mh_gspphot_upper(mh_gspphot_upper), phot_g_mean_mag(phot_g_mean_mag), ra(ra), ra_error(ra_error), dec(dec), dec_error(dec_error), pmra(pmra), pmra_error(pmra_error), pmdec(pmdec), pmdec_error(pmdec_error), radial_velocity(radial_velocity), radial_velocity_error(radial_velocity_error), source_id(source_id), age_flame(age_flame), age_flame_lower(age_flame_lower), age_flame_upper(age_flame_upper), mass_flame(mass_flame), mass_flame_lower(mass_flame_lower), mass_flame_upper(mass_flame_upper), lum_flame(lum_flame), lum_flame_lower(lum_flame_lower), lum_flame_upper(lum_flame_upper), mh_gspspec(mh_gspspec), mh_gspspec_lower(mh_gspspec_lower), mh_gspspec_upper(mh_gspspec_upper), U(U), V(V), W(W), age_mean(age_mean), age_std(age_std), met_mean(met_mean), met_std(met_std), error_U(error_U), error_V(error_V), error_W(error_W) {}

    void addData(double parallax, double parallax_error, double distance_gspphot, double distance_gspphot_lower, double distance_gspphot_upper, double teff_gspphot, double teff_gspphot_lower, double teff_gspphot_upper, double logg_gspphot, double logg_gspphot_lower, double logg_gspphot_upper, double bp_rp, double mh_gspphot, double mh_gspphot_lower, double mh_gspphot_upper, double phot_g_mean_mag, double ra, double ra_error, double dec, double dec_error, double pmra, double pmra_error, double pmdec, double pmdec_error, double radial_velocity, double radial_velocity_error, long long int source_id, double age_flame, double age_flame_lower, double age_flame_upper, double mass_flame, double mass_flame_lower, double mass_flame_upper, double lum_flame, double lum_flame_lower, double lum_flame_upper, double mh_gspspec, double mh_gspspec_lower, double mh_gspspec_upper, double U, double V, double W, double age_mean, double age_std, double met_mean, double met_std, double error_U, double error_V, double error_W)
    {
        this->parallax.push_back(parallax);
        this->parallax_error.push_back(parallax_error);
        this->distance_gspphot.push_back(distance_gspphot);
        this->distance_gspphot_lower.push_back(distance_gspphot_lower);
        this->distance_gspphot_upper.push_back(distance_gspphot_upper);
        this->teff_gspphot.push_back(teff_gspphot);
        this->teff_gspphot_lower.push_back(teff_gspphot_lower);
        this->teff_gspphot_upper.push_back(teff_gspphot_upper);
        this->logg_gspphot.push_back(logg_gspphot);
        this->logg_gspphot_lower.push_back(logg_gspphot_lower);
        this->logg_gspphot_upper.push_back(logg_gspphot_upper);
        this->bp_rp.push_back(bp_rp);
        this->mh_gspphot.push_back(mh_gspphot);
        this->mh_gspphot_lower.push_back(mh_gspphot_lower);
        this->mh_gspphot_upper.push_back(mh_gspphot_upper);
        this->phot_g_mean_mag.push_back(phot_g_mean_mag);
        this->ra.push_back(ra);
        this->ra_error.push_back(ra_error);
        this->dec.push_back(dec);
        this->dec_error.push_back(dec_error);
        this->pmra.push_back(pmra);
        this->pmra_error.push_back(pmra_error);
        this->pmdec.push_back(pmdec);
        this->pmdec_error.push_back(pmdec_error);
        this->radial_velocity.push_back(radial_velocity);
        this->radial_velocity_error.push_back(radial_velocity_error);
        this->source_id.push_back(source_id);
        this->age_flame.push_back(age_flame);
        this->age_flame_lower.push_back(age_flame_lower);
        this->age_flame_upper.push_back(age_flame_upper);
        this->mass_flame.push_back(mass_flame);
        this->mass_flame_lower.push_back(mass_flame_lower);
        this->mass_flame_upper.push_back(mass_flame_upper);
        this->lum_flame.push_back(lum_flame);
        this->lum_flame_lower.push_back(lum_flame_lower);
        this->lum_flame_upper.push_back(lum_flame_upper);
        this->mh_gspspec.push_back(mh_gspspec);
        this->mh_gspspec_lower.push_back(mh_gspspec_lower);
        this->mh_gspspec_upper.push_back(mh_gspspec_upper);
        this->U.push_back(U);
        this->V.push_back(V);
        this->W.push_back(W);
        this->age_mean.push_back(age_mean);
        this->age_std.push_back(age_std);
        this->met_mean.push_back(met_mean);
        this->met_std.push_back(met_std);
        this->error_U.push_back(error_U);
        this->error_V.push_back(error_V);
        this->error_W.push_back(error_W);
    }
    void loadData(const std::string& filename);
    std::vector<size_t> filterData(double logTe, double logL, bool useAllStars);
    std::vector<std::pair<double, double>> extractAgeAndMetallicity(const std::vector<size_t>& indices);
    
};

void GAIA::loadData(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line); // Read the header line

    // Parse the header line to get column indices
    std::unordered_map<std::string, int> columnIndices;
    std::stringstream headerStream(line);
    std::string columnName;
    int index = 0;
    bool firstSourceId = true;
    while (std::getline(headerStream, columnName, ',')) {
        if (columnName == "source_id") {
            if (firstSourceId) {
                firstSourceId = false;
            } else {
                columnName = "source_id_2"; // Rename the second source_id column
            }
        }
        columnIndices[columnName] = index++;
    }

    // Read the data lines
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<std::string> values;

        // Split the line into values
        while (std::getline(ss, value, ',')) {
            values.push_back(value);
        }

        // Check if the number of values matches the number of columns
        if (values.size() != columnIndices.size()) {
            while (values.size() < columnIndices.size()) {
                values.push_back(""); // Add empty strings for missing values
            }
        }

        try {
            parallax.push_back(values[columnIndices["parallax"]].empty() ? std::nan("") : std::stod(values[columnIndices["parallax"]]));
            parallax_error.push_back(values[columnIndices["parallax_error"]].empty() ? std::nan("") : std::stod(values[columnIndices["parallax_error"]]));
            distance_gspphot.push_back(values[columnIndices["distance_gspphot"]].empty() ? std::nan("") : std::stod(values[columnIndices["distance_gspphot"]]));
            distance_gspphot_lower.push_back(values[columnIndices["distance_gspphot_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["distance_gspphot_lower"]]));
            distance_gspphot_upper.push_back(values[columnIndices["distance_gspphot_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["distance_gspphot_upper"]]));
            teff_gspphot.push_back(values[columnIndices["teff_gspphot"]].empty() ? std::nan("") : std::stod(values[columnIndices["teff_gspphot"]]));
            teff_gspphot_lower.push_back(values[columnIndices["teff_gspphot_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["teff_gspphot_lower"]]));
            teff_gspphot_upper.push_back(values[columnIndices["teff_gspphot_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["teff_gspphot_upper"]]));
            logg_gspphot.push_back(values[columnIndices["logg_gspphot"]].empty() ? std::nan("") : std::stod(values[columnIndices["logg_gspphot"]]));
            logg_gspphot_lower.push_back(values[columnIndices["logg_gspphot_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["logg_gspphot_lower"]]));
            logg_gspphot_upper.push_back(values[columnIndices["logg_gspphot_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["logg_gspphot_upper"]]));
            bp_rp.push_back(values[columnIndices["bp_rp"]].empty() ? std::nan("") : std::stod(values[columnIndices["bp_rp"]]));
            mh_gspphot.push_back(values[columnIndices["mh_gspphot"]].empty() ? std::nan("") : std::stod(values[columnIndices["mh_gspphot"]]));
            mh_gspphot_lower.push_back(values[columnIndices["mh_gspphot_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["mh_gspphot_lower"]]));
            mh_gspphot_upper.push_back(values[columnIndices["mh_gspphot_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["mh_gspphot_upper"]]));
            phot_g_mean_mag.push_back(values[columnIndices["phot_g_mean_mag"]].empty() ? std::nan("") : std::stod(values[columnIndices["phot_g_mean_mag"]]));
            ra.push_back(values[columnIndices["ra"]].empty() ? std::nan("") : std::stod(values[columnIndices["ra"]]));
            ra_error.push_back(values[columnIndices["ra_error"]].empty() ? std::nan("") : std::stod(values[columnIndices["ra_error"]]));
            dec.push_back(values[columnIndices["dec"]].empty() ? std::nan("") : std::stod(values[columnIndices["dec"]]));
            dec_error.push_back(values[columnIndices["dec_error"]].empty() ? std::nan("") : std::stod(values[columnIndices["dec_error"]]));
            pmra.push_back(values[columnIndices["pmra"]].empty() ? std::nan("") : std::stod(values[columnIndices["pmra"]]));
            pmra_error.push_back(values[columnIndices["pmra_error"]].empty() ? std::nan("") : std::stod(values[columnIndices["pmra_error"]]));
            pmdec.push_back(values[columnIndices["pmdec"]].empty() ? std::nan("") : std::stod(values[columnIndices["pmdec"]]));
            pmdec_error.push_back(values[columnIndices["pmdec_error"]].empty() ? std::nan("") : std::stod(values[columnIndices["pmdec_error"]]));
            radial_velocity.push_back(values[columnIndices["radial_velocity"]].empty() ? std::nan("") : std::stod(values[columnIndices["radial_velocity"]]));
            radial_velocity_error.push_back(values[columnIndices["radial_velocity_error"]].empty() ? std::nan("") : std::stod(values[columnIndices["radial_velocity_error"]]));
            source_id.push_back(values.at(columnIndices.at("source_id")).empty() ? -1 : std::stoll(values.at(columnIndices.at("source_id"))));
            age_flame.push_back(values[columnIndices["age_flame"]].empty() ? std::nan("") : std::stod(values[columnIndices["age_flame"]]));
            age_flame_lower.push_back(values[columnIndices["age_flame_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["age_flame_lower"]]));
            age_flame_upper.push_back(values[columnIndices["age_flame_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["age_flame_upper"]]));
            mass_flame.push_back(values[columnIndices["mass_flame"]].empty() ? std::nan("") : std::stod(values[columnIndices["mass_flame"]]));
            mass_flame_lower.push_back(values[columnIndices["mass_flame_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["mass_flame_lower"]]));
            mass_flame_upper.push_back(values[columnIndices["mass_flame_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["mass_flame_upper"]]));
            lum_flame.push_back(values[columnIndices["lum_flame"]].empty() ? std::nan("") : std::stod(values[columnIndices["lum_flame"]]));
            lum_flame_lower.push_back(values[columnIndices["lum_flame_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["lum_flame_lower"]]));
            lum_flame_upper.push_back(values[columnIndices["lum_flame_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["lum_flame_upper"]]));
            mh_gspspec.push_back(values[columnIndices["mh_gspspec"]].empty() ? std::nan("") : std::stod(values[columnIndices["mh_gspspec"]]));
            mh_gspspec_lower.push_back(values[columnIndices["mh_gspspec_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["mh_gspspec_lower"]]));
            mh_gspspec_upper.push_back(values[columnIndices["mh_gspspec_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["mh_gspspec_upper"]]));
            U.push_back(values[columnIndices["U"]].empty() ? std::nan("") : std::stod(values[columnIndices["U"]]));
            V.push_back(values[columnIndices["V"]].empty() ? std::nan("") : std::stod(values[columnIndices["V"]]));
            W.push_back(values[columnIndices["W"]].empty() ? std::nan("") : std::stod(values[columnIndices["W"]]));
            age_mean.push_back(values[columnIndices["age_mean"]].empty() ? std::nan("") : std::stod(values[columnIndices["age_mean"]]));
            age_std.push_back(values[columnIndices["age_std"]].empty() ? std::nan("") : std::stod(values[columnIndices["age_std"]]));
            met_mean.push_back(values[columnIndices["met_mean"]].empty() ? std::nan("") : std::stod(values[columnIndices["met_mean"]]));
            met_std.push_back(values[columnIndices["met_std"]].empty() ? std::nan("") : std::stod(values[columnIndices["met_std"]]));
            // Load error_U, error_V, error_W if present
            if (columnIndices.count("error_U")) {
                error_U.push_back(values[columnIndices["error_U"]].empty() ? std::nan("") : std::stod(values[columnIndices["error_U"]]));
            }
            if (columnIndices.count("error_V")) {
                error_V.push_back(values[columnIndices["error_V"]].empty() ? std::nan("") : std::stod(values[columnIndices["error_V"]]));
            }
            if (columnIndices.count("error_W")) {
                error_W.push_back(values[columnIndices["error_W"]].empty() ? std::nan("") : std::stod(values[columnIndices["error_W"]]));
            }
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument: " << e.what() << " in line: " << line << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range: " << e.what() << " in line: " << line << std::endl;
        }
    }

    file.close();
}

double converttodistance(double parallax) {
    if (parallax == 0) {
        return std::numeric_limits<double>::infinity(); // Handle zero parallax
    } //maybe add filter for small parallaxes here as well e.g. <=0.2?

    return 1000.0 / parallax; //convert parallax in milliarcseconds to distance in parsecs
}

std::tuple<double, double, double> compute3dspacevelocity(
    double ra, double dec, double radial_velocity, double pmra, double pmdec, double parallax) {
    
    try {
        py::module astropy_coords = py::module::import("astropy.coordinates");
        py::module astropy_units = py::module::import("astropy.units");
        py::module numpy = py::module::import("numpy");
      

        double distance = converttodistance(parallax);

        auto skycoord = astropy_coords.attr("SkyCoord")(
            py::arg("ra") = py::float_(ra) * astropy_units.attr("deg"),
            py::arg("dec") = py::float_(dec) * astropy_units.attr("deg"),
            py::arg("distance") = py::float_(distance) * astropy_units.attr("pc"),
            py::arg("pm_ra_cosdec") = py::float_(pmra) * astropy_units.attr("mas") / astropy_units.attr("yr"),
            py::arg("pm_dec") = py::float_(pmdec) * astropy_units.attr("mas") / astropy_units.attr("yr"),
            py::arg("radial_velocity") = py::float_(radial_velocity) * astropy_units.attr("km") / astropy_units.attr("s")
        );

        auto galactocentric = astropy_coords.attr("Galactocentric")();
        auto gc = skycoord.attr("transform_to")(galactocentric);
        auto v_xyz = gc.attr("velocity").attr("d_xyz").attr("to")("km/s");

        py::array_t<double> v_xyz_array = v_xyz.cast<py::array_t<double>>();
        auto buf = v_xyz_array.request();
        double* ptr = static_cast<double*>(buf.ptr);

        return std::make_tuple(ptr[0], ptr[1], ptr[2]);

    } catch (const py::error_already_set& e) {
        std::cerr << "Python exception: " << e.what() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "C++ exception: " << e.what() << std::endl;
    }

    // Return NaNs if something went wrong
    return std::make_tuple(std::nan(""), std::nan(""), std::nan(""));
}

// Compute uncertainties in U, V, W using error propagation (Johnson & Soderblom 1987)
// Inputs: ra, dec, parallax, pmra, pmdec, radial_velocity and their errors
// Returns: std::tuple<double, double, double> (sigma_U, sigma_V, sigma_W)
std::tuple<double, double, double> compute3dspacevelocityerror(
    double ra, double dec, double parallax, double pmra, double pmdec, double radial_velocity,
    double ra_error, double dec_error, double parallax_error, double pmra_error, double pmdec_error, double radial_velocity_error)
{
    // Convert angles to radians
    double ra_rad = ra * M_PI / 180.0;
    double dec_rad = dec * M_PI / 180.0;

    // Convert parallax to distance in parsec
    double d = converttodistance(parallax);
    double d_err = (parallax != 0) ? (1000.0 * parallax_error) / (parallax * parallax) : std::nan("");

    // Proper motions in arcsec/yr
    double pmra_sec = pmra / 1000.0 / 3600.0;
    double pmdec_sec = pmdec / 1000.0 / 3600.0;
    double pmra_err_sec = pmra_error / 1000.0 / 3600.0;
    double pmdec_err_sec = pmdec_error / 1000.0 / 3600.0;

    // Conversion factor: 4.74047 = 1 AU/yr in km/s
    double k = 4.74047;

    // Space velocity components in equatorial system
    double v_alpha = k * pmra / parallax;
    double v_delta = k * pmdec / parallax;
    double v_r = radial_velocity;

    // Partial derivatives for error propagation (see Johnson & Soderblom 1987)
    // For U, V, W, we use the transformation matrix (see their eq. 1-3)
    double cos_ra = std::cos(ra_rad), sin_ra = std::sin(ra_rad);
    double cos_dec = std::cos(dec_rad), sin_dec = std::sin(dec_rad);

    // Transformation matrix elements
    double l1 = -sin_ra;
    double l2 = cos_ra;
    double l3 = 0.0;
    double m1 = -sin_dec * cos_ra;
    double m2 = -sin_dec * sin_ra;
    double m3 = cos_dec;
    double n1 = cos_dec * cos_ra;
    double n2 = cos_dec * sin_ra;
    double n3 = sin_dec;

    // U, V, W (right-handed, U toward Galactic center)
    double U = l1 * v_alpha + m1 * v_delta + n1 * v_r;
    double V = l2 * v_alpha + m2 * v_delta + n2 * v_r;
    double W = l3 * v_alpha + m3 * v_delta + n3 * v_r;

    // Partial derivatives for error propagation
    // dU/dpmra, dU/dpmdec, dU/dparallax, dU/drv
    double dU_dpmra = l1 * k / parallax;
    double dU_dpmdec = m1 * k / parallax;
    double dU_dparallax = -k * (l1 * pmra + m1 * pmdec) / (parallax * parallax);
    double dU_drv = n1;

    double dV_dpmra = l2 * k / parallax;
    double dV_dpmdec = m2 * k / parallax;
    double dV_dparallax = -k * (l2 * pmra + m2 * pmdec) / (parallax * parallax);
    double dV_drv = n2;

    double dW_dpmra = l3 * k / parallax;
    double dW_dpmdec = m3 * k / parallax;
    double dW_dparallax = -k * (l3 * pmra + m3 * pmdec) / (parallax * parallax);
    double dW_drv = n3;

    // Error propagation (neglecting covariances)
    double sigma_U = std::sqrt(
        std::pow(dU_dpmra * pmra_error, 2) +
        std::pow(dU_dpmdec * pmdec_error, 2) +
        std::pow(dU_dparallax * parallax_error, 2) +
        std::pow(dU_drv * radial_velocity_error, 2)
    );
    double sigma_V = std::sqrt(
        std::pow(dV_dpmra * pmra_error, 2) +
        std::pow(dV_dpmdec * pmdec_error, 2) +
        std::pow(dV_dparallax * parallax_error, 2) +
        std::pow(dV_drv * radial_velocity_error, 2)
    );
    double sigma_W = std::sqrt(
        std::pow(dW_dpmra * pmra_error, 2) +
        std::pow(dW_dpmdec * pmdec_error, 2) +
        std::pow(dW_dparallax * parallax_error, 2) +
        std::pow(dW_drv * radial_velocity_error, 2)
    );

    return std::make_tuple(sigma_U, sigma_V, sigma_W);
}

void saveCSVWithUVW(const std::string& inputFile, const std::string& outputFile,
                    const std::vector<std::tuple<double, double, double>>& uvw) {
    
    std::ifstream in(inputFile);
    std::ofstream out(outputFile);

    std::string line;
    std::getline(in, line);  // header
    out << line << ",U,V,W\n";  // write new header

    size_t i = 0;
    while (std::getline(in, line)) {
        out << line;

        if (i < uvw.size()) {
            auto [U, V, W] = uvw[i++];
            out << "," << U << "," << V << "," << W;
        }

        out << "\n";
    }

    in.close();
    out.close();
}

std::vector<std::tuple<double, double, double>> printUVW(const GAIA& gaiaData, size_t max_stars = 100) {
    std::vector<std::tuple<double, double, double>> uvw_results;

    size_t n = std::min(max_stars, gaiaData.parallax.size());
    for (size_t i = 0; i < n; ++i) {
        double parallax = gaiaData.parallax[i];
        double radial_velocity = gaiaData.radial_velocity[i];
        double pmra = gaiaData.pmra[i];
        double pmdec = gaiaData.pmdec[i];
        double ra = gaiaData.ra[i];
        double dec = gaiaData.dec[i];

        // Skip if any required value is NaN
        if (std::isnan(parallax) || std::isnan(radial_velocity) ||
            std::isnan(pmra) || std::isnan(pmdec)) {
            uvw_results.emplace_back(NAN, NAN, NAN);
            continue;
        }

        try {
            auto [U, V, W] = compute3dspacevelocity(
                ra, dec, radial_velocity, pmra, pmdec, parallax
            );
            std::cout << "Star " << i << ": U=" << U << " V=" << V << " W=" << W << std::endl;
            uvw_results.emplace_back(U, V, W);
        } catch (const std::exception& e) {
            std::cerr << "Error computing space velocity for star " << i
                      << " (ra=" << ra << ", dec=" << dec
                      << ", parallax=" << parallax
                      << ", pmra=" << pmra
                      << ", pmdec=" << pmdec
                      << ", rv=" << radial_velocity << "): "
                      << e.what() << std::endl;
            uvw_results.emplace_back(NAN, NAN, NAN);
        }
    }

    return uvw_results;
}

using DispersionMap = std::map<std::pair<int, int>, DispersionStats>;

DispersionMap computeVelocityDispersionBinned(
    const GAIA& gaiaData,
    double param1_start, double param1_end, size_t num_param1_bins,
    double param2_start, double param2_end, size_t num_param2_bins,
    std::function<double(size_t)> get_param1,
    std::function<double(size_t)> get_param2,
    const std::string& output_filename  // optional output file
) {
    // Map: (bin1, bin2) -> vector of (U, V, W)
    std::map<std::pair<int, int>, std::vector<std::tuple<double, double, double>>> bin_map;

    double param1_step = (param1_end - param1_start) / num_param1_bins;
    double param2_step = (param2_end - param2_start) / num_param2_bins;

    size_t n = gaiaData.U.size();

    for (size_t i = 0; i < n; ++i) {
        if (std::isnan(gaiaData.U[i]) || std::isnan(gaiaData.V[i]) || std::isnan(gaiaData.W[i]))
            continue;

        double val1 = get_param1(i);
        double val2 = get_param2(i);
        if (std::isnan(val1) || std::isnan(val2))
            continue;

        int bin1 = static_cast<int>((val1 - param1_start) / param1_step);
        int bin2 = static_cast<int>((val2 - param2_start) / param2_step);

        if (bin1 < 0 || bin1 >= static_cast<int>(num_param1_bins) ||
            bin2 < 0 || bin2 >= static_cast<int>(num_param2_bins)) {
            continue;
        }

        bin_map[{bin1, bin2}].emplace_back(gaiaData.U[i], gaiaData.V[i], gaiaData.W[i]);
    }

    // Print how many stars fall within each bin
    std::cout << "Stars per bin (bin1, bin2): count" << std::endl;
    size_t total_star_count = 0;
    std::map<std::pair<int, int>, size_t> bin_counts;
    for (const auto& [bin, velocities] : bin_map) {
        std::cout << "(" << bin.first << ", " << bin.second << "): " << velocities.size() << std::endl;
        total_star_count += velocities.size();
        bin_counts[bin] = velocities.size();
    }
    std::cout << "Total stars in all bins: " << total_star_count << std::endl;

    // We'll store the uncertainty as a separate map: (bin1, bin2) -> (delta_sigma_U, delta_sigma_V, delta_sigma_W, delta_sigma_total)
    std::map<std::pair<int, int>, std::tuple<double, double, double, double>> uncertainty_map;

    DispersionMap dispersion_result;
    for (const auto& [bin, velocities] : bin_map) {
        double sum_U = 0, sum_V = 0, sum_W = 0;
        size_t count = velocities.size();

        for (const auto& [U, V, W] : velocities) {
            sum_U += U;
            sum_V += V;
            sum_W += W;
        }

        double mean_U = sum_U / count;
        double mean_V = sum_V / count;
        double mean_W = sum_W / count;

        double sq_sum_U = 0, sq_sum_V = 0, sq_sum_W = 0;
        for (const auto& [U, V, W] : velocities) {
            sq_sum_U += std::pow(U - mean_U, 2);
            sq_sum_V += std::pow(V - mean_V, 2);
            sq_sum_W += std::pow(W - mean_W, 2);
        }

        double sigma_U = std::sqrt(sq_sum_U / count);
        double sigma_V = std::sqrt(sq_sum_V / count);
        double sigma_W = std::sqrt(sq_sum_W / count);
        double sigma_total = std::sqrt(sigma_U * sigma_U + sigma_V * sigma_V + sigma_W * sigma_W);

        // Uncertainty: delta_sigma = sigma / sqrt(2(N-1))
        double delta_sigma_U = (count > 1) ? sigma_U / std::sqrt(2.0 * (count - 1)) : std::nan("");
        double delta_sigma_V = (count > 1) ? sigma_V / std::sqrt(2.0 * (count - 1)) : std::nan("");
        double delta_sigma_W = (count > 1) ? sigma_W / std::sqrt(2.0 * (count - 1)) : std::nan("");
        double delta_sigma_total = (count > 1) ? sigma_total / std::sqrt(2.0 * (count - 1)) : std::nan("");

        dispersion_result[bin] = {sigma_U, sigma_V, sigma_W, sigma_total};
        uncertainty_map[bin] = {delta_sigma_U, delta_sigma_V, delta_sigma_W, delta_sigma_total};
    }
    if (!output_filename.empty()) {
        std::filesystem::path out_path(output_filename);
        if (!std::filesystem::exists(out_path.parent_path())) {
            std::filesystem::create_directories(out_path.parent_path());
        }
        std::ofstream ofs(output_filename);
        if (!ofs) {
            std::cerr << "Error opening file for writing: " << output_filename << "\n";
        } else {
            // Add bin boundaries columns to the header, and uncertainty columns, and star count
            ofs << "bin1,bin2,param1_lower,param1_upper,param1_center,"
                "param2_lower,param2_upper,param2_center,"
                "sigma_U,sigma_V,sigma_W,sigma_total,"
                "delta_sigma_U,delta_sigma_V,delta_sigma_W,delta_sigma_total,"
                "n_stars\n";

            for (const auto& [bin, stats] : dispersion_result) {
                int bin1 = bin.first;
                int bin2 = bin.second;

                double param1_lower = param1_start + bin1 * param1_step;
                double param1_upper = param1_lower + param1_step;
                double param1_center = param1_lower + 0.5 * param1_step;

                double param2_lower = param2_start + bin2 * param2_step;
                double param2_upper = param2_lower + param2_step;
                double param2_center = param2_lower + 0.5 * param2_step;

                auto [delta_sigma_U, delta_sigma_V, delta_sigma_W, delta_sigma_total] = uncertainty_map.at(bin);

                size_t n_stars = bin_counts[bin];

                ofs << bin1 << "," << bin2 << ","
                    << param1_lower << "," << param1_upper << "," << param1_center << ","
                    << param2_lower << "," << param2_upper << "," << param2_center << ","
                    << stats.sigma_U << "," << stats.sigma_V << ","
                    << stats.sigma_W << "," << stats.sigma_total << ","
                    << delta_sigma_U << "," << delta_sigma_V << ","
                    << delta_sigma_W << "," << delta_sigma_total << ","
                    << n_stars << "\n";
            }
            ofs.close();
        }
    }
    return dispersion_result;
}

// Compute velocity dispersion in 1D bins (either age or metallicity), with star counts and uncertainties
// Now also propagates errors on UVW if available (error_U, error_V, error_W)
// param: "age" or "feh" (metallicity)
// Returns: map bin_index -> DispersionStats, and saves CSV if output_filename is not empty
std::map<int, DispersionStats> computeVelocityDispersionBinned1D(
    const GAIA& gaiaData,
    double param_start, double param_end, size_t num_bins,
    std::function<double(size_t)> get_param,
    const std::string& param_name, // "age" or "feh"
    const std::string& output_filename = "",
    const std::vector<double>* error_U_ptr = nullptr,
    const std::vector<double>* error_V_ptr = nullptr,
    const std::vector<double>* error_W_ptr = nullptr,
    bool require_age_30pct = true, //only accept stars with age uncertainty <= 30%
    bool require_errorUVW_lt2 = true, //only accept stars with error_U/V/W < 2 km/s
    bool require_parallax_over_error_gt10 = true,
    bool require_vtan_ltval = true,//only accept stars with parallax/parallax_error > 10
    bool filter_teff_lum = true, //enable Teff/Luminosity filter
    double teff_min = 5000.0,    //minimum Teff (K)
    double lum_min = 0.1,        //minimum Luminosity
    bool filter_colour_mag = true, //enable colour/magnitude filter
    double bp_rp_min = 0.5,
    double bp_rp_max = 1.5,
    double M_G_min = 4.0,
    double M_G_max = 7.0,
    double max_age = 8.0, // Gyr
    bool require_total_velocity_lt150 = true, // NEW: filter on total velocity
    bool filtermaxage = true //filters out old stars 
) {
    std::map<int, std::vector<std::tuple<double, double, double, double, double, double>>> bin_map;
    double param_step = (param_end - param_start) / num_bins;
    size_t n = gaiaData.age_flame.size();
    std::cout << "Age_flame_lower size: " << gaiaData.age_flame_lower.size() << "\n";
    std::cout << "Age_flame_upper size: " << gaiaData.age_flame_upper.size() << "\n";
    std::cout << "Age_flame size: " << gaiaData.age_flame.size() << "\n";
    std::cout << "total stars: " << gaiaData.source_id.size() << "\n";
    std::cout << "parallax size: " << gaiaData.parallax.size() << "\n";
    std::cout << "parallax_error size: " << gaiaData.parallax_error.size() << "\n";

    bool use_errors = (error_U_ptr && error_V_ptr && error_W_ptr &&
                       error_U_ptr->size() == n && error_V_ptr->size() == n && error_W_ptr->size() == n);

    // Counters for diagnostics
    static size_t total_checked = 0, total_missing = 0, total_failed_uncertainty = 0, total_passed = 0;
    static size_t total_checked_err = 0, total_failed_err = 0, total_passed_err = 0;
    static size_t total_checked_parallax = 0, total_failed_parallax = 0, total_passed_parallax = 0;
    static size_t total_checked_vtan = 0, total_failed_vtan = 0, total_passed_vtan = 0;
    static size_t total_checked_tefflum = 0, total_failed_tefflum = 0, total_passed_tefflum = 0;
    static size_t total_checked_colmag = 0, total_failed_colmag = 0, total_passed_colmag = 0;
    static size_t total_checked_totalvel = 0, total_failed_totalvel = 0, total_passed_totalvel = 0;
    static size_t total_checked_maxage = 0, total_failed_maxage = 0, total_passed_maxage = 0;

    for (size_t i = 0; i < n; ++i) {
        if (std::isnan(gaiaData.U[i]) || std::isnan(gaiaData.V[i]) || std::isnan(gaiaData.W[i]))
            continue;
        double val = get_param(i);
        if (std::isnan(val)) continue;


        // Filter: only use stars with age < max_age (e.g., 8 Gyr) if enabled
        double max_age = 8.0; // Gyr
        bool filter_max_age = true; // set to false to disable this filter
        if (filter_max_age) {
            ++total_checked_maxage;
            double age = gaiaData.age_flame[i];
            if (std::isnan(age) || age >= max_age) {
            ++total_failed_maxage;
            continue;
            }
            ++total_passed_maxage;
        }
        // NEW: Teff/Luminosity filter
        if (filter_teff_lum) {
            ++total_checked_tefflum;
            bool fail = false;
            if (i >= gaiaData.teff_gspphot.size() || i >= gaiaData.lum_flame.size()) {
                fail = true;
            } else {
                double teff = gaiaData.teff_gspphot[i];
                double lum = gaiaData.lum_flame[i];
                if (std::isnan(teff) || std::isnan(lum) || teff < teff_min || lum < lum_min) {
                    fail = true;
                }
            }
            if (fail) {
                ++total_failed_tefflum;
                continue;
            }
            ++total_passed_tefflum;
        }

        // NEW: Colour/Magnitude filter
        if (filter_colour_mag) {
            ++total_checked_colmag;
            bool fail = false;
            if (i >= gaiaData.bp_rp.size() || i >= gaiaData.phot_g_mean_mag.size() || i >= gaiaData.parallax.size()) {
                fail = true;
            } else {
                double bp_rp = gaiaData.bp_rp[i];
                double phot_g = gaiaData.phot_g_mean_mag[i];
                double parallax = gaiaData.parallax[i];
                if (std::isnan(bp_rp) || std::isnan(phot_g) || std::isnan(parallax) || parallax <= 0) {
                    fail = true;
                } else {
                    double M_G = phot_g + 5.0 * std::log10(parallax / 1000.0);
                    if (bp_rp <= bp_rp_min || bp_rp >= bp_rp_max ||
                        M_G <= M_G_min || M_G >= M_G_max) {
                        fail = true;
                    }
                }
            }
            if (fail) {
                ++total_failed_colmag;
                continue;
            }
            ++total_passed_colmag;
        }

        // Parallax uncertainty cut
        if (require_parallax_over_error_gt10) {
            ++total_checked_parallax;
            if (i >= gaiaData.parallax.size() || i >= gaiaData.parallax_error.size()) {
                ++total_failed_parallax;
                continue;
            }
            double parallax = gaiaData.parallax[i];
            double parallax_error = gaiaData.parallax_error[i];
            double parallax_over_error = parallax / parallax_error;

            if (std::isnan(parallax) || std::isnan(parallax_error) || parallax_error == 0) {
                ++total_failed_parallax;
                continue;
            }
            if (i<10) {
                std::cout << "DEBUG [" << i << "]: parallax = " << parallax
                          << ", parallax_error = " << parallax_error << ", parallax_over_error = " << parallax_over_error << "\n";
            }

            if (parallax_over_error <= 100.0) {
                ++total_failed_parallax;
                continue;
            }
            ++total_passed_parallax;
        }

        // If require_age_30pct is set, check age uncertainty
        if (require_age_30pct) {
            ++total_checked;

            bool missing_age = (i >= gaiaData.age_flame.size());
            bool missing_lower = (i >= gaiaData.age_flame_lower.size());
            bool missing_upper = (i >= gaiaData.age_flame_upper.size());

            double age = missing_age ? std::nan("") : gaiaData.age_flame[i];
            double age_lower = missing_lower ? std::nan("") : gaiaData.age_flame_lower[i];
            double age_upper = missing_upper ? std::nan("") : gaiaData.age_flame_upper[i];

            // Enhanced DEBUG PRINT: print for first N failures
            static size_t debug_print_count = 0;
            constexpr size_t max_debug_prints = 20;
            if (debug_print_count < max_debug_prints) {
                std::cout << "DEBUG [" << i << "]: "
                    << "age_flame=" << (missing_age ? "MISSING" : std::to_string(age))
                    << ", age_flame_lower=" << (missing_lower ? "MISSING" : std::to_string(age_lower))
                    << ", age_flame_upper=" << (missing_upper ? "MISSING" : std::to_string(age_upper))
                    << ", isNaN(age)=" << std::isnan(age)
                    << ", isNaN(lower)=" << std::isnan(age_lower)
                    << ", isNaN(upper)=" << std::isnan(age_upper)
                    << ", age <= 0: " << (age <= 0)
                    << "\n";
                ++debug_print_count;
            }

            if (missing_age || missing_lower || missing_upper ||
            std::isnan(age) || std::isnan(age_lower) || std::isnan(age_upper) || age <= 0) {
            ++total_missing;
            if (missing_age) std::cout << "  [DEBUG] Missing age_flame at index " << i << "\n";
            if (missing_lower) std::cout << "  [DEBUG] Missing age_flame_lower at index " << i << "\n";
            if (missing_upper) std::cout << "  [DEBUG] Missing age_flame_upper at index " << i << "\n";
            if (!missing_age && std::isnan(age)) std::cout << "  [DEBUG] age_flame is NaN at index " << i << "\n";
            if (!missing_lower && std::isnan(age_lower)) std::cout << "  [DEBUG] age_flame_lower is NaN at index " << i << "\n";
            if (!missing_upper && std::isnan(age_upper)) std::cout << "  [DEBUG] age_flame_upper is NaN at index " << i << "\n";
            if (!missing_age && age <= 0) std::cout << "  [DEBUG] age_flame <= 0 at index " << i << "\n";
            continue;
            }

            double age_err = 0.5 * (std::abs(age - age_lower) + std::abs(age_upper - age));
            double frac_err = age_err / age;

            if (frac_err > 0.3) {
            ++total_failed_uncertainty;
            std::cout << "  [DEBUG] frac_err > 0.3 at index " << i << " (frac_err=" << frac_err << ")\n";
            continue;
            }

            ++total_passed;
        }

        // If require_errorUVW_lt2 is set, check error_U/V/W < 2 km/s
        if (require_errorUVW_lt2 && use_errors) {
            ++total_checked_err;
            double eU = (*error_U_ptr)[i];
            double eV = (*error_V_ptr)[i];
            double eW = (*error_W_ptr)[i];
            if (std::isnan(eU) || std::isnan(eV) || std::isnan(eW) ||
                eU >= 2.0 || eV >= 2.0 || eW >= 2.0) {
                ++total_failed_err;
                continue;
            }
            ++total_passed_err;
        }

        //vtan cut
        if (require_vtan_ltval) {
            ++total_checked_vtan;
            if (i >= gaiaData.parallax.size() || i >= gaiaData.parallax_error.size()) {
                ++total_failed_vtan;
                continue;
            }
            double parallax = gaiaData.parallax[i];
            double parallax_error = gaiaData.parallax_error[i];
            double vtan = 4.74047 * std::sqrt(gaiaData.pmra[i] * gaiaData.pmra[i] + gaiaData.pmdec[i] * gaiaData.pmdec[i]) / parallax;
            
            // DEBUG PRINT
            if (i < 10) { // Print only first 10 for debugging
                std::cout << "DEBUG [" << i << "]: v_tan = " << vtan
                          << ", parallax = " << parallax
                          << ", parallax_error = " << parallax_error << "\n";
            }

            if (std::isnan(parallax) || std::isnan(parallax_error) || parallax_error == 0) {
                ++total_failed_vtan;
                continue;
            }

            if (vtan > 100.0) { // Arbitrary cut, can be adjusted
                ++total_failed_vtan;
                continue;
            }
            ++total_passed_vtan;
        }

        // NEW: total velocity filter
        if (require_total_velocity_lt150) {
            ++total_checked_totalvel;
            double Ulsr = gaiaData.U[i] + 11.10;
            double Vlsr = gaiaData.V[i] + 12.24 - 232.0;
            double Wlsr = gaiaData.W[i] + 7.25;
            double total_vel = std::sqrt(Ulsr * Ulsr + Vlsr * Vlsr + Wlsr * Wlsr);

            if (i < 10) {
                std::cout << "DEBUG [" << i << "]: Ulsr = " << Ulsr
                          << ", Vlsr = " << Vlsr
                          << ", Wlsr = " << Wlsr
                          << ", total_vel = " << total_vel << "\n";
            }

            if (std::isnan(total_vel) || total_vel > 150.0) {
                ++total_failed_totalvel;
                continue;
            }
            ++total_passed_totalvel;
        }

        int bin = static_cast<int>((val - param_start) / param_step);
        if (bin < 0 || bin >= static_cast<int>(num_bins)) continue;
        double eU = use_errors ? (*error_U_ptr)[i] : std::nan("");
        double eV = use_errors ? (*error_V_ptr)[i] : std::nan("");
        double eW = use_errors ? (*error_W_ptr)[i] : std::nan("");
        bin_map[bin].emplace_back(gaiaData.U[i], gaiaData.V[i], gaiaData.W[i], eU, eV, eW);
    }

    std::map<int, DispersionStats> dispersion_result;
    std::map<int, size_t> bin_counts;
    std::map<int, std::tuple<double, double, double, double>> uncertainty_map;
    std::map<int, std::tuple<double, double, double, double>> propagated_error_map;

    for (const auto& [bin, velocities] : bin_map) {
        double sum_U = 0, sum_V = 0, sum_W = 0;
        size_t count = velocities.size();
        for (const auto& [U, V, W, eU, eV, eW] : velocities) {
            sum_U += U; sum_V += V; sum_W += W;
        }
        double mean_U = sum_U / count, mean_V = sum_V / count, mean_W = sum_W / count;
        double sq_sum_U = 0, sq_sum_V = 0, sq_sum_W = 0;
        for (const auto& [U, V, W, eU, eV, eW] : velocities) {
            sq_sum_U += std::pow(U - mean_U, 2);
            sq_sum_V += std::pow(V - mean_V, 2);
            sq_sum_W += std::pow(W - mean_W, 2);
        }
        double sigma_U = std::sqrt(sq_sum_U / count);
        double sigma_V = std::sqrt(sq_sum_V / count);
        double sigma_W = std::sqrt(sq_sum_W / count);
        double sigma_total = std::sqrt(sigma_U * sigma_U + sigma_V * sigma_V + sigma_W * sigma_W);
        double delta_sigma_U = (count > 1) ? sigma_U / std::sqrt(2.0 * (count - 1)) : std::nan("");
        double delta_sigma_V = (count > 1) ? sigma_V / std::sqrt(2.0 * (count - 1)) : std::nan("");
        double delta_sigma_W = (count > 1) ? sigma_W / std::sqrt(2.0 * (count - 1)) : std::nan("");
        double delta_sigma_total = (count > 1) ? sigma_total / std::sqrt(2.0 * (count - 1)) : std::nan("");
        dispersion_result[bin] = {sigma_U, sigma_V, sigma_W, sigma_total};
        uncertainty_map[bin] = {delta_sigma_U, delta_sigma_V, delta_sigma_W, delta_sigma_total};
        bin_counts[bin] = count;

        // Propagate errors on UVW if available
        if (use_errors) {
            double sum_eU2 = 0, sum_eV2 = 0, sum_eW2 = 0;
            for (const auto& [U, V, W, eU, eV, eW] : velocities) {
                if (!std::isnan(eU)) sum_eU2 += eU * eU;
                if (!std::isnan(eV)) sum_eV2 += eV * eV;
                if (!std::isnan(eW)) sum_eW2 += eW * eW;
            }
            double mean_eU2 = sum_eU2 / count;
            double mean_eV2 = sum_eV2 / count;
            double mean_eW2 = sum_eW2 / count;
            double error_sigma_U = std::sqrt((sigma_U * sigma_U + mean_eU2) / (2.0 * (count - 1)));
            double error_sigma_V = std::sqrt((sigma_V * sigma_V + mean_eV2) / (2.0 * (count - 1)));
            double error_sigma_W = std::sqrt((sigma_W * sigma_W + mean_eW2) / (2.0 * (count - 1)));
            double error_sigma_total = std::sqrt((sigma_total * sigma_total +
                                                mean_eU2 + mean_eV2 + mean_eW2) / (2.0 * (count - 1)));
            propagated_error_map[bin] = {error_sigma_U, error_sigma_V, error_sigma_W, error_sigma_total};
        }
    }

    if (!output_filename.empty()) {
        std::filesystem::path out_path(output_filename);
        if (!std::filesystem::exists(out_path.parent_path())) {
            std::filesystem::create_directories(out_path.parent_path());
        }
        std::ofstream ofs(output_filename);
        if (!ofs) {
            std::cerr << "Error opening file for writing: " << output_filename << "\n";
        } else {
            ofs << "bin,param_lower,param_upper,param_center,"
                   "sigma_U,sigma_V,sigma_W,sigma_total,"
                   "delta_sigma_U,delta_sigma_V,delta_sigma_W,delta_sigma_total,"
                   "n_stars";
            if (use_errors) {
                ofs << ",error_sigma_U,error_sigma_V,error_sigma_W,error_sigma_total";
            }
            ofs << "\n";
            for (const auto& [bin, stats] : dispersion_result) {
                double param_lower = param_start + bin * param_step;
                double param_upper = param_lower + param_step;
                double param_center = param_lower + 0.5 * param_step;
                auto [delta_sigma_U, delta_sigma_V, delta_sigma_W, delta_sigma_total] = uncertainty_map.at(bin);
                size_t n_stars = bin_counts[bin];
                ofs << bin << "," << param_lower << "," << param_upper << "," << param_center << ","
                    << stats.sigma_U << "," << stats.sigma_V << "," << stats.sigma_W << "," << stats.sigma_total << ","
                    << delta_sigma_U << "," << delta_sigma_V << "," << delta_sigma_W << "," << delta_sigma_total << ","
                    << n_stars;
                if (use_errors) {
                    auto [errU, errV, errW, errTot] = propagated_error_map.at(bin);
                    ofs << "," << errU << "," << errV << "," << errW << "," << errTot;
                }
                ofs << "\n";
            }
            ofs.close();
        }
    }
    if (filter_teff_lum) {
        std::cout << "Teff/Luminosity filtering diagnostics:\n";
        std::cout << "  Checked stars: " << total_checked_tefflum << "\n";
        std::cout << "  Passed: " << total_passed_tefflum << "\n";
        std::cout << "  Failed (teff < " << teff_min << " or lum < " << lum_min << " or NaN): " << total_failed_tefflum << "\n";
    }
    if (filter_colour_mag) {
        std::cout << "Colour/Magnitude filtering diagnostics:\n";
        std::cout << "  Checked stars: " << total_checked_colmag << "\n";
        std::cout << "  Passed: " << total_passed_colmag << "\n";
        std::cout << "  Failed (bp_rp not in (" << bp_rp_min << "," << bp_rp_max << ") or M_G not in (" << M_G_min << "," << M_G_max << ") or NaN): " << total_failed_colmag << "\n";
    }
    if (require_parallax_over_error_gt10) {
        std::cout << "Parallax/Parallax_error filtering diagnostics:\n";
        std::cout << "  Checked stars: " << total_checked_parallax << "\n";
        std::cout << "  Passed: " << total_passed_parallax << "\n";
        std::cout << "  Failed parallax/parallax_error <= 10 or NaN: " << total_failed_parallax << "\n";
    }
    if (require_age_30pct) {
        std::cout << "Age filtering diagnostics:\n";
        std::cout << "  Checked stars: " << total_checked << "\n";
        std::cout << "  Passed: " << total_passed << "\n";
        std::cout << "  Failed uncertainty > 30%: " << total_failed_uncertainty << "\n";
        std::cout << "  Missing or invalid data: " << total_missing << "\n";
    }
    if (require_errorUVW_lt2 && use_errors) {
        std::cout << "UVW error filtering diagnostics:\n";
        std::cout << "  Checked stars: " << total_checked_err << "\n";
        std::cout << "  Passed: " << total_passed_err << "\n";
        std::cout << "  Failed error_U/V/W >= 2 km/s or NaN: " << total_failed_err << "\n";
    }
    if (require_vtan_ltval) {
        std::cout << "Vtan filtering diagnostics:\n";
        std::cout << "  Checked stars: " << total_checked_vtan << "\n";
        std::cout << "  Passed: " << total_passed_vtan << "\n";
        std::cout << "  Failed vtan > 100 km/s or NaN: " << total_failed_vtan << "\n";
    }
    if (require_total_velocity_lt150) {
        std::cout << "Total velocity filtering diagnostics:\n";
        std::cout << "  Checked stars: " << total_checked_totalvel << "\n";
        std::cout << "  Passed: " << total_passed_totalvel << "\n";
        std::cout << "  Failed total velocity > 150 km/s or NaN: " << total_failed_totalvel << "\n";
    }
    if (filtermaxage) {
        std::cout << "Max age filtering diagnostics:\n";
        std::cout << "  Checked stars: " << total_checked_maxage << "\n";
        std::cout << "  Passed: " << total_passed_maxage << "\n";
        std::cout << "  Failed (age >= " << max_age << " or NaN): " << total_failed_maxage << "\n";
    }
    return dispersion_result;
}
int main() {
    py::scoped_interpreter guard{};

    std::string folderpath = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets";        
    std::string input_csv = folderpath + "/GAIA10parallaxastrometry.csv";
        
    //GAIA gaiaInput;
    //gaiaInput.loadData(input_csv);

    //calculate UVW for every star in the file
    //std::vector<std::tuple<double, double, double>> uvw;
    //size_t n = gaiaInput.parallax.size();
    //for (size_t i = 0; i < n; ++i) {
    //    double parallax = gaiaInput.parallax[i];
    //    double radial_velocity = gaiaInput.radial_ve//locity[i];
    //    double pmra = gaiaInput.pmra[i];
    //    double pmdec = gaiaInput.pmdec[i];
    //    double ra = gaiaInput.ra[i];
    //    double dec = gaiaInput.dec[i];

            // Skip if any required value is NaN
    //    if (std::isnan(parallax) || std::isnan(radial_velocity) ||
    //        std::isnan(pmra) || std::isnan(pmdec)) {
    //        uvw.emplace_back(NAN, NAN, NAN);
    //        continue;
    //    }

    //    std::cout << "Calculating UVW for star " << i
    //              << " (ra=" << ra << ", dec=" << dec
    //              << ", parallax=" << parallax
    //              << ", pmra=" << pmra
    //              << ", pmdec=" << pmdec
    //              << ", radial_velocity=" << radial_velocity << ")"             
    //              << std::endl;
    //    auto [U, V, W] = compute3dspacevelocity(ra, dec, radial_velocity, pmra, pmdec, parallax);
    //    uvw.emplace_back(U, V, W);
    //}
    //saveCSVWithUVW(input_csv, output_csv, uvw);

    // Load both files
    // Load input files

    std::string uvw_csv = folderpath + "/GAIA10_UVW_appended.csv";
    std::string fallback_csv = folderpath + "/GAIA10_with_UVW.csv";
    std::string uvw_errors_csv = folderpath + "/GAIA10_with_UVW_errors.csv";
    GAIA gaiaUVW;
    gaiaUVW.loadData(uvw_csv);

    GAIA gaiaFallback;
    gaiaFallback.loadData(fallback_csv);

    GAIA gaiaErrors;
    gaiaErrors.loadData(uvw_errors_csv);

    // Build a set of source_ids already processed in appended CSV
    std::unordered_set<long long> processed_ids(gaiaUVW.source_id.begin(), gaiaUVW.source_id.end());

    // Combined vectors
    std::vector<double> combined_U, combined_V, combined_W, combined_age;
    std::vector<double> combined_error_U, combined_error_V, combined_error_W;
    std::vector<long long> combined_source_id;
    std::vector<double> combined_age_flame;
    std::vector<double> combined_age_flame_lower;
    std::vector<double> combined_age_flame_upper;
    std::vector<double> combined_parallax;
    std::vector<double> combined_parallax_error;
    std::vector<double> combined_pmra;
    std::vector<double> combined_pmdec;
    std::vector<double> combined_teff_gspphot;
    std::vector<double> combined_lum_flame;
    std::vector<double> combined_mh_gspspec;
    std::vector<double> combined_mh_gspspec_lower;
    std::vector<double> combined_mh_gspspec_upper;


    //mh vectors from GAIA10_with_UVW_errors.csv
    std::vector<double> mh_gspspec;
    std::vector<double> mh_gspspec_lower;
    std::vector<double> mh_gspspec_upper;
    std::vector<double> phot_g_mean_mag;
    std::vector<double> bp_rp;
    std::vector<double> U;
    std::vector<double> V;
    std::vector<double> W;
     
    if (gaiaErrors.mh_gspspec.size() > 0) {
        mh_gspspec = gaiaErrors.mh_gspspec;
        mh_gspspec_lower = gaiaErrors.mh_gspspec_lower;
        mh_gspspec_upper = gaiaErrors.mh_gspspec_upper;
    } else {
        std::cerr << "[WARNING] No metallicity data found in GAIA10_with_UVW_errors.csv" << std::endl;
    }   

    // Fallback age_flame + errors map
    std::unordered_map<long long, double> fallback_ageflame;
    std::unordered_map<long long, std::tuple<double, double, double>> fallback_errors;

    for (size_t i = 0; i < gaiaFallback.source_id.size(); ++i) {
        if (i < gaiaFallback.age_flame.size())
            fallback_ageflame[gaiaFallback.source_id[i]] = gaiaFallback.age_flame[i];
    }

    for (size_t i = 0; i < gaiaErrors.source_id.size(); ++i) {
        fallback_errors[gaiaErrors.source_id[i]] = std::make_tuple(
            gaiaErrors.error_U.size() > i ? gaiaErrors.error_U[i] : std::nan(""),
            gaiaErrors.error_V.size() > i ? gaiaErrors.error_V[i] : std::nan(""),
            gaiaErrors.error_W.size() > i ? gaiaErrors.error_W[i] : std::nan("")
        );
    }

    // Pass 1: From GAIA10_UVW_appended.csv
    for (size_t i = 0; i < gaiaUVW.source_id.size(); ++i) {
        double age = std::nan("");
        double age_flame = std::nan(""), age_lower = std::nan(""), age_upper = std::nan("");

        if (i < gaiaUVW.age_flame.size()) age_flame = gaiaUVW.age_flame[i];
        else std::cerr << "[WARNING] Missing age_flame at index " << i << std::endl;

        if (i < gaiaUVW.age_flame_lower.size()) age_lower = gaiaUVW.age_flame_lower[i];
        else std::cerr << "[WARNING] Missing age_flame_lower at index " << i << std::endl;

        if (i < gaiaUVW.age_flame_upper.size()) age_upper = gaiaUVW.age_flame_upper[i];
        else std::cerr << "[WARNING] Missing age_flame_upper at index " << i << std::endl;

        if (i < gaiaUVW.age_mean.size() && !std::isnan(gaiaUVW.age_mean[i])) {
            age = gaiaUVW.age_mean[i];
        } else {
            auto it = fallback_ageflame.find(gaiaUVW.source_id[i]);
            if (it != fallback_ageflame.end() && !std::isnan(it->second)) {
                age = it->second;
            }
        }

        auto err_it = fallback_errors.find(gaiaUVW.source_id[i]);
        double eU = std::nan(""), eV = std::nan(""), eW = std::nan("");
        if (err_it != fallback_errors.end()) {
            std::tie(eU, eV, eW) = err_it->second;
        }

        if (i < gaiaUVW.U.size() && i < gaiaUVW.V.size() && i < gaiaUVW.W.size() &&
            !std::isnan(age) &&
            !std::isnan(gaiaUVW.U[i]) &&
            !std::isnan(gaiaUVW.V[i]) &&
            !std::isnan(gaiaUVW.W[i])) {

            combined_U.push_back(gaiaUVW.U[i]);
            combined_V.push_back(gaiaUVW.V[i]);
            combined_W.push_back(gaiaUVW.W[i]);
            combined_age.push_back(age);
            combined_source_id.push_back(gaiaUVW.source_id[i]);
            combined_error_U.push_back(eU);
            combined_error_V.push_back(eV);
            combined_error_W.push_back(eW);
            combined_age_flame.push_back(age_flame);
            combined_age_flame_lower.push_back(age_lower);
            combined_age_flame_upper.push_back(age_upper);
            combined_parallax.push_back(gaiaUVW.parallax[i]);
            combined_parallax_error.push_back(gaiaUVW.parallax_error[i]);
            combined_pmra.push_back(gaiaUVW.pmra[i]);
            combined_pmdec.push_back(gaiaUVW.pmdec[i]);
            combined_teff_gspphot.push_back(gaiaUVW.teff_gspphot[i]);
            combined_lum_flame.push_back(gaiaUVW.lum_flame[i]);
            combined_mh_gspspec.push_back(gaiaUVW.mh_gspspec[i]);
            combined_mh_gspspec_lower.push_back(gaiaUVW.mh_gspspec_lower[i]);
            combined_mh_gspspec_upper.push_back(gaiaUVW.mh_gspspec_upper[i]);
        }
    }

    // Pass 2: From fallback file (those not in GAIA10_UVW_appended.csv)
    for (size_t i = 0; i < gaiaFallback.source_id.size(); ++i) {
        if (processed_ids.count(gaiaFallback.source_id[i]) == 0) {
            if (i < gaiaFallback.age_flame.size() &&
                i < gaiaFallback.U.size() &&
                i < gaiaFallback.V.size() &&
                i < gaiaFallback.W.size() &&
                !std::isnan(gaiaFallback.age_flame[i]) &&
                !std::isnan(gaiaFallback.U[i]) &&
                !std::isnan(gaiaFallback.V[i]) &&
                !std::isnan(gaiaFallback.W[i])) {

                combined_U.push_back(gaiaFallback.U[i]);
                combined_V.push_back(gaiaFallback.V[i]);
                combined_W.push_back(gaiaFallback.W[i]);
                combined_age.push_back(gaiaFallback.age_flame[i]);
                combined_source_id.push_back(gaiaFallback.source_id[i]);

                auto err_it = fallback_errors.find(gaiaFallback.source_id[i]);
                double eU = std::nan(""), eV = std::nan(""), eW = std::nan("");
                if (err_it != fallback_errors.end()) {
                    std::tie(eU, eV, eW) = err_it->second;
                }

                combined_error_U.push_back(eU);
                combined_error_V.push_back(eV);
                combined_error_W.push_back(eW);

                // Copy flame fields from fallback with checks
                double age_flame = std::nan(""), age_lower = std::nan(""), age_upper = std::nan("");
                if (i < gaiaFallback.age_flame.size()) age_flame = gaiaFallback.age_flame[i];
                else std::cerr << "[WARNING] Missing fallback age_flame at index " << i << std::endl;

                if (i < gaiaFallback.age_flame_lower.size()) age_lower = gaiaFallback.age_flame_lower[i];
                else std::cerr << "[WARNING] Missing fallback age_flame_lower at index " << i << std::endl;

                if (i < gaiaFallback.age_flame_upper.size()) age_upper = gaiaFallback.age_flame_upper[i];
                else std::cerr << "[WARNING] Missing fallback age_flame_upper at index " << i << std::endl;

                combined_age_flame.push_back(age_flame);
                combined_age_flame_lower.push_back(age_lower);
                combined_age_flame_upper.push_back(age_upper);
                combined_parallax.push_back(gaiaFallback.parallax[i]);
                combined_parallax_error.push_back(gaiaFallback.parallax_error[i]);
                combined_pmra.push_back(gaiaFallback.pmra[i]);
                combined_pmdec.push_back(gaiaFallback.pmdec[i]);
                combined_teff_gspphot.push_back(gaiaFallback.teff_gspphot[i]);
                combined_lum_flame.push_back(gaiaFallback.lum_flame[i]);
                combined_mh_gspspec.push_back(gaiaFallback.mh_gspspec[i]);
                combined_mh_gspspec_lower.push_back(gaiaFallback.mh_gspspec_lower[i]);
                combined_mh_gspspec_upper.push_back(gaiaFallback.mh_gspspec_upper[i]);
            }
        }
    }

    // Define binning
    double age_start = 0.1;
    double age_end = 14.0;
    size_t num_age_bins = 100;

    double mh_start = -1.5;
    double mh_end = 0.42;
    size_t num_mh_bins = 100;

    auto get_combined_age = [&combined_age](size_t i) { return combined_age[i]; };
    auto get_combined_mh = [&combined_mh_gspspec](size_t i) { return combined_mh_gspspec[i]; };

    auto get_mh = [&mh_gspspec](size_t i) {
        if (i < mh_gspspec.size()) {
            return mh_gspspec[i];
        }
        return std::nan("");
    };

    // Wrap results into GAIA object
    GAIA combinedGAIA;
    combinedGAIA.U = std::move(combined_U);
    combinedGAIA.V = std::move(combined_V);
    combinedGAIA.W = std::move(combined_W);
    combinedGAIA.age_flame = std::move(combined_age_flame);
    combinedGAIA.age_flame_lower = std::move(combined_age_flame_lower);
    combinedGAIA.age_flame_upper = std::move(combined_age_flame_upper);
    combinedGAIA.parallax = std::move(combined_parallax);
    combinedGAIA.parallax_error = std::move(combined_parallax_error);
    combinedGAIA.pmra = std::move(combined_pmra);
    combinedGAIA.pmdec = std::move(combined_pmdec);
    combinedGAIA.teff_gspphot = std::move(combined_teff_gspphot);
    combinedGAIA.lum_flame = std::move(combined_lum_flame);
    combinedGAIA.mh_gspspec = std::move(combined_mh_gspspec);
    combinedGAIA.mh_gspspec_lower = std::move(combined_mh_gspspec_lower);
    combinedGAIA.mh_gspspec_upper = std::move(combined_mh_gspspec_upper);

    // Compute dispersion with error propagation
    auto dispersion_result_1d = computeVelocityDispersionBinned1D(
        gaiaErrors,
        mh_start, mh_end, num_mh_bins,
        get_mh, //get_combined_age
        "feh",
        folderpath + "/outputdata/velocitydispersion/sigma_mh_binned_agecut_maxage" + std::to_string(num_age_bins) + ".csv",
        &combined_error_U,
        &combined_error_V,
        &combined_error_W,
        true,  // require_age_30pct
        false,  // require_errorUVW_lt2
        false,   // require_parallax_over_error_gt10
        true,   // require_vtan_ltval
        false,  // filter_teff_lum
        5000.0, // teff_min
        0.1,   // lum_min
        false,   // require_bp_rp
        0.3,   // bp_rp_min
        2.5,   // bp_rp_max
        3.0,  // M_G_min
        10.0,  // M_G_max
        8.0, //max age
        true, //filter total velocity
        true  //filter max age
    );

    std::cout << "Average no. of stars per age bin: " << combined_age.size() / num_age_bins << std::endl;

    std::string output_disp_csv = folderpath + "/outputdata/velocitydispersion/sigma_agemean_metmean_binned.csv";
    //std::string output_disp_csv = folderpath + "/outputdata/velocitydispersion/sigma_age_mh_binned_" +
    //                              std::to_string(num_age_bins) + "x" + std::to_string(num_feh_bins) + ".csv";
    //auto dispersion_result = computeVelocityDispersionBinned(
    //    gaiaUVW,
    //    age_start, age_end, num_age_bins,
    //    feh_start, feh_end, num_feh_bins,
    //    get_age, get_feh,
    //    output_disp_csv
    //);

    //std::cout << "Average no. of stars per age bin: " << gaiaUVW.age_flame.size() / num_age_bins << std::endl;
    //std::cout << "Average no. of stars per metallicity bin: " << gaiaUVW.mh_gspspec.size() / num_feh_bins << std::endl;
    //std::cout << "2D velocity dispersion by age and metallicity bin written to: " << output_disp_csv << std::endl;

    //auto dispersion_result_1d = computeVelocityDispersionBinned1D(
    //    gaiaUVW,
    //    age_start, age_end, num_age_bins,
    //    get_age,
    //    "age",
    //    folderpath + "/outputdata/velocitydispersion/combined_sigma_age_binned" + std::to_string(num_age_bins) + ".csv"
        //output_disp_csv
    //);
    //std::cout << "Average no. of stars per age bin: " << gaiaUVW.age_mean.size() / num_age_bins << std::endl;

    //std::string uvw_input_csv = folderpath + "/GAIA10_with_UVW.csv";
    //GAIA gaiaUVWErr;
    //gaiaUVWErr.loadData(uvw_input_csv);

    //std::vector<std::tuple<double, double, double>> uvw_errors;
    //size_t n = gaiaUVWErr.parallax.size();

    // Check that all required vectors are non-empty and have the same size
    //if (n == 0 ||
    //    gaiaUVWErr.ra.size() != n || gaiaUVWErr.dec.size() != n ||
    //    gaiaUVWErr.pmra.size() != n || gaiaUVWErr.pmdec.size() != n ||
    //    gaiaUVWErr.radial_velocity.size() != n || gaiaUVWErr.ra_error.size() != n ||
    //    gaiaUVWErr.dec_error.size() != n || gaiaUVWErr.parallax_error.size() != n ||
    //    gaiaUVWErr.pmra_error.size() != n || gaiaUVWErr.pmdec_error.size() != n ||
    //    gaiaUVWErr.radial_velocity_error.size() != n) {
    //    std::cerr << "Error: GAIA data vectors are empty or have inconsistent sizes. Aborting UVW error calculation." << std::endl;
    //   std::cerr << "Sizes: "
    //              << "ra=" << gaiaUVWErr.ra.size() << ", "
    //              << "dec=" << gaiaUVWErr.dec.size() << ", "
    //              << "pmra=" << gaiaUVWErr.pmra.size() << ", "
    //              << "pmdec=" << gaiaUVWErr.pmdec.size() << ", "
    //              << "radial_velocity=" << gaiaUVWErr.radial_velocity.size() << ", "
    //              << "ra_error=" << gaiaUVWErr.ra_error.size() << ", "
    //              << "dec_error=" << gaiaUVWErr.dec_error.size() << ", "
    //              << "parallax_error=" << gaiaUVWErr.parallax_error.size() << ", "
    //              << "pmra_error=" << gaiaUVWErr.pmra_error.size() << ", "
    //              << "pmdec_error=" << gaiaUVWErr.pmdec_error.size() << ", "
    //              << "radial_velocity_error=" << gaiaUVWErr.radial_velocity_error.size() << std::endl;
    //    return 1;
    //}

    //for (size_t i = 0; i < n; ++i) {
    //    double ra = gaiaUVWErr.ra[i];
    //    double dec = gaiaUVWErr.dec[i];
    //    double parallax = gaiaUVWErr.parallax[i];
    //    double pmra = gaiaUVWErr.pmra[i];
    //    double pmdec = gaiaUVWErr.pmdec[i];
    //    double radial_velocity = gaiaUVWErr.radial_velocity[i];
    //    double ra_error = gaiaUVWErr.ra_error[i];
    //    double dec_error = gaiaUVWErr.dec_error[i];
    //    double parallax_error = gaiaUVWErr.parallax_error[i];
    //    double pmra_error = gaiaUVWErr.pmra_error[i];
    //    double pmdec_error = gaiaUVWErr.pmdec_error[i];
    //    double radial_velocity_error = gaiaUVWErr.radial_velocity_error[i];

    //    auto [sigma_U, sigma_V, sigma_W] = compute3dspacevelocityerror(
    //        ra, dec, parallax, pmra, pmdec, radial_velocity,
    //        ra_error, dec_error, parallax_error, pmra_error, pmdec_error, radial_velocity_error
    //    );
    //    uvw_errors.emplace_back(sigma_U, sigma_V, sigma_W);
    //}

    // Append errors to the CSV
    //std::ifstream in(uvw_input_csv);
    //std::ofstream out(folderpath + "/GAIA10_with_UVW_errors.csv");
    //std::string line;
    //std::getline(in, line);
    //out << line << ",error_U,error_V,error_W\n";
    //size_t i = 0;
    //while (std::getline(in, line)) {
    //    out << line;
    //    if (i < uvw_errors.size()) {
    //        auto [sigma_U, sigma_V, sigma_W] = uvw_errors[i++];
    //        out << "," << sigma_U << "," << sigma_V << "," << sigma_W;
    //    }
    //    out << "\n";
    //}
    //in.close();
    //out.close();
    //std::cout << "UVW uncertainties appended to GAIA10_with_UVW_errors.csv\n";
    return 0;
}


