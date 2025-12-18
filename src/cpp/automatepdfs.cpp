#include <iostream>
#include <fstream>
#include <string>
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
#include "gaiaisogrid.h"
#include "gaialikelihood.h"
#include "isochronetypes.h"

// Add pybind11 namespace alias
namespace py = pybind11;

// Utility to read a 2D matrix from file
std::vector<std::vector<double>> load_posterior(const std::string& filename, size_t rows, size_t cols) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::vector<std::vector<double>> P(rows, std::vector<double>(cols));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (!(infile >> P[i][j])) {
                throw std::runtime_error("Error reading posterior matrix at " + std::to_string(i) + "," + std::to_string(j));
            }
        }
    }

    return P;

}
std::vector<double> gradient(const std::vector<double>& x) {
    size_t n = x.size();
    std::vector<double> grad(n, 0.0);
    if (n < 2) return grad;
    grad[0] = x[1] - x[0];
    for (size_t i = 1; i < n - 1; ++i) {
        grad[i] = 0.5 * (x[i + 1] - x[i - 1]);
    }
    grad[n - 1] = x[n - 1] - x[n - 2];
    return grad;
}

int countProcessedStars(const std::string& appended_csv_path) {
    std::ifstream appended_file(appended_csv_path);
    if (!appended_file.is_open()) return 0;

    std::string line;
    std::getline(appended_file, line); // skip header
    int count = 0;

    while (std::getline(appended_file, line)) {
        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string field;
        while (std::getline(ss, field, ',')) {
            fields.push_back(field);
        }

        // Check if age_mean (should be 5th from end) is non-empty
        if (fields.size() >= 4 && !fields[fields.size() - 4].empty()) {
            ++count;
        } else {
            break;  // stop at first unprocessed row
        }
    }

    return count;
}

int main() {
    py::scoped_interpreter guard{};
    auto start_time = std::chrono::high_resolution_clock::now();

    std::string folderPath = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets";
    std::string isochronefolderpath = folderPath + "/isochrones2";
    std::vector<Isochrone> isochrones = readIsochroneFilesIso(isochronefolderpath);
    size_t numIsochrones = isochrones.size();
    size_t numDataPoints = isochrones[0].initialMass.size();

    std::vector<std::vector<double>> weights(numIsochrones, std::vector<double>(numDataPoints, 0.0));
    ParameterSet paramSet = logTe_logLLo;

    // Manual entry option
    bool manualEntry = true;
    std::cout << "Do you want to manually enter star parameters? (y/n): ";
    std::string manualInput;
    std::getline(std::cin, manualInput);
    if (!manualInput.empty() && (manualInput[0] == 'y' || manualInput[0] == 'Y')) {
        manualEntry = true;
    }

    std::vector<std::tuple<double, double, double, double>> star_params;
    std::vector<std::string> lines;

    if (manualEntry) {
        double specifiedlogTe, specifiedlogLLo, specifiedGbpGrp, specifiedG;
        std::cout << "Enter logTeff: ";
        std::cin >> specifiedlogTe;
        std::cout << "Enter logLLo: ";
        std::cin >> specifiedlogLLo;
        std::cout << "Enter GbpGrp: ";
        std::cin >> specifiedGbpGrp;
        std::cout << "Enter G: ";
        std::cin >> specifiedG;
        star_params.emplace_back(specifiedlogTe, specifiedlogLLo, specifiedG, specifiedGbpGrp);
        lines.push_back("manual_entry");
    } else {
        // ... (original CSV reading and selection code here)
        // [Paste the original CSV reading and selection code from your previous main]
        // For brevity, not repeated here.
        // Make sure to keep the rest of your code unchanged.
        std::string gaia_csv = folderPath + "/GAIA10_with_UVW.csv";
        std::ifstream infile(gaia_csv);
        if (!infile) {
            std::cerr << "Could not open GAIA10_with_UVW.csv for reading.\n";
            return 1;
        }
        std::string header, line;
        std::getline(infile, header);
        std::vector<std::string> columns;
        std::stringstream header_ss(header);
        std::string col;
        while (std::getline(header_ss, col, ',')) {
            columns.push_back(col);
        }
        int idx_teff = -1, idx_lum = -1, idx_g = -1, idx_bp_rp = -1;
        for (size_t i = 0; i < columns.size(); ++i) {
            if (columns[i] == "teff_gspphot") idx_teff = i;
            if (columns[i] == "lum_flame") idx_lum = i;
            if (columns[i] == "phot_g_mean_mag") idx_g = i;
            if (columns[i] == "bp_rp") idx_bp_rp = i;
        }
        if (idx_teff == -1 || idx_lum == -1 || idx_g == -1 || idx_bp_rp == -1) {
            std::cerr << "Could not find required columns in GAIA10_with_UVW.csv\n";
            return 1;
        }
        std::vector<std::string> lines_csv;
        lines_csv.push_back(header);
        while (std::getline(infile, line)) {
            lines_csv.push_back(line);
            std::stringstream ss(line);
            std::vector<std::string> fields;
            std::string field;
            while (std::getline(ss, field, ',')) {
                fields.push_back(field);
            }
            if (fields.size() > std::max({idx_teff, idx_lum, idx_g, idx_bp_rp})) {
                try {
                    std::string teff_str = fields[idx_teff];
                    std::string lum_str  = fields[idx_lum];
                    std::string g_str    = fields[idx_g];
                    std::string bp_rp_str = fields[idx_bp_rp];

                    if (teff_str.empty() || lum_str.empty() || g_str.empty() || bp_rp_str.empty()) {
                        std::cerr << "Skipping line number: " << (lines_csv.size()) << "/" << (star_params.size() + 1) << " with missing data.\n";
                        continue;
                    }

                    double specifiedlogTe = std::log10(std::stod(teff_str));
                    double specifiedlogLLo = std::log10(std::stod(lum_str));
                    double specifiedG = std::stod(g_str);
                    double specifiedGbpGrp = std::stod(bp_rp_str);

                    star_params.emplace_back(specifiedlogTe, specifiedlogLLo, specifiedG, specifiedGbpGrp);
                } catch (const std::exception& e) {
                    std::cerr << "Error parsing numeric fields:\n";
                    std::cerr << "Exception: " << e.what() << std::endl;
                    std::abort();
                }
            }
        }
        infile.close();
        lines = lines_csv;
        // ... (continue with binning and selection as in your original code)
    }

    // Setup output file and write header
    std::string gaia_appended_csv = folderPath + "/GAIA10_UVW_appended.csv";
    std::string new_header = manualEntry ? "manual_entry,age_mean,age_std,met_mean,met_std" : lines[0] + ",age_mean,age_std,met_mean,met_std";
    bool fileExists = std::filesystem::exists(gaia_appended_csv);
    std::ofstream outfile(gaia_appended_csv, std::ios::app);
    if (!outfile) {
        std::cerr << "Failed to open " << gaia_appended_csv << " for appending.\n";
        return 1;
    }
    if (!fileExists) {
        outfile << new_header << "\n";
        outfile.flush();
    }

    int startStar = 0;
    int endStar = static_cast<int>(star_params.size()) - 1;
    size_t starsCount = endStar - startStar + 1;

    std::ostringstream param1Stream, param2Stream;
    if (paramSet == logTe_logLLo) {
        param1Stream << "dummy1";
        param2Stream << "dummy2";
    }

    std::string weightsFilename = isochronefolderpath + "/weights/weights_" + (paramSet == logTe_logLLo ? "logTe_logLLo" : "Gbp_G") + ".csv";
    loadWeights(weightsFilename, isochrones, weights, paramSet);

    std::string prior_file = folderPath + "/outputdata/prior_values.txt";
    std::vector<std::vector<double>> prior;

    GAIA gaiaData;
    gaiaData.loadData(folderPath + "/GAIA10_with_UVW.csv");

    for (size_t idx = startStar; idx <= static_cast<size_t>(endStar); ++idx) {
        double specifiedlogTe, specifiedlogLLo, specifiedG, specifiedGbpGrp;
        std::tie(specifiedlogTe, specifiedlogLLo, specifiedG, specifiedGbpGrp) = star_params[idx];

        std::ostringstream param1StreamStar, param2StreamStar;
        if (paramSet == logTe_logLLo) {
            param1StreamStar << std::fixed << std::setprecision(2) << specifiedlogTe;
            param2StreamStar << std::fixed << std::setprecision(2) << specifiedlogLLo;
        }
        std::string probabilitiesFilename = folderPath + "/probabilities/probabilities_" 
            + (paramSet == logTe_logLLo ? "logTe" : "G_bp") 
            + param1StreamStar.str() + "_" 
            + (paramSet == logTe_logLLo ? "logLLo" : "G") 
            + param2StreamStar.str() + ".txt";

        std::vector<std::vector<double>> jointPDF = binningisogrid(isochrones, weights, specifiedlogTe, specifiedlogLLo, specifiedGbpGrp, specifiedG, probabilitiesFilename, paramSet);

        if (prior.empty()) {
            prior = load_posterior(prior_file, jointPDF.size(), jointPDF[0].size());
        }

        std::vector<size_t> filteredIndices = gaiaData.filterData(specifiedlogTe, specifiedlogLLo, false);
        std::vector<size_t> validFilteredIndices = matchStarsWithIsochrones(gaiaData, isochrones, filteredIndices);

        std::vector<double> totalLikelihoods;
        calculateResidualsAndLikelihoods(gaiaData, isochrones, validFilteredIndices, totalLikelihoods);

        BinningResult binningResult = binning(gaiaData, validFilteredIndices, totalLikelihoods, specifiedlogTe, specifiedlogLLo, true);
        std::vector<std::vector<double>> smoothedCppArray = binningResult.smoothedLikelihoods;

        size_t n_age = jointPDF.size();
        size_t n_met = jointPDF[0].size();
        std::vector<std::vector<double>> posterior(n_age, std::vector<double>(n_met, 0.0));
        for (size_t i = 0; i < n_age; ++i)
            for (size_t j = 0; j < n_met; ++j)
                posterior[i][j] = prior[i][j] * jointPDF[i][j] * smoothedCppArray[i][j];

        std::ostringstream posterior_filename;
        if (manualEntry) {
            posterior_filename << folderPath << "/outputdata/posteriors/posterior_manual_logTe"
                << std::fixed << std::setprecision(2) << specifiedlogTe
                << "_logLLo" << std::fixed << std::setprecision(2) << specifiedlogLLo
                << "_GbpGrp" << std::fixed << std::setprecision(2) << specifiedGbpGrp
                << "_G" << std::fixed << std::setprecision(2) << specifiedG << ".txt";
        } else {
            posterior_filename << folderPath << "/outputdata/posteriors/posterior_logTe"
                << std::fixed << std::setprecision(2) << specifiedlogTe
                << "_logLLo" << std::fixed << std::setprecision(2) << specifiedlogLLo
                << "_GbpGrp" << std::fixed << std::setprecision(2) << specifiedGbpGrp
                << "_G" << std::fixed << std::setprecision(2) << specifiedG << ".txt";
        }

        std::ofstream posterior_file(posterior_filename.str());
        if (posterior_file) {
            for (size_t i = 0; i < n_age; ++i) {
                for (size_t j = 0; j < n_met; ++j) {
                    posterior_file << posterior[i][j];
                    if (j + 1 < n_met) posterior_file << " ";
                }
                posterior_file << "\n";
            }
            posterior_file.close();
        }

        std::vector<double> age(n_age), metallicity(n_met);
        double age_min = 0.1, age_max = 14.0;
        double met_min = -1.5, met_max = 0.42;
        for (size_t i = 0; i < n_age; ++i)
            age[i] = age_min + i * (age_max - age_min) / (n_age - 1);
        for (size_t j = 0; j < n_met; ++j)
            metallicity[j] = met_min + j * (met_max - met_min) / (n_met - 1);

        std::vector<double> dt = gradient(age);
        std::vector<double> dZ = gradient(metallicity);

        std::vector<double> posterior_age(n_age, 0.0), posterior_met(n_met, 0.0);
        for (size_t i = 0; i < n_age; ++i)
            for (size_t j = 0; j < n_met; ++j)
                posterior_age[i] += posterior[i][j] * dZ[j];

        for (size_t j = 0; j < n_met; ++j)
            for (size_t i = 0; i < n_age; ++i)
                posterior_met[j] += posterior[i][j] * dt[i];

        double area_age = 0.0, area_met = 0.0;
        for (size_t i = 0; i < n_age; ++i) area_age += posterior_age[i] * dt[i];
        for (size_t j = 0; j < n_met; ++j) area_met += posterior_met[j] * dZ[j];
        for (double& p : posterior_age) p /= area_age;
        for (double& p : posterior_met) p /= area_met;

        double age_mean = 0.0, age_var = 0.0;
        for (size_t i = 0; i < n_age; ++i) {
            age_mean += age[i] * posterior_age[i] * dt[i];
            age_var += age[i] * age[i] * posterior_age[i] * dt[i];
        }
        double age_std = std::sqrt(age_var - age_mean * age_mean);

        double met_mean = 0.0, met_var = 0.0;
        for (size_t j = 0; j < n_met; ++j) {
            met_mean += metallicity[j] * posterior_met[j] * dZ[j];
            met_var += metallicity[j] * metallicity[j] * posterior_met[j] * dZ[j];
        }
        double met_std = std::sqrt(met_var - met_mean * met_mean);

        std::ostringstream outrow;
        if (manualEntry) {
            outrow << "manual_entry";
        } else {
            outrow << lines[idx + 1];
        }
        if (std::isnan(age_mean) || std::isnan(age_std) || std::isnan(met_mean) || std::isnan(met_std)) {
            outrow << ",nan,nan,nan,nan";
        } else {
            outrow << "," << age_mean << "," << age_std << "," << met_mean << "," << met_std;
        }
        outfile << outrow.str() << "\n";

        std::cout << "Processed star " << (idx - startStar + 1) << " / " << starsCount << "\n";
        std::cout << "Star parameters: logTe = " << specifiedlogTe
                  << ", logLLo = " << specifiedlogLLo
                  << ", GbpGrp = " << specifiedGbpGrp
                  << ", G = " << specifiedG << "\n";
        std::cout << "Posterior Age mean: " << age_mean << ", std: " << age_std
                  << ", Met mean: " << met_mean << ", std: " << met_std << "\n";
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    double total_seconds = elapsed.count();
    double avg_seconds = (starsCount > 0) ? total_seconds / starsCount : 0.0;

    std::cout << "\nTotal elapsed time: " << total_seconds << " seconds\n";
    std::cout << "Average time per star: " << avg_seconds << " seconds\n";
    return 0;
}
