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


namespace fs = std::filesystem;
namespace plt = matplotlibcpp;
namespace py = pybind11;

class Isochrone { //stores header values and arrays for data points
public:
    double age; //located in header
    double MH;
    double Z;
    std::vector<double> isoweight; //located in table
    std::vector<double> deltaM;
    std::vector<double> deltaMH;
    std::vector<double> deltaT;
    std::vector<double> initialMass;
    std::vector<double> mass;
    std::vector<double> logLLo;
    std::vector<double> logTe;
    std::vector<double> G;
    std::vector<double> Gbp;
    std::vector<double> Grp;
    Isochrone() : age(0.0), MH(0.0), Z(0.0) {} //default constructor
    Isochrone (double age, double MH, double Z)
        : age(age), MH(MH), Z(Z), isoweight(0.0), deltaM(0.0), deltaMH(0.0), deltaT(0.0) {}

    void addData(double initialMass, double mass, double logLLo, double logTe, double G, double Gbp, double Grp) {
        this->initialMass.push_back(initialMass);
        this->mass.push_back(mass);
        this->logLLo.push_back(logLLo);
        this->logTe.push_back(logTe);
        this->G.push_back(G);
        this->Gbp.push_back(Gbp);
        this->Grp.push_back(Grp);
        this->isoweight.push_back(0.0); //initialises isoweight etc to 0.0 by adding 0.0 to the vector
        this->deltaM.push_back(0.0);
        this->deltaMH.push_back(0.0);
        this->deltaT.push_back(0.0);
    } //adds data points to the arrays
    //this-> accesses members of current object
    //push_back adds the value to the end of the vector
    //e.g. this->initialMass.push_back(initialMass) adds the value of initialMass to the end of the initialMass vector of the current isochrone object

};

enum ParameterSet {
    logTe_logLLo,
    Gbp_G
};

//struct to store file metadata for sorting
struct IsochroneFileMetadata {
    std::string filePath;
    double age; //age in Gyr
    double Z;   //metallicity
};

//helper function to parse filename for age and Z
bool parseIsochroneFilename(const std::string& filename, double& age, double& Z) {
    // Example: 14000z0363071y298P00O0D0E0.isc_gaia-dr3
    std::regex pattern(R"((\d+)z(\d+)y\d+P\d+O\d+D\d+E\d+\.isc_gaia-dr3)");
    std::smatch matches;

    if (std::regex_search(filename, matches, pattern) && matches.size() == 3) {
        age = std::stod(matches[1].str()) / 1000.0; // Age in Gyr (e.g., 14000 → 14.0)
        Z = std::stod("0." + matches[2].str());     // Z as decimal (e.g., 0363071 → 0.0363071)
        return true;
    }
    return false;
}

//modified function to read isochrones in order
std::vector<Isochrone> readIsochroneFiles(const std::string& folderPath) {
    std::vector<Isochrone> isochrones; // Stores all isochrones
    std::vector<IsochroneFileMetadata> fileMetadataList; // Temporary container for sorting

    //step 1:iterate through files and extract metadata
    for (const auto& entry : fs::directory_iterator(folderPath)) {
        std::cout << "Found file: " << entry.path() << std::endl;
        if (entry.path().extension() == ".isc_gaia-dr3") { // Check for the correct file extension
            double age, Z;
            std::string filename = entry.path().filename().string();
            std::cout << "Checking filename: " << filename << std::endl;

            if (parseIsochroneFilename(filename, age, Z)) {
                std::cout << "Parsed age: " << age << ", Z: " << Z << std::endl;
                fileMetadataList.push_back({entry.path().string(), age, Z});
            } else {
                std::cout << "Filename did not match pattern: " << filename << std::endl;
            }
        }
    }

    //step 2:sort files by age and Z
    std::sort(fileMetadataList.begin(), fileMetadataList.end(), [](const IsochroneFileMetadata& a, const IsochroneFileMetadata& b) {
        if (a.age == b.age) {
            return a.Z < b.Z; //sort by Z if ages are equal
        }
        return a.age < b.age; //otherwise, sort by age
    });

    //step 3:read and process isochrones in sorted order
    for (const auto& metadata : fileMetadataList) {
        std::ifstream file(metadata.filePath);
        if (file.is_open()) {
            Isochrone isochrone(0.0, 0.0, 0.0); //initialize isochrone object
            std::string line;
            bool ageFound = false, MHFound = false, ZFound = false;

            //read header lines to extract age, MH, Z
            for (int j = 0; j < 8; ++j) {
                std::getline(file, line);
    
                std::smatch match;

                // Match Age (Myr)
                if (!ageFound && std::regex_search(line, match, std::regex(R"(Age \(Myr\) = ([\d\.]+))"))) {
                    isochrone.age = std::stod(match[1].str()) / 1000.0; // Convert to Gyr
                    ageFound = true;
            }

                // Match [M/H]
                if (!MHFound && std::regex_search(line, match, std::regex(R"(\[M/H\] = ([\d\.\-]+))"))) {
                    isochrone.MH = std::stod(match[1].str());
                    MHFound = true;
            }

                // Match Z
                if (!ZFound && std::regex_search(line, match, std::regex(R"(Z = ([\d\.]+))"))) {
                    isochrone.Z = std::stod(match[1].str());
                    ZFound = true;
            }
}

            // Read data points
            while (std::getline(file, line)) {
                std::istringstream iss(line);
                double initialMass, mass, logLLo, logTe, G, Gbp, Grp;
                if (iss >> initialMass >> mass >> logLLo >> logTe >> G >> Gbp >> Grp) {
                    isochrone.addData(initialMass, mass, logLLo, logTe, G, Gbp, Grp);
                } else {
                    std::cerr << "Invalid data format in file: " << metadata.filePath << std::endl;
                }
            }
            file.close();
            isochrones.push_back(isochrone);
        } else {
            std::cerr << "Could not open file: " << metadata.filePath << std::endl;
        }
    }

    return isochrones; //return the vector of isochrones
}

class GAIA {
public:
    std::vector<double> parallax;
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

    GAIA() = default;
    GAIA (std::vector<double> parallax, std::vector<double> distance_gspphot, std::vector<double> distance_gspphot_lower, std::vector<double> distance_gspphot_upper, std::vector<double> teff_gspphot, std::vector<double>  teff_gspphot_lower, std::vector<double> teff_gspphot_upper, std::vector<double> logg_gspphot, std::vector<double> logg_gspphot_lower, std::vector<double> logg_gspphot_upper, std::vector<double> bp_rp, std::vector<double> mh_gspphot, std::vector<double> mh_gspphot_lower, std::vector<double> mh_gspphot_upper, std::vector<double> phot_g_mean_mag, std::vector<long long int> source_id, std::vector<double> age_flame, std::vector<double> age_flame_lower, std::vector<double> age_flame_upper, std::vector<double> mass_flame, std::vector<double> mass_flame_lower, std::vector<double> mass_flame_upper, std::vector<double> lum_flame, std::vector<double> lum_flame_lower, std::vector<double> lum_flame_upper, std::vector<double> mh_gspspec, std::vector<double> mh_gspspec_lower, std::vector<double> mh_gspspec_upper)
        : parallax(parallax), distance_gspphot(distance_gspphot), distance_gspphot_lower(distance_gspphot_lower), distance_gspphot_upper(distance_gspphot_upper), teff_gspphot(teff_gspphot), teff_gspphot_lower(teff_gspphot_lower), teff_gspphot_upper(teff_gspphot_upper), logg_gspphot(logg_gspphot), logg_gspphot_lower(logg_gspphot_lower), logg_gspphot_upper(logg_gspphot_upper), bp_rp(bp_rp), mh_gspphot(mh_gspphot), mh_gspphot_lower(mh_gspphot_lower), mh_gspphot_upper(mh_gspphot_upper), phot_g_mean_mag(phot_g_mean_mag), source_id(source_id), age_flame(age_flame), age_flame_lower(age_flame_lower), age_flame_upper(age_flame_upper), mass_flame(mass_flame), mass_flame_lower(mass_flame_lower), mass_flame_upper(mass_flame_upper), lum_flame(lum_flame), lum_flame_lower(lum_flame_lower), lum_flame_upper(lum_flame_upper), mh_gspspec(mh_gspspec), mh_gspspec_lower(mh_gspspec_lower), mh_gspspec_upper(mh_gspspec_upper) {}

    void addData(double parallax, double distance_gspphot, double distance_gspphot_lower, double distance_gspphot_upper, double teff_gspphot, double teff_gspphot_lower, double teff_gspphot_upper, double logg_gspphot, double logg_gspphot_lower, double logg_gspphot_upper, double bp_rp, double mh_gspphot, double mh_gspphot_lower, double mh_gspphot_upper, double phot_g_mean_mag, long long int source_id, double age_flame, double age_flame_lower, double age_flame_upper, double mass_flame, double mass_flame_lower, double mass_flame_upper, double lum_flame, double lum_flame_lower, double lum_flame_upper, double mh_gspspec, double mh_gspspec_lower, double mh_gspspec_upper)
     {
        this->parallax.push_back(parallax);
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
    }
    void loadData(const std::string& filename);
    std::vector<size_t> filterData(double logTe, double logL, bool useAllStars);
    std::vector<std::pair<double, double>> extractAgeAndMetallicity(const std::vector<size_t>& indices) const;
    
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
    // Debug: Print the column indices
    //std::cout << "Column Indices:" << std::endl;
    //for (const auto& pair : columnIndices) {
    //    std::cout << pair.first << ": " << pair.second << std::endl;
    //}

    // Read the data lines
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<std::string> values;

        // Split the line into values
        while (std::getline(ss, value, ',')) {
            values.push_back(value);
        }
        // Debug: Print the values
        // std::cout << "Values:" << std::endl;
        // for (const auto& val : values) {
        //    std::cout << val << " ";
        //}
        //std::cout << std::endl;

        // Check if the number of values matches the number of columns
        if (values.size() != columnIndices.size()) {
            while (values.size() < columnIndices.size()) {
                values.push_back(""); // Add empty strings for missing values
            }
        }


        try {
            // Assign values to the appropriate vectors
            //std::cerr << "Parsing parallax: '" << values[columnIndices["parallax"]] << "'" << std::endl;
            parallax.push_back(values[columnIndices["parallax"]].empty() ? std::nan("") : std::stod(values[columnIndices["parallax"]]));

            //std::cerr << "Parsing distance_gspphot: '" << values[columnIndices["distance_gspphot"]] << "'" << std::endl;
            distance_gspphot.push_back(values[columnIndices["distance_gspphot"]].empty() ? std::nan("") : std::stod(values[columnIndices["distance_gspphot"]]));

            //std::cerr << "Parsing distance_gspphot_lower: '" << values[columnIndices["distance_gspphot_lower"]] << "'" << std::endl;
            distance_gspphot_lower.push_back(values[columnIndices["distance_gspphot_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["distance_gspphot_lower"]]));

            //std::cerr << "Parsing distance_gspphot_upper: '" << values[columnIndices["distance_gspphot_upper"]] << "'" << std::endl;
            distance_gspphot_upper.push_back(values[columnIndices["distance_gspphot_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["distance_gspphot_upper"]]));

            //std::cerr << "Parsing teff_gspphot: '" << values[columnIndices["teff_gspphot"]] << "'" << std::endl;
            teff_gspphot.push_back(values[columnIndices["teff_gspphot"]].empty() ? std::nan("") : std::stod(values[columnIndices["teff_gspphot"]]));

            //std::cerr << "Parsing teff_gspphot_lower: '" << values[columnIndices["teff_gspphot_lower"]] << "'" << std::endl;
            teff_gspphot_lower.push_back(values[columnIndices["teff_gspphot_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["teff_gspphot_lower"]]));

            //std::cerr << "Parsing teff_gspphot_upper: '" << values[columnIndices["teff_gspphot_upper"]] << "'" << std::endl;
            teff_gspphot_upper.push_back(values[columnIndices["teff_gspphot_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["teff_gspphot_upper"]]));

            //std::cerr << "Parsing logg_gspphot: '" << values[columnIndices["logg_gspphot"]] << "'" << std::endl;
            logg_gspphot.push_back(values[columnIndices["logg_gspphot"]].empty() ? std::nan("") : std::stod(values[columnIndices["logg_gspphot"]]));

            //std::cerr << "Parsing logg_gspphot_lower: '" << values[columnIndices["logg_gspphot_lower"]] << "'" << std::endl;
            logg_gspphot_lower.push_back(values[columnIndices["logg_gspphot_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["logg_gspphot_lower"]]));

            //std::cerr << "Parsing logg_gspphot_upper: '" << values[columnIndices["logg_gspphot_upper"]] << "'" << std::endl;
            logg_gspphot_upper.push_back(values[columnIndices["logg_gspphot_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["logg_gspphot_upper"]]));

            //std::cerr << "Parsing bp_rp: '" << values[columnIndices["bp_rp"]] << "'" << std::endl;
            bp_rp.push_back(values[columnIndices["bp_rp"]].empty() ? std::nan("") : std::stod(values[columnIndices["bp_rp"]]));

            //std::cerr << "Parsing mh_gspphot: '" << values[columnIndices["mh_gspphot"]] << "'" << std::endl;
            mh_gspphot.push_back(values[columnIndices["mh_gspphot"]].empty() ? std::nan("") : std::stod(values[columnIndices["mh_gspphot"]]));

            //std::cerr << "Parsing mh_gspphot_lower: '" << values[columnIndices["mh_gspphot_lower"]] << "'" << std::endl;
            mh_gspphot_lower.push_back(values[columnIndices["mh_gspphot_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["mh_gspphot_lower"]]));

            //std::cerr << "Parsing mh_gspphot_upper: '" << values[columnIndices["mh_gspphot_upper"]] << "'" << std::endl;
            mh_gspphot_upper.push_back(values[columnIndices["mh_gspphot_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["mh_gspphot_upper"]]));

            //std::cerr << "Parsing phot_g_mean_mag: '" << values[columnIndices["phot_g_mean_mag"]] << "'" << std::endl;
            phot_g_mean_mag.push_back(values[columnIndices["phot_g_mean_mag"]].empty() ? std::nan("") : std::stod(values[columnIndices["phot_g_mean_mag"]]));

            //std::cerr << "Parsing source_id: '" << values[columnIndices["source_id"]] << "'" << std::endl;
            source_id.push_back(values[columnIndices["source_id"]].empty() ? -1 : std::stoll(values[columnIndices["source_id"]]));

            //std::cerr << "Parsing age_flame: '" << values[columnIndices["age_flame"]] << "'" << std::endl;
            age_flame.push_back(values[columnIndices["age_flame"]].empty() ? std::nan("") : std::stod(values[columnIndices["age_flame"]]));

            //std::cerr << "Parsing age_flame_lower: '" << values[columnIndices["age_flame_lower"]] << "'" << std::endl;
            age_flame_lower.push_back(values[columnIndices["age_flame_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["age_flame_lower"]]));

            //std::cerr << "Parsing age_flame_upper: '" << values[columnIndices["age_flame_upper"]] << "'" << std::endl;
            age_flame_upper.push_back(values[columnIndices["age_flame_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["age_flame_upper"]]));

            //std::cerr << "Parsing mass_flame: '" << values[columnIndices["mass_flame"]] << "'" << std::endl;
            mass_flame.push_back(values[columnIndices["mass_flame"]].empty() ? std::nan("") : std::stod(values[columnIndices["mass_flame"]]));

            //std::cerr << "Parsing mass_flame_lower: '" << values[columnIndices["mass_flame_lower"]] << "'" << std::endl;
            mass_flame_lower.push_back(values[columnIndices["mass_flame_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["mass_flame_lower"]]));

            //std::cerr << "Parsing mass_flame_upper: '" << values[columnIndices["mass_flame_upper"]] << "'" << std::endl;
            mass_flame_upper.push_back(values[columnIndices["mass_flame_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["mass_flame_upper"]]));

            //std::cerr << "Parsing lum_flame: '" << values[columnIndices["lum_flame"]] << "'" << std::endl;
            lum_flame.push_back(values[columnIndices["lum_flame"]].empty() ? std::nan("") : std::stod(values[columnIndices["lum_flame"]]));

            //std::cerr << "Parsing lum_flame_lower: '" << values[columnIndices["lum_flame_lower"]] << "'" << std::endl;
            lum_flame_lower.push_back(values[columnIndices["lum_flame_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["lum_flame_lower"]]));

            //std::cerr << "Parsing lum_flame_upper: '" << values[columnIndices["lum_flame_upper"]] << "'" << std::endl;
            lum_flame_upper.push_back(values[columnIndices["lum_flame_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["lum_flame_upper"]]));

            //std::cerr << "Parsing mh_gspspec: '" << values[columnIndices["mh_gspspec"]] << "'" << std::endl;
            mh_gspspec.push_back(values[columnIndices["mh_gspspec"]].empty() ? std::nan("") : std::stod(values[columnIndices["mh_gspspec"]]));

            //std::cerr << "Parsing mh_gspspec_lower: '" << values[columnIndices["mh_gspspec_lower"]] << "'" << std::endl;
            mh_gspspec_lower.push_back(values[columnIndices["mh_gspspec_lower"]].empty() ? std::nan("") : std::stod(values[columnIndices["mh_gspspec_lower"]]));

            //std::cerr << "Parsing mh_gspspec_upper: '" << values[columnIndices["mh_gspspec_upper"]] << "'" << std::endl;
            mh_gspspec_upper.push_back(values[columnIndices["mh_gspspec_upper"]].empty() ? std::nan("") : std::stod(values[columnIndices["mh_gspspec_upper"]]));
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument: " << e.what() << " in line: " << line << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range: " << e.what() << " in line: " << line << std::endl;
        }

    }

    file.close();
}

//function to print isochrones
void print(const std::vector<Isochrone>& isochrones) {
    for (const auto& isochrone : isochrones) {
        //std::cout << "Isochrone: age=" << isochrone.age << ", MH=" << isochrone.MH << ", Z=" << isochrone.Z << std::endl;
    }
}

//function to print GAIA data
void print(const GAIA& gaiaData) {
    size_t limit = std::min(gaiaData.parallax.size(), static_cast<size_t>(10));
    for (size_t i = 0; i < limit; ++i) {
        std::cout << "GAIA Data: parallax=" << gaiaData.parallax[i]
                  << ", distance_gspphot=" << gaiaData.distance_gspphot[i]
                  << ", teff_gspphot=" << gaiaData.teff_gspphot[i]
                  << ", logg_gspphot=" << gaiaData.logg_gspphot[i]
                  << ", bp_rp=" << gaiaData.bp_rp[i]
                  << ", mh_gspphot=" << gaiaData.mh_gspphot[i]
                  << ", phot_g_mean_mag=" << gaiaData.phot_g_mean_mag[i]
                  << ", source_id=" << gaiaData.source_id[i]
                  << ", age_flame=" << gaiaData.age_flame[i]
                  << ", mass_flame=" << gaiaData.mass_flame[i]
                  << ", lum_flame=" << gaiaData.lum_flame[i]
                  << ", mh_gspspec=" << gaiaData.mh_gspspec[i]
                  << std::endl;
    }
}

std::vector<size_t> GAIA::filterData(double logTe, double logL, bool useAllStars) {
    std::vector<size_t> filteredIndices;

    for (size_t i = 0; i < teff_gspphot.size(); ++i) {
        if (teff_gspphot[i] > 0 && lum_flame[i] > 0) { // ensure valid values

            if (useAllStars) {
                // Include all valid stars, skip filtering
                filteredIndices.push_back(i);
                continue;
            }

            // Filtering mode
            double logTeff = std::log10(teff_gspphot[i]);
            double logLum = std::log10(lum_flame[i]);
            if (std::abs(logTeff - logTe) < 0.2 && std::abs(logLum - logL) < 0.1) {
                filteredIndices.push_back(i);
            }
        }
    }

    //std::cout << "Filtered Indices:" << std::endl;
    //for (size_t index : filteredIndices) {
    //    std::cout << "Index: " << index 
    //              << ", Teff: " << teff_gspphot[index] 
    //              << ", Lum: " << lum_flame[index] 
    //              << ", Age: " << age_flame[index]
    //              << ", Metallicity: " << mh_gspspec[index]
    //              << std::endl;
    //}
    //std::cout << std::endl;

    return filteredIndices;
}

std::vector<std::pair<double, double>> GAIA::extractAgeAndMetallicity(const std::vector<size_t>& indices) const {
    std::vector<std::pair<double, double>> ageAndMetallicity;
    for (size_t index : indices) {
        double age = age_flame[index];
        double metallicity = mh_gspspec[index];
        ageAndMetallicity.push_back(std::make_pair(age, metallicity));
    }
    return ageAndMetallicity;
}

Isochrone findClosestIsochrone(const std::vector<Isochrone>& isochrones, double age, double metallicity) {
    Isochrone closestIsochrone;
    double minDistance = std::numeric_limits<double>::max();
    for (const auto& isochrone : isochrones) {
        double distance = std::sqrt(std::pow(isochrone.age - age, 2) + std::pow(isochrone.MH - metallicity, 2));
        if (distance < minDistance) {
            minDistance = distance;
            closestIsochrone = isochrone;
        }
    }
    return closestIsochrone;
}

size_t findClosestPointIndex(const Isochrone& isochrone, double logTe, double logL) {
    size_t closestIndex = 0;
    double minDistance = std::numeric_limits<double>::max();
    for (size_t j = 0; j < isochrone.logTe.size(); ++j) {
        double distance = std::sqrt(std::pow(isochrone.logTe[j] - logTe, 2) + std::pow(isochrone.logLLo[j] - logL, 2));
        if (distance < minDistance) {
            minDistance = distance;
            closestIndex = j;
        }
    }
    return closestIndex;
}

std::vector<size_t> matchStarsWithIsochrones(const GAIA& gaiaData, const std::vector<Isochrone>& isochrones, const std::vector<size_t>& filteredIndices) {
    size_t totalStars = filteredIndices.size(); //total stars within range
    size_t validStars = 0; //counter for stars with no NaN values and valid isochrone data
    std::vector<size_t> validFilteredIndices; //stores valid filtered indices
    std::vector<std::pair<double, double>> ageAndMetallicity = gaiaData.extractAgeAndMetallicity(filteredIndices);
    for (size_t i = 0; i < filteredIndices.size(); ++i) {
        double age = ageAndMetallicity[i].first;
        double metallicity = ageAndMetallicity[i].second;
        //check for valid metallicity
        if (std::isnan(metallicity)) {
            //std::cerr << "Error: Metallicity is NaN for star at index " << filteredIndices[i] << std::endl;
            continue;
        }
        //check for valid age
        if (std::isnan(age)) {
            //std::cerr << "Error: Age is NaN for star at index " << filteredIndices[i] << std::endl;
            continue;
        }
        Isochrone closestIsochrone = findClosestIsochrone(isochrones, age, metallicity);
        //check if the closest isochrone has valid data
        if (closestIsochrone.logTe.empty() || closestIsochrone.logLLo.empty()) {
            //std::cerr << "Error: Closest isochrone has no data for star at index " << filteredIndices[i] << std::endl;
            continue;
        }
        double logTe = std::log10(gaiaData.teff_gspphot[filteredIndices[i]]);
        double logL = std::log10(gaiaData.lum_flame[filteredIndices[i]]);
        size_t closestPointIndex = findClosestPointIndex(closestIsochrone, logTe, logL);
        //check if the closestPointIndex is valid
        if (closestPointIndex >= closestIsochrone.logTe.size() || closestPointIndex >= closestIsochrone.logLLo.size()) {
            std::cerr << "Error: Closest point index is out of bounds for star at index " << filteredIndices[i] << std::endl;
            continue;
        }
        //increment valid stars counter
        validStars++;
        validFilteredIndices.push_back(filteredIndices[i]);
        //std::cout << "Star with Age: " << age << ", Metallicity: " << metallicity
        //          << ", logTe: " << logTe << ", logL: " << logL
        //          << " matched with Isochrone point with Age: " << closestIsochrone.age << ", Metallicity: " << closestIsochrone.MH << ", logTe: " << closestIsochrone.logTe[closestPointIndex]
        //          << ", logL: " << closestIsochrone.logLLo[closestPointIndex] << std::endl;
    }
    //print summary
    //std::cout << "Total stars within range: " << totalStars << std::endl;
    //std::cout << "Total valid stars matched: " << validStars << std::endl;
    //std::cout << "After matchStarsWithIsochrones:" << std::endl;
    //std::cout << "Filtered Indices Size: " << filteredIndices.size() << std::endl;
    //std::cout << "Valid Filtered Indices Size: " << validFilteredIndices.size() << std::endl;
    //std::cout << "Isochrones Size: " << isochrones.size() << std::endl;
    //std::cout << "GAIA Data Size: " << gaiaData.teff_gspphot.size() << std::endl;

    return validFilteredIndices;
}

double calculateResidualLikelihood(double observed, double expected, double sigma) {
    if (sigma <= 0) {
        std::cerr << "Error: Standard deviation (sigma) must be positive." << std::endl;
        return 0.0;
    }
    double residual = observed - expected;
    return std::exp(-0.5 * std::pow(residual / sigma, 2.0)) / (std::sqrt(2.0 * M_PI) * sigma);
}
double convertlogg(double mass, double logL, double logTe) {
    const double logT_sun = std::log10(5777.0);
    const double logg_sun = 4.438;
    //calculate logg using mass, logL, and logTe
    double logg = logg_sun + std::log10(mass) - logL + 4.0 * logTe - 4.0
     * logT_sun;
    return logg;
}

double convertdistancetoMG(double parallax, double phot_g_mean_mag) {
    double MG = phot_g_mean_mag - 5 * std::log10(1000/parallax) + 5;
    return MG;
}

double convertMGtoMV(double MG, double bp_rp) {
    double MV = MG -0.017 -( 0.006 * bp_rp) - (0.176 * bp_rp * bp_rp);
    return MV;
}

double convertcolours(double VI) {
    double bp_rp = -0.0162 + (1.274 * VI ) - (0.08143 * VI * VI);
    return bp_rp;
}

double defineerrors(double upper, double lower, double threshold = 0.1) {
    if (std::isnan(upper) || std::isnan(lower)) {
        std::cerr << "Error: Upper or lower error value is NaN. Assigning default sigma." << std::endl;
        return threshold; //default sigma value (adjust)
    }
    double sigma = (upper - lower) / 2.0;
    if (sigma <= 0) {
        std::cerr << "Error: Calculated sigma is non-positive. Assigning default sigma." << std::endl;
        return threshold; //default sigma value
    }
    if (sigma < threshold) {
        sigma = std::sqrt(sigma * sigma + threshold * threshold);
    }
    return sigma;
}
//add 0.1 in quadrature for small values of sigma

void calculateResidualsAndLikelihoods(const GAIA& gaiaData, const std::vector<Isochrone>& isochrones, const std::vector<size_t>& filteredIndices, std::vector<double>& totalLikelihoods) {
    size_t validStars = 0; //counter for stars with valid data
    //std::cout << "Entering calculateResidualsAndLikelihoods..." << std::endl;
    //safety checks on input sizes
    if (filteredIndices.empty() || isochrones.empty() || gaiaData.teff_gspphot.empty()) {
        //std::cerr << "Error: One or more required datasets are empty. Exiting function." << std::endl;
        return;
    }

    //std::cout << "Filtered Indices Size: " << filteredIndices.size() << std::endl;

    if (isochrones.empty()) {
        //std::cerr << "Error: Isochrones vector is empty." << std::endl;
        return;
    }
    //std::cout << "Isochrones Size: " << isochrones.size() << std::endl;
    if (!isochrones.empty()) {
        //std::cout << "First Isochrone Age: " << isochrones[0].age << ", MH: " << isochrones[0].MH << std::endl;
    }

    //std::cout << "GAIA Data Size: " << gaiaData.teff_gspphot.size() << std::endl;
    if (!gaiaData.teff_gspphot.empty()) {
        //std::cout << "First GAIA Teff: " << gaiaData.teff_gspphot[0] << std::endl;
    }

    std::vector<std::pair<double, double>> ageAndMetallicity = gaiaData.extractAgeAndMetallicity(filteredIndices);
    //std::cout << "Age and Metallicity Size: " << ageAndMetallicity.size() << std::endl;
    if (ageAndMetallicity.size() != filteredIndices.size()) {
        //std::cerr << "Error: Mismatch between filteredIndices size (" << filteredIndices.size()
        //      << ") and ageAndMetallicity size (" << ageAndMetallicity.size() << ")." << std::endl;
    return;
    }

    totalLikelihoods.clear();

    for (size_t i = 0; i < filteredIndices.size(); ++i) {
        //std::cout << "Processing star at filtered index: " << filteredIndices[i] << std::endl;
        //std::cout << "Age: " << ageAndMetallicity[i].first << ", Metallicity: " << ageAndMetallicity[i].second << std::endl;
        double totalLikelihoodForStar = 0.0; //initialize total likelihood for the star
        //ensure index is within bounds
        if (filteredIndices[i] >= gaiaData.teff_gspphot.size() || 
            filteredIndices[i] >= gaiaData.logg_gspphot.size() || 
            filteredIndices[i] >= gaiaData.lum_flame.size()) {
            //std::cerr << "Error: Index " << filteredIndices[i] << " is out of bounds for GAIA data arrays." << std::endl;
            continue;
        }

        double age = ageAndMetallicity[i].first;
        double metallicity = ageAndMetallicity[i].second;

        //validate GAIA data field
        if (std::isnan(gaiaData.logg_gspphot[filteredIndices[i]]) ||
            std::isnan(gaiaData.teff_gspphot[filteredIndices[i]]) ||
            std::isnan(gaiaData.lum_flame[filteredIndices[i]])) {
            //std::cerr << "Error: GAIA data contains NaN for star at index " << filteredIndices[i] << std::endl;
            continue;
        }

        //std::cout << "Finding closest isochrone for Age: " << age << ", Metallicity: " << metallicity << std::endl;
        Isochrone closestIsochrone = findClosestIsochrone(isochrones, age, metallicity);
        //std::cout << "Closest Isochrone found: Age=" << closestIsochrone.age << ", MH=" << closestIsochrone.MH << std::endl;    

        //check if the closest isochrone has valid data
        if (closestIsochrone.logTe.empty() || closestIsochrone.logLLo.empty()) {
            //std::cerr << "Error: Closest isochrone has no data for star at index " << filteredIndices[i] << std::endl;
            continue;
        }

        double logTe = std::log10(gaiaData.teff_gspphot[filteredIndices[i]]);
        double logL = std::log10(gaiaData.lum_flame[filteredIndices[i]]);
        double observedlogg = gaiaData.logg_gspphot[filteredIndices[i]];
        double totalLikelihoodForPoint = 0.0;
        //increment valid stars counter
        validStars++;

        //calculate residuals and likelihoods for each parameter over each isochrone
        for (size_t j = 0; j < closestIsochrone.logTe.size(); ++j){
            //isochrone data for point j
            double isochroneLogTe = closestIsochrone.logTe[j];
            double isochroneLogL = closestIsochrone.logLLo[j];
            double isochroneLogg = convertlogg(
                closestIsochrone.mass[j],
                closestIsochrone.logLLo[j],
                closestIsochrone.logTe[j]
            );
            //calculate residuals and likelihoods for each parameter at each point
            // 1. logg
            double loggLikelihood = calculateResidualLikelihood(
                observedlogg,
                isochroneLogg,
                defineerrors(
                    gaiaData.logg_gspphot_upper[filteredIndices[i]],
                    gaiaData.logg_gspphot_lower[filteredIndices[i]]
                )
            );
            // 2. magnitude (MG and MV)
            double gaiaMG = convertdistancetoMG(
                gaiaData.parallax[filteredIndices[i]],
                gaiaData.phot_g_mean_mag[filteredIndices[i]]
            );

            // double gaiaMV = convertMGtoMV(gaiaMG, gaiaData.bp_rp[filteredIndices[i]]); we have gaiaMG in the isochrones in this case
            double mvLikelihood = calculateResidualLikelihood(
                gaiaMG,
                closestIsochrone.G[j],
                0.1 //example uncertainty
                );
            // 3. colour (BP-RP)
            double gaiaBP_RP = gaiaData.bp_rp[filteredIndices[i]];
            double isochroneBP_RP = closestIsochrone.Gbp[j] - closestIsochrone.Grp[j];
            double colorLikelihood = calculateResidualLikelihood(
                gaiaBP_RP,
                isochroneBP_RP,
                0.05 //example uncertainty for BP-RP
            );
            // 4. logTe
            double logTeLikelihood = calculateResidualLikelihood(
                logTe,
                closestIsochrone.logTe[j],
                defineerrors(
                    gaiaData.teff_gspphot_upper[filteredIndices[i]],
                    gaiaData.teff_gspphot_lower[filteredIndices[i]]
                )
            );
            // 5. logL
            double logLLikelihood = calculateResidualLikelihood(
                logL,
                closestIsochrone.logLLo[j],
                defineerrors(
                    gaiaData.lum_flame_upper[filteredIndices[i]],
                    gaiaData.lum_flame_lower[filteredIndices[i]]
                )
            );
            //combine likelihoods
            double totalLikelihoodForPoint = loggLikelihood * mvLikelihood * colorLikelihood * logTeLikelihood * logLLikelihood;
            //add to total likelihoods, accumulate likelihood
            totalLikelihoodForStar += totalLikelihoodForPoint;
        }
        //store total likelihood for the star
        totalLikelihoods.push_back(totalLikelihoodForStar);
        //print summary
        //std::cout << "Star at index " << filteredIndices[i]
        //          << " with Age: " << age
        //          << ", Metallicity: " << metallicity
        //          << ", logTe: " << logTe
        //          << ", logL: " << logL
        //          << ", matched with isochrone with age: " << closestIsochrone.age
        //          << ", MH: " << closestIsochrone.MH
        //          << ", Total Likelihood: " << totalLikelihoodForStar
        //          << "\n" 
        //          << std::endl;
    }

    //print summary
    //std::cout << "Total valid stars with residuals and likelihoods calculated: " << validStars << std::endl;
}

void saveSmoothedProbabilitiesToText(const std::vector<std::vector<double>>& finalArray, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        //std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    for (const auto& row : finalArray) {
        for (const auto& value : row) {
            outFile << value << " ";
        }
        outFile << "\n";
    }

    outFile.close();
}

struct BinningResult {
    std::vector<std::vector<double>> binnedLikelihoods;
    std::vector<std::vector<double>> smoothedLikelihoods;
};

BinningResult binning(const GAIA& gaiaData, const std::vector<size_t>& filteredIndices, std::vector<double>& totalLikelihoods, double logTe, double logL, bool Logscale = false) {
    //py::scoped_interpreter guard{}; //start the interpreter and keep it alive
    double numAgeBins = 50;
    double numMHBins = 48;
    double minAge = 0.1;
    double maxAge = 14.0;
    double minMH = -1.5;
    double maxMH = 0.42;
    double ageBinWidth = (maxAge - minAge) / numAgeBins;
    double MHBinWidth = (maxMH - minMH) / numMHBins;

    //initialize 2D grid for binned likelihoods and star counts
    std::vector<std::vector<double>> binnedLikelihoods(numAgeBins, std::vector<double>(numMHBins, 0.0));
    std::vector<std::vector<double>> averageLikelihood(numAgeBins, std::vector<double>(numMHBins, 0.0));
    std::vector<std::vector<int>> starCounts(numAgeBins, std::vector<int>(numMHBins, 0));
    std::vector<std::vector<int>> totalStarCounts(numAgeBins, std::vector<int>(numMHBins, 0));

    // Iterate over all stars in the dataset to calculate total star counts per bin
    size_t totalStarsProcessed = 0; 
    for (size_t i = 0; i < gaiaData.age_flame.size(); ++i) {
        double age = gaiaData.age_flame[i];
        double metallicity = gaiaData.mh_gspspec[i];

    // Skip invalid ages or metallicities
        if (std::isnan(age) || std::isnan(metallicity)) {
            continue;
        }

    // Calculate bin indices
        int ageBinIndex = static_cast<int>((age - minAge) / ageBinWidth);
        int MHBinIndex = static_cast<int>((metallicity - minMH) / MHBinWidth);

    // Ensure indices are within bounds
        if (ageBinIndex >= 0 && ageBinIndex < numAgeBins && MHBinIndex >= 0 && MHBinIndex   < numMHBins) {
        totalStarCounts[ageBinIndex][MHBinIndex]++; // Increment total star count for the bin
        totalStarsProcessed++; // Increment total star count
        } else {
        // Debug: Print stars that fall outside the binning range
            //std::cerr << "Star with Age: " << age << ", Metallicity: " << metallicity
            //      << " falls outside the binning range." << std::endl;
        }
    }

    //iterate over filtered indices
    for (size_t i = 0; i < filteredIndices.size(); ++i) {
        double age = gaiaData.age_flame[filteredIndices[i]];
        double metallicity = gaiaData.mh_gspspec[filteredIndices[i]];
        double likelihood = totalLikelihoods[i];

        //skip invalid ages or metallicities
        if (std::isnan(age) || std::isnan(metallicity)) {
            //std::cerr << "Skipping star with invalid age or metallicity at index " << filteredIndices[i] << std::endl;
            continue;
        }

        //std::cout << "Star at index " << filteredIndices[i]
        //          << " with Age: " << age
        //          << ", Metallicity: " << metallicity
        //          << ", Likelihood: " << likelihood << std::endl;

        //calculate bin indices
        int ageBinIndex = static_cast<int>((age - minAge) / ageBinWidth);
        int MHBinIndex = static_cast<int>((metallicity - minMH) / MHBinWidth);
        

        //ensure indices are within bounds
        if (ageBinIndex < 0 || ageBinIndex >= numAgeBins || MHBinIndex < 0 || MHBinIndex >= numMHBins) {
            //std::cerr << "Skipping star with out-of-bounds age or metallicity at index " << filteredIndices[i] << std::endl;
            continue;
        }

        //accumulate likelihood into the corresponding bin (idk if this is correct method)
        binnedLikelihoods[ageBinIndex][MHBinIndex] += likelihood;
        starCounts[ageBinIndex][MHBinIndex]++; //increment star count for the bin
    }
    //calculate average likelihood for each bin
    //for (int ageBin = 0; ageBin < numAgeBins; ++ageBin) {
    //    for (int MHBin = 0; MHBin < numMHBins; ++MHBin) {
    //        if (starCounts[ageBin][MHBin] > 0) {
    //            binnedLikelihoods[ageBin][MHBin] /= starCounts[ageBin][MHBin]; //average likelihood
    //        }
    //    }
    //}
    //normalize the likelihoods so that the sum of all likelihoods equals 1
    double totalLikelihoodSum = 0.0;
    //calculate the total sum of all likelihoods
    for (const auto& row : binnedLikelihoods) {
        for (double value : row) {
            totalLikelihoodSum += value;
        }
    }
    //normalize each bin's likelihood
    if (totalLikelihoodSum > 0.0) {
        for (auto& row : binnedLikelihoods) {
            for (double& value : row) {
                value /= totalLikelihoodSum; // Normalize to make the sum of all bins equal to 1
            }
        }
    }
    //debug, printing total star counts per bin
    //std::cout << "\nTotal star counts per bin:\n" << std::endl;
    for (int ageBin = 0; ageBin < numAgeBins; ++ageBin) {
        for (int MHBin = 0; MHBin < numMHBins; ++MHBin) {
            if (totalStarCounts[ageBin][MHBin] > 0) {
                double ageBinStart = minAge + ageBin * ageBinWidth;
                double ageBinEnd = ageBinStart + ageBinWidth;
                double MHBinStart = minMH + MHBin * MHBinWidth;
                double MHBinEnd = MHBinStart + MHBinWidth;
                //std::cout << "Bin (Age: [" << ageBinStart << ", " 
                //          << ageBinEnd   << "], " <<  "Metallicity: [" 
                //          << MHBinStart << ", " << MHBinEnd << "]) "
                //          << "Total stars: " << totalStarCounts[ageBin][MHBin] << ", "
                //          << "Likelihood: " << binnedLikelihoods[ageBin][MHBin]
                //          << std::endl;
            }
        }
    }
    //std::cout << "\nTotal stars in dataset with ages:" << totalStarsProcessed;
    //convert totalstarcounts to a format suitable for plotting
    std::vector<std::vector<double>> starCountsForPlot(numAgeBins, std::vector<double>(numMHBins, 0.0));
    for (int ageBin = 0; ageBin < numAgeBins; ++ageBin) {
        for (int MHBin = 0; MHBin < numMHBins; ++MHBin) {
            starCountsForPlot[ageBin][MHBin] = static_cast<double>(totalStarCounts[ageBin][MHBin]);
        }
    }
    //import matplotlib and pyplot and colors
    py::module_ plt = py::module_::import("matplotlib.pyplot");
    py::module_ font_manager = py::module_::import("matplotlib.font_manager");
    py::module_ colors = py::module_::import("matplotlib.colors");
    
    //plot the total star counts using matplotlib
    //plt.attr("figure")();
    //find the smallest nonzero value in starCountsForPlot
    if (Logscale) { // for logscale we want the 0 values to be filled with a small value
        double minNonZero = std::numeric_limits<double>::max();
        for (const auto& row : starCountsForPlot) {
            for (double val : row) {
                if (val > 0.0 && val < minNonZero) {
                    minNonZero = val;
                }
            }
        }
        //set zeros to a dex lower than the smallest nonzero value
        double fillValue = (minNonZero < std::numeric_limits<double>::max()) ? minNonZero / 10.0 : 1e-2;
        for (auto& row : starCountsForPlot) {
            for (auto& val : row) {
                if (val == 0.0) val = fillValue;
            }
        }
    }
    else {//do nothing 
    }

    py::array_t<double> numpyArrayStarCounts = py::cast(starCountsForPlot);
    std::string outputDir = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/plots/";
    //plt.attr("imshow")(numpyArrayStarCounts, 
                       //py::arg("extent") = py::make_tuple(minMH, maxMH, minAge, maxAge), 
                       //py::arg("origin") = "lower", 
                       //py::arg("aspect") = "auto", 
                       //py::arg("cmap") = "plasma"; // Use a different color map for distinction
    //plt.attr("xlabel")("Metallicity [MH]");
    //plt.attr("ylabel")("Age [Gyr]");
    //plt.attr("title")("Total Star Counts per Bin");
    //plt.attr("colorbar")();
    std::ostringstream starCountFilenameStream;
    starCountFilenameStream << outputDir << "likelihoods/total_star_counts_GAIA10.png";
    std::string starCountFilename = starCountFilenameStream.str();
    //plt.attr("savefig")(starCountFilename);
    //plt.attr("close")();

    //data for plotting
    std::vector<double> x, y, z; // x = mh, y = AGE ,z = likelihood
    for (int ageBin = 0; ageBin < numAgeBins; ++ageBin) {
        for (int MHBin = 0; MHBin < numMHBins; ++MHBin) {
            double ageCenter = minAge + (ageBin + 0.5) * ageBinWidth;
            double MHCenter = minMH + (MHBin + 0.5) * MHBinWidth;
            x.push_back(MHCenter);
            y.push_back(ageCenter);
            z.push_back(binnedLikelihoods[ageBin][MHBin]);
        }
    }
    if (Logscale) { // for logscale we want the 0 values to be filled with a small value
        double minNonZero = std::numeric_limits<double>::max();
        for (const auto& row : binnedLikelihoods) {
            for (double val : row) {
                if (val > 0.0 && val < minNonZero) {
                    minNonZero = val;
                }
            }
        }
        //set zeros to a dex lower than the smallest nonzero value
        double fillValue = (minNonZero < std::numeric_limits<double>::max()) ? minNonZero / 10.0 : 1e-2;
        for (auto& row : binnedLikelihoods) {
            for (auto& val : row) {
                if (val == 0.0) val = fillValue;
            }
        }
    }
    else {//do nothing 
    }
    //py::array_t<double> numpyArray = py::cast(binnedLikelihoods);
    //add kernel smoothing to the binnedLikelihoods
    // Step 1: Convert to log (add small epsilon to avoid log(0))
    py::module_ np = py::module_::import("numpy");
    py::module_ scipy_ndimage = py::module_::import("scipy.ndimage");

    double epsilon = 1e-40;
    py::array_t<double> numpyArray = py::cast(binnedLikelihoods);
    py::object safeArray = np.attr("maximum")(numpyArray, epsilon);
    py::object logArray = np.attr("log")(safeArray);

    // Step 2: Smooth in log space
    double sigma_age = 1.0; //3.6 based on bin widths
    double sigma_mh = 1.0; //2.5
    py::object logSmoothed = scipy_ndimage.attr("gaussian_filter")(
        logArray, py::make_tuple(sigma_age, sigma_mh)
    );

    // Step 3: Exponentiate
    py::object smoothedArray = np.attr("exp")(logSmoothed);

    // Step 4: Renormalize
    double total = py::float_(np.attr("sum")(smoothedArray));
    if (total > 0.0) {
        smoothedArray = smoothedArray.attr("__truediv__")(total);
    }

    // Optional: Cast back to py::array_t<double> if needed
    py::array_t<double> finalArray = smoothedArray.cast<py::array_t<double>>();

    // Step 5: Convert finalArray to std::vector<std::vector<double>> for saving
    auto buf = finalArray.unchecked<2>();
    ssize_t numRows = buf.shape(0);
    ssize_t numCols = buf.shape(1);

    std::vector<std::vector<double>> smoothedCppArray(numRows, std::vector<double>(numCols));

    for (ssize_t i = 0; i < numRows; ++i) {
        for (ssize_t j = 0; j < numCols; ++j) {
            smoothedCppArray[i][j] = buf(i, j);
        }
    }


    //plot the data using matplotlib-cpp
    if (Logscale) {
        //plt.attr("imshow")(finalArray, py::arg("extent") = py::make_tuple(minMH, maxMH, minAge, maxAge), py::arg("origin") = "lower", py::arg("aspect") = "auto", py::arg("norm") = colors.attr("LogNorm")());
    } else {
        //plt.attr("imshow")(finalArray, py::arg("extent") = py::make_tuple(minMH, maxMH, minAge, maxAge), py::arg("origin") = "lower", py::arg("aspect") = "auto");
    }
    //plt.attr("xlabel")("Metallicity [MH]");
    //plt.attr("ylabel")("Age [Gyr]");
    std::ostringstream title;
    title << "Age-Metallicity Observational Likelihoods for stars at \n LogTe = " << logTe
      << " and LogL/Lo = " << logL;
    //plt.attr("title")(title.str());
    //plt.attr("colorbar")();
    std::filesystem::create_directory(outputDir);
    std::ostringstream filenameStream;
    filenameStream << outputDir << "gaia/likelihoods/likelihood_"
               << (Logscale ? "log_" : "linear_")
               << "logL" << logL << "_logTe" << logTe << ".png";
    std::string filename = filenameStream.str();
    //plt.attr("savefig")(filename);
    //plt::show();

    // Debug: Count the number of stars and calculate the average likelihood for each bin
    //std::cout << "\nBin Summary (Age, Metallicity):\n";
    for (int ageBin = 0; ageBin < numAgeBins; ++ageBin) {
        for (int MHBin = 0; MHBin < numMHBins; ++MHBin) {
            if (binnedLikelihoods[ageBin][MHBin] > 0.0) { // Ignore bins with 0 likelihood
                double ageBinStart = minAge + ageBin * ageBinWidth;
                double ageBinEnd = ageBinStart + ageBinWidth;
                double MHBinStart = minMH + MHBin * MHBinWidth;
                double MHBinEnd = MHBinStart + MHBinWidth;
            // Count the number of stars contributing to this bin
                int starCount = 0;
                for (size_t i = 0; i < filteredIndices.size(); ++i) {
                    double age = gaiaData.age_flame[filteredIndices[i]];
                    double metallicity = gaiaData.mh_gspspec[filteredIndices[i]];

                    int starAgeBinIndex = static_cast<int>((age - minAge) / ageBinWidth);
                    int starMHBinIndex = static_cast<int>((metallicity - minMH) / MHBinWidth);

                    if (starAgeBinIndex == ageBin && starMHBinIndex == MHBin) {
                    ++starCount;
                    }
                }

                // Calculate the average likelihood for the bin
                double averageLikelihood = binnedLikelihoods[ageBin][MHBin] / starCount;

                // Print the bin details
                //std::cout << "Bin (Age: [" << ageBinStart << ", " 
                //          << ageBinEnd   << "], " <<  "Metallicity: [" 
                //          << MHBinStart << ", " << MHBinEnd << "]) "
                //          << "Stars: " << starCount
                //          << ", Average Likelihood: " << averageLikelihood << std::endl;
            }
        }
    }
    return {binnedLikelihoods, smoothedCppArray};
}

void saveProbabilitiesToText(const std::vector<std::vector<double>>& averageLikelihood, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        //std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    for (const auto& row : averageLikelihood) {
        for (const auto& value : row) {
            outFile << value << " ";
        }
        outFile << "\n";
    }
    outFile.close();
}
#ifdef TEST_GAIALIKELIHOOD_MAIN
int main() {
    //py::scoped_interpreter guard{}; // Initialize Python interpreter for pybind11
    //REMOVE^^IF DOING AUTOMATEPDFS
    std::string folderpath = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets";
    GAIA gaiaData;
    gaiaData.loadData(folderpath + "/GAIA10parallax.csv");
    std::vector<Isochrone> isochrones = readIsochroneFiles(folderpath + "/isochrones2");
    //print(isochrones);
    //print(gaiaData);
    double logTe = 3.81;
    double logL = 0.64; //example values of star wanting to be assessed
    std::vector<size_t> filteredIndices = gaiaData.filterData(logTe, logL, false); //false for filtering mode
    std::vector<size_t> validFilteredIndices = matchStarsWithIsochrones(gaiaData, isochrones, filteredIndices);
    std::vector<double> totalLikelihoods;  
    
    calculateResidualsAndLikelihoods(gaiaData, isochrones, validFilteredIndices, totalLikelihoods);
    BinningResult binningResult = binning(gaiaData, validFilteredIndices, totalLikelihoods, logTe, logL, true);
    std::vector<std::vector<double>> averageLikelihood = binningResult.binnedLikelihoods;
    std::vector<std::vector<double>> smoothedCppArray = binningResult.smoothedLikelihoods;
    //metallicityfit(gaiaData);
    std::ostringstream probabilityFilenameStream;
    probabilityFilenameStream << "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/outputdata/likelihoodjointpdfs/"
               << std::fixed << std::setprecision(2) << logTe << "_" << logL << ".txt";
    //saveProbabilitiesToText(averageLikelihood, probabilityFilenameStream.str());
    saveSmoothedProbabilitiesToText(smoothedCppArray, probabilityFilenameStream.str());

    return 0;
}
#endif // TEST_GAIALIKELIHOOD_MAIN -DTEST_GAIALIKELIHOOD_MAIN to compile with main function
