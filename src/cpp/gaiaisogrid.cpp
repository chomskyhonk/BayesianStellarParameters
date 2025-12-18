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
bool parseIsochroneFilenameIso(const std::string& filename, double& age, double& Z) {
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
std::vector<Isochrone> readIsochroneFilesIso(const std::string& folderPath) {
    std::vector<Isochrone> isochrones; // Stores all isochrones
    std::vector<IsochroneFileMetadata> fileMetadataList; // Temporary container for sorting

    //step 1:iterate through files and extract metadata
    for (const auto& entry : fs::directory_iterator(folderPath)) {
        std::cout << "Found file: " << entry.path() << std::endl;
        if (entry.path().extension() == ".isc_gaia-dr3") { // Check for the correct file extension
            double age, Z;
            std::string filename = entry.path().filename().string();
            std::cout << "Checking filename: " << filename << std::endl;

            if (parseIsochroneFilenameIso(filename, age, Z)) {
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

void printIsochroneData(const std::vector<Isochrone>& isochrones) {
    for (const auto& isochrone : isochrones) {
        std::cout << "Age: " << isochrone.age << ", [M/H]: " << isochrone.MH << ", Z: " << isochrone.Z << std::endl;
        size_t numValuesToPrint = std::min(isochrone.initialMass.size(), static_cast<size_t>(10));
        for (size_t i = 0; i < numValuesToPrint; ++i) { //just prints the first ten values for ^^ example
            std::cout << "Initial Mass: " << isochrone.initialMass[i]
                      << ", Mass: " << isochrone.mass[i]
                      << ", logLLo: " << isochrone.logLLo[i]
                      << ", logTe: " << isochrone.logTe[i]
                      << ", G: " << isochrone.G[i]
                      << ", Gbp: " << isochrone.Gbp[i]
                      << ", Grp: " << isochrone.Grp[i] << std::endl;
        }
    }
}

void calculateWeights(std::vector<Isochrone>& isochrones, std::vector<std::vector<double>>& weights, const std::string& filename, ParameterSet paramSet) { //do not need const because we are modifying the isochrones vector
    size_t numIsochrones = isochrones.size();
    if (numIsochrones == 0) return;//failsafe 

    size_t numIsochronePoints = isochrones[0].initialMass.size();

    std::ofstream outFile(filename);
    outFile << "Age, MH, Z, initialMass, Mass, ";
    if (paramSet == logTe_logLLo) {
        outFile << "logLLo, logTe, ";
    } else {
        outFile << "G, Gbp, ";
    }
    outFile << "isoweight\n";

//i indexeds the isochrones, j indexes the data points on the isochrones
    for (size_t i = 0; i < numIsochrones; ++i) { //iterating over all isochrones
        for (size_t j = 0; j < numIsochronePoints; ++j) { //iterating over all data points on current isochrone
            double ageValue = isochrones[i].age; //stores age value for current isochrone
            double MHvalue = isochrones[i].MH; //stores value of MH for the current isochrone
            double param1Value, param2Value;
            if (paramSet == logTe_logLLo) {
                param1Value = isochrones[i].logTe[j];
                param2Value = isochrones[i].logLLo[j];
            } else {
                param1Value = isochrones[i].Gbp[j];
                param2Value = isochrones[i].G[j];
            }
            double massBefore = (j > 0) ? isochrones[i].initialMass[j - 1] : isochrones[i].initialMass[j]; //handles edge case for the first point
            double massAfter = (j < numIsochronePoints - 1) ? isochrones[i].initialMass[j + 1] : isochrones[i].initialMass[j]; //handles edge case for the last point
            double deltaM = 0.5 * std::abs(massAfter - massBefore); //calculates deltaM for the current data point

            double minAgeDistance = std::numeric_limits<double>::max(); //initialises minDistance to the maximum value of a double
            double minMHDistance = std::numeric_limits<double>::max(); //initialises minDistance to the maximum value of a double
            double nearestAge = 0.0; //initialises nearestAge to 0.0
            double nearestMH = 0.0; //initialises nearestMH to 0.0
            
            for (size_t l = 0; l < numIsochrones; ++l) { //iterating over all isochrones again to compare paramspace distances
                if (l != i && isochrones[l].MH != MHvalue && isochrones[l].age != ageValue) { //checks if the current isochrone is not the same as the ith isochrone and has a different metallicity/age
                //this still compares isochrones with potentially equal metallicities
                //but its fine because the weight values are 0 and wont contribute to the sum
                    double ageDistance;
                    if (paramSet == logTe_logLLo) {
                        ageDistance = std::sqrt(
                            std::pow(isochrones[l].logTe[j] - param1Value, 2) + //add sigma values for logTe and logLLo
                            std::pow(isochrones[l].logLLo[j] - param2Value, 2)
                        );
                    } else {
                        ageDistance = std::sqrt(
                            std::pow(isochrones[l].Gbp[j] - param1Value, 2) + //add sigma values for BV and Mv
                            std::pow(isochrones[l].G[j] - param2Value, 2)
                        );
                    }
                    //calculates the distance between the current data point and the data point at index l on the kth isochrone

                    if (ageDistance < minAgeDistance) { //checks if the age distance is less than the minimum age distance
                        minAgeDistance = ageDistance; //if it is, sets the minimum age distance to the age distance
                        nearestAge = isochrones[l].age; //sets the nearest age to the age of the kth isochrone
                        nearestMH = isochrones[l].MH;                          
                    }
                }
            }
                
            double deltaT = 0.5 * std::abs(ageValue - nearestAge);
            double deltaMH = 0.5 * std::abs(MHvalue - nearestMH);
            
            //update delta values to fit isochrone array structure
            isochrones[i].deltaT[j] = deltaT;
            isochrones[i].deltaMH[j] = deltaMH;
            isochrones[i].deltaM[j] = deltaM;
            
            //calculate the weight for the current data point
            weights[i][j] = deltaT * deltaMH * deltaM;
            
            //printing everything for debugging
            std::cout << "Isochrone " << i << ", Point " << j << ": deltaT = " << deltaT
                      << ", deltaMH = " << deltaMH << " (| " << MHvalue << " - " << nearestMH << " |)"
                      << ", deltaM = " << deltaM
                      << ", isoweight = " << weights[i][j] << std::endl;
            //write to file
            outFile << ageValue << ", " << MHvalue << ", " << isochrones[i].Z << ", "
                    << isochrones[i].initialMass[j] << ", " << isochrones[i].mass[j] << ", ";
            if (paramSet == logTe_logLLo) {
                outFile << isochrones[i].logLLo[j] << ", " << isochrones[i].logTe[j] << ", ";
            } else {
                outFile << isochrones[i].G[j] << ", " << isochrones[i].Gbp[j] << ", ";
            }
            outFile << weights[i][j] << "\n";
        }
    }
    outFile.close();
}

bool loadWeights(const std::string& filename, std::vector<Isochrone>& isochrones, std::vector<std::vector<double>>& weights, ParameterSet paramSet) {
    constexpr double epsilon = 1e-4;

    auto areEqual = [](double a, double b) {
        return std::abs(a - b) < epsilon;
    };

    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }

    std::string line;
    std::getline(inFile, line); // Skip header line

    size_t isoIndex = 0;
    size_t weightIndex = 0;

    while (std::getline(inFile, line)) {

        //std::cout << "Reading line: " << line << std::endl;  //debug output
        //line.erase(0, line.find_first_not_of(" \t\r\n")); // Trim left
        //line.erase(line.find_last_not_of(" \t\r\n") + 1); // Trim right
        std::istringstream iss(line);
        double age, MH, Z, initialMass, mass, param1, param2, weight;
        char comma;
        if (paramSet == logTe_logLLo) {
            if (!(iss >> age >> comma >> MH >> comma >> Z >> comma >> initialMass >> comma >> mass 
                >> comma >> param2 >> comma >> param1 >> comma >> weight)) {
                //std::cerr << "Error parsing line: " << line << std::endl;
                continue;
            }
        } else {
            if (!(iss >> age >> comma >> MH >> comma >> Z >> comma >> initialMass >> comma >> mass 
                >> comma >> param2 >> comma >> param1 >> comma >> weight)) {
                //std::cerr << "Error parsing line: " << line << std::endl;
                continue;
            }
        }

        // Iterate through isochrones in order
        while (isoIndex < isochrones.size() &&
               (!areEqual(isochrones[isoIndex].age, age) || 
                !areEqual(isochrones[isoIndex].MH, MH) || 
                !areEqual(isochrones[isoIndex].Z, Z))) {
            isoIndex++;
            weightIndex = 0; // Reset weight index for the new isochrone
        }

        if (isoIndex < isochrones.size()) {
            // Check for matching initialMass
            if (weightIndex < isochrones[isoIndex].initialMass.size() &&
                areEqual(isochrones[isoIndex].initialMass[weightIndex], initialMass)) {
                // Assign the weight
                weights[isoIndex][weightIndex] = weight;
                //std::cout << "Assigned weight=" << weight << " to weights[" << isoIndex << "][" << weightIndex << "]" << std::endl;
                weightIndex++;
            } else {
                std::cerr << "Warning: Data point not found in isochrone at age=" << age
                          << ", MH=" << MH << ", Z=" << Z
                          << ", initialMass=" << initialMass << std::endl;
            }
        } else {
            std::cerr << "Error: No matching isochrone found for age=" << age
                      << ", MH=" << MH << ", Z=" << Z << std::endl;
            break;
        }
    }

    inFile.close();
    return true;
}

double gaussian(double param1Residual, double param2Residual, double param1sigma, double param2sigma) {
    return (1.0 / (2 * M_PI * param1sigma * param2sigma)) * std::exp(-0.5 * std::pow(param1Residual / (param1sigma), 2) -0.5 * std::pow(param2Residual / (param2sigma), 2));
}

double gaussian2(double param1Residual, double param2Residual, double param3Residual, double param4Residual, double param1sigma, double param2sigma, double param3sigma, double param4sigma) {
    return (1.0 / (std::pow(2 * M_PI, 2) * param1sigma * param2sigma * param3sigma * param4sigma)) * std::exp(-0.5 * std::pow(param1Residual / (param1sigma), 2) -0.5 * std::pow(param2Residual / (param2sigma), 2) -0.5 * std::pow(param3Residual / (param3sigma), 2) -0.5 * std::pow(param4Residual / (param4sigma), 2));
}

int findBinIndex(double value, double minValue, double binWidth, int numBins) {
    int binIndex = static_cast<int>((value - minValue) / binWidth);
    if (binIndex < 0) binIndex = 0;
    if (binIndex >= numBins) binIndex = numBins - 1;
    return binIndex;
}//for a given value, this function calculates the bin index it belongs to
//value-minValue gives the distance from the minimum value, divide by binwidth 
//to determine which bin value falls under. converts to int to floor the value and get bin index

std::vector<std::vector<double>> binningisogrid(const std::vector<Isochrone>& isochrones, const std::vector<std::vector<double>>& weights, double param1, double param2, double param3, double param4, const std::string& filename, ParameterSet paramSet) {
    size_t numIsochrones = isochrones.size();
    size_t numIsochronePoints = isochrones[0].initialMass.size();
    const int numAgeBins = 50;
    const int numMHBins = 48;
    std::vector<std::vector<double>> jointPDF(numAgeBins, std::vector<double>(numMHBins, 0.0));//2D array to store joint PDF values, size of agebins*mhbins, initialises all values to 0.0
    std::vector<std::vector<bool>> binPopulated(numAgeBins, std::vector<bool>(numMHBins, false)); //2D array to track populated bins
    std::vector<std::set<double>> distinctAgesInBins(numAgeBins);
    std::vector<std::set<double>> distinctMHInBins(numMHBins);
    //bin ranges:
    double minAge = 0.1;
    double maxAge = 14.0;
    double minMH = -1.5;
    double maxMH = 0.42;
    //bin width:
    double ageBinWidth = (maxAge - minAge) / numAgeBins;
    double MHBinWidth = (maxMH - minMH) / numMHBins;

    std::random_device rd;//randomly selects points for debugging
    std::mt19937 gen(rd());

    //ensure each bin is populated
    for (double age = minAge; age <= maxAge; age += ageBinWidth) {
        for (double MH = minMH; MH <= maxMH;  MH += MHBinWidth) {
            int ageBinIndex = findBinIndex(age, minAge, ageBinWidth, numAgeBins);
            int MHBinIndex = findBinIndex(MH, minMH, MHBinWidth, numMHBins); //finds the bin index for the current MH value
            binPopulated[ageBinIndex][MHBinIndex] = true;
        }
    }//corresponding bin in the binpopulated 2d array is marked as true/false 
    //depending on whether it is populated or not

    //open the file before the loop - saving probabilities to txt file for debugging
    std::ofstream probFile(filename);
    probFile << "Probability\n";  //add a header before writing values
    if (!probFile) {
        //std::cerr << "Error opening probabilities.txt for writing" << std::endl;
    }
    
    //iterating over all isochrones to bin age/MH
    for (size_t i = 0; i < numIsochrones; ++i) {
        int ageBinIndex = findBinIndex(isochrones[i].age, minAge, ageBinWidth, numAgeBins);
        int MHBinIndex = findBinIndex(isochrones[i].MH, minMH, MHBinWidth, numMHBins);
        //std::cout << "Isochrone " << i << ": ageBinIndex = " << ageBinIndex << ", MHBinIndex = " << MHBinIndex << std::endl;
        
        //mark bin as populated
        binPopulated[ageBinIndex][MHBinIndex] = true;
        //track distinct age bins
        distinctAgesInBins[ageBinIndex].insert(isochrones[i].age);
        distinctMHInBins[MHBinIndex].insert(isochrones[i].MH);
        //generate a random selection of 5 points for debugging purposes
        std::vector<size_t> indices(numIsochronePoints);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), gen);
        indices.resize(std::min<size_t>(5, numIsochronePoints));

        
        for (size_t j = 0; j < numIsochronePoints; ++j) {
            double param1Value, param2Value, param3Value, param4Value;
            if (paramSet == logTe_logLLo) {
                param1Value = isochrones[i].logTe[j];
                param2Value = isochrones[i].logLLo[j];
            //} else {
                param3Value = isochrones[i].Gbp[j] - isochrones[i].Grp[j]; //colour = Gbp - Grp
                param4Value = isochrones[i].G[j];
            }

            double weight = weights[i][j];

            double param1Residual = param1 - param1Value;
            double param2Residual = param2 - param2Value;
            double param3Residual = param3 - param3Value;
            double param4Residual = param4 - param4Value;


            double param1Sigma = 0.05; //B-V / logTe / set ~0.10 and ~0.02-0.03 respectively
            double param2Sigma = 0.05; //Mv / logLLo / set  and ~0.03-0.06 respectively
            double param3Sigma = 0.1; //Gbp-Grp 0.1
            double param4Sigma = 0.2; //M_G 0.1
            //read in and add baseline uncertainty

            double gaussianValue = gaussian2(param1Residual, param2Residual, param3Residual, param4Residual, param1Sigma, param2Sigma, param3Sigma, param4Sigma);
            //double gaussianValue = gaussian(param1Residual, param2Residual, param1Sigma, param2Sigma);
            double probability = weight * gaussianValue;
            
            jointPDF[ageBinIndex][MHBinIndex] += probability; //represents a 2D histogram where each bin contains the cumulative probability of all isochrone points within its age-MH range

            //write probability to the file
            if (probFile) {
                probFile << probability << "\n";
            }

            //printing residuals for debugging
            if (std::find(indices.begin(), indices.end(), j) != indices.end()) {
                //std::cout << "Isochrone " << i << ", Point " << j << ": Parameter 1 Residual = " << param1Residual
                //          << ", Parameter 2 Resiudal = " << param2Residual
                //          << ", Parameter 3 Residual = " << param3Residual
                //          << ", Parameter 4 Residual = " << param4Residual
                //          << ", gaussianValue = " << gaussianValue
                //          << ", probability = " << probability << std::endl;
            }
        }
    }
    //close the file after the loop
    probFile.close();

    //normalise the joint PDF
    double totalSum = 0.0;
    for (int ageBinIndex = 0; ageBinIndex < numAgeBins; ++ageBinIndex) {
        for (int MHBinIndex = 0; MHBinIndex < numMHBins; ++MHBinIndex) {
            totalSum += jointPDF[ageBinIndex][MHBinIndex];
        } 
    }
    if (totalSum > 0) {
        for (int ageBinIndex = 0; ageBinIndex < numAgeBins; ++ageBinIndex) {
            for (int MHBinIndex = 0; MHBinIndex < numMHBins; ++MHBinIndex) {
                jointPDF[ageBinIndex][MHBinIndex] /= totalSum;
            }
        }
    }

    double normalizedSum = 0.0;
    //std::cout << "Normalized Joint PDF:" << std::endl;
    for (int ageBinIndex = 0; ageBinIndex < numAgeBins; ++ageBinIndex) {
        for (int MHBinIndex = 0; MHBinIndex < numMHBins; ++MHBinIndex) {
            //std::cout << "[Age Bin " << ageBinIndex << ", MH Bin " << MHBinIndex << ": " << jointPDF[ageBinIndex][MHBinIndex] << "] ";
            normalizedSum += jointPDF[ageBinIndex][MHBinIndex];
        }
        //std::cout << std::endl;
    }
    //std::cout << "Sum of all probabilities: " << normalizedSum << std::endl;


    //print the joint PDF for debugging
    //std::cout << "Joint PDF:" << std::endl;
    //for (int ageBinIndex = 0; ageBinIndex < numAgeBins; ++ageBinIndex) {
    //    for (int MHBinIndex = 0; MHBinIndex < numMHBins; ++MHBinIndex) {
    //        std::cout << "[Age Bin " << ageBinIndex << ", MH Bin " << MHBinIndex << ": " << jointPDF[ageBinIndex][MHBinIndex] << "]";
    //    }
    //    std::cout << std::endl;
    //}

    //print unpopulated bins
    //std::cout << "Unpopulated Bins:" << std::endl;
    for (int ageBinIndex = 0; ageBinIndex < numAgeBins; ++ageBinIndex) {
        for (int MHBinIndex = 0; MHBinIndex < numMHBins; ++MHBinIndex) {
            if (!binPopulated[ageBinIndex][MHBinIndex]) {
                //std::cout << "Unpopulated Bin - Age Bin " << ageBinIndex << ", MH Bin " << MHBinIndex << std::endl;
            }
        }
    }

   //std::cout << "Number of distinct isochrone metallicities contributing to each MH bin:" << std::endl;
    for (int MHBinIndex = 0; MHBinIndex < numMHBins; ++MHBinIndex) {
        //std::cout << "MH Bin " << MHBinIndex << ": " << distinctMHInBins[MHBinIndex].size() << " distinct metallicities" << std::endl;
    }
    //std::cout << "Number of distinct isochrone ages contributing to each age bin:" << std::endl;
    for (int ageBinIndex = 0; ageBinIndex < numAgeBins; ++ageBinIndex) {
        //std::cout << "Age Bin " << ageBinIndex << ": " << distinctAgesInBins[ageBinIndex].size() << " distinct ages" << std::endl;
    }
    
    return jointPDF;
}

void plotJointPDF(const std::vector<std::vector<double>>& jointPDF, double minAge, double maxAge, double minMH, double maxMH, double param1, double param2, double param3, double param4, ParameterSet paramSet, bool logscale = false) {
    py::scoped_interpreter guard{}; //start the interpreter and keep it alive
    
    const int numAgeBins = jointPDF.size();//no. of elements in outer vector in jointPDF - no. of age bins
    const int numMHBins = jointPDF[0].size();//no of elements in first inner vector of jointPDF - no. of MH bins
    //debugging this ^^
    //std::cout << "numAgeBins: " << numAgeBins << ", numMHBins: " << numMHBins << std::endl;
    //flatten the 2D jointPDF vector into a 1D vector of floats - required for numpy array
    std::vector<float> flattenedJointPDF; 
    for (const auto& row : jointPDF) {
        for (double value : row) {
            flattenedJointPDF.push_back(static_cast<float>(value));
        }
    }
    //debugging flattenedjointpdf - should equal numAgeBins*numMHBins
    //std::cout << "Flattened jointPDF size: " << flattenedJointPDF.size() << std::endl;
    //convert flattenedJointPDF to a numpy array
    py::array_t<float> numpyArray(flattenedJointPDF.size(), flattenedJointPDF.data());
    //reshape the numpy array to the correct dimensions
    numpyArray = numpyArray.attr("reshape")(numAgeBins, numMHBins);

    //import matplotlib and pyplot and colors
    py::module_ plt = py::module_::import("matplotlib.pyplot");
    py::module_ font_manager = py::module_::import("matplotlib.font_manager");
    py::module_ colors = py::module_::import("matplotlib.colors");

    double minNonzero = std::numeric_limits<double>::max();
    for (const auto& row : jointPDF) {
        for (double v : row) {
            if (v > 0.0 && v < minNonzero) minNonzero = v;
        }
    }
    

    //create the plot
    if (logscale) {
        plt.attr("imshow")(numpyArray, py::arg("extent") = py::make_tuple(minMH, maxMH, minAge, maxAge), py::arg("origin") = "lower", py::arg("aspect") = "auto", py::arg("norm") = colors.attr("LogNorm")(py::arg("vmin") = 1e-50));
    } else {
        plt.attr("imshow")(numpyArray, py::arg("extent") = py::make_tuple(minMH, maxMH, minAge, maxAge), py::arg("origin") = "lower", py::arg("aspect") = "auto");
    }
    plt.attr("colorbar")();
    plt.attr("ylabel")("Age (Gyr)");
    plt.attr("xlabel")("Metallicity [M/H]");
    plt.attr("title")("Joint PDF of Age and Metallicity");

    //determine parameter names based on enum parameterset
    std::string param1Name, param2Name, param3Name, param4Name;
    std::string folderName;
    if (paramSet == logTe_logLLo) {
        param1Name = "logTe";
        param2Name = "logLLo";
        folderName = "JointPDF_logTe_logLLo";
    //} else {
        param3Name = "G_bp - G_rp";
        param4Name = "G-mag";
        //folderName = "JointPDF_BV_Mv";
    }
    //legend for specified star
    std::ostringstream legendStream;
    legendStream << std::fixed << std::setprecision(2); 
    legendStream << "Specified Star:\n" << param1Name << " = " << param1 << "\n" << param2Name << " = " << param2 << "\n" << param3Name << " = " << param3 << "\n" << param4Name << " = " << param4;
    std::string legendText = legendStream.str();
    py::object font_properties = font_manager.attr("FontProperties")(py::arg("size") = 8); //set font size to 8
    //0.98, 0.98 for top right corner, 0.23, 0.98 for top left corner
    plt.attr("text")(0.29, 0.98, legendText, py::arg("transform") = plt.attr("gca")().attr("transAxes"), py::arg("fontsize") = 8, py::arg("verticalalignment") = "top", py::arg("horizontalalignment") = "right", py::arg("bbox") = py::dict(py::arg("facecolor") = "white", py::arg("alpha") = 0.75));
    std::string outputDir = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/plots/gaia/" + folderName;
    std::filesystem::create_directory(outputDir);
    std::ostringstream filenameStream;
    filenameStream << outputDir << "/jointpdf_" << param1Name << param1 << "_" << param2Name << param2 << "_" << param3Name << param3 << "_" << param4Name << param4 << (logscale ? "_logscale" : "_linearscale") << ".png";
    // Add logscale info to the plot title
    std::string scaleText = logscale ? " (log scale)" : " (linear scale)";
    plt.attr("title")(("Joint PDF of Age and Metallicity" + scaleText));
    std::string filename = filenameStream.str();
    plt.attr("savefig")(filename);
    plt.attr("show")();
}

void saveProbabilitiesToTextIso(const std::vector<std::vector<double>>& jointPDF, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    for (const auto& row : jointPDF) {
        for (const auto& value : row) {
            outFile << value << " ";
        }
        outFile << "\n";
    }
    outFile.close();
}

#ifdef TEST_GAIAISOGRID_MAIN
int main() {
    std::string folderPath = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/isochrones2";
    std::vector<Isochrone> isochrones = readIsochroneFilesIso(folderPath);//stores result of this function call as the variable isochrones
    size_t numIsochrones = isochrones.size(); //no. of isochrone objects in the isochrones vector
    std::cout << "Number of isochrones: " << numIsochrones << std::endl;
    size_t numDataPoints = isochrones[0].initialMass.size(); //no. of data points in the initialmass vector of the isochrone with index 0
    std::cout << "Number of data points in a single isochrone: " << numDataPoints << std::endl;
    //printIsochroneData(isochrones); //calls the isochrones variable as defined above

    std::vector<std::vector<double>> weights(numIsochrones, std::vector<double>(numDataPoints, 0.0));
    ParameterSet paramSet = logTe_logLLo;  //change to Mv_BV for plotting Mv vs B-V

    double specifiedlogTe = 3.82; //^^
    double specifiedlogLLo = 1.80; //not great data but still use from GAIA
    double specifiedG = -1.50; //G
    double specifiedGbpGrp = 0.60; //Gbp - Grp

    //create a stringstream to format logTe and logLLo or Mv and B-V. set precision to 3sf
    std::ostringstream param1Stream, param2Stream, param3Stream, param4Stream;
    if (paramSet == logTe_logLLo) {
        param1Stream << std::fixed << std::setprecision(2) << specifiedlogTe;
        param2Stream << std::fixed << std::setprecision(2) << specifiedlogLLo;
    //} else {
        param3Stream << std::fixed << std::setprecision(2) << specifiedGbpGrp;
        param4Stream << std::fixed << std::setprecision(2) << specifiedG;
    }

    std::string weightsFilename = folderPath + "/weights_" + (paramSet == logTe_logLLo ? "logTe_logLLo" : "Gbp_G") + ".csv";
    std::string probabilitiesFilename = folderPath + "/probabilities/probabilities_" + (paramSet == logTe_logLLo ? "logTe" : "G_bp") + param1Stream.str() + "_" + (paramSet == logTe_logLLo ? "logLLo" : "G") + param2Stream.str() + ".txt";
    if (!loadWeights(weightsFilename, isochrones, weights, paramSet)) {
        calculateWeights(isochrones, weights, weightsFilename, paramSet);
    }
    //calculateWeights(isochrones, weights, weightsFilename, paramSet);
     std::vector<std::vector<double>> jointPDF;
    if (paramSet == logTe_logLLo) {
        jointPDF = binningisogrid(isochrones, weights, specifiedlogTe, specifiedlogLLo, specifiedGbpGrp, specifiedG, probabilitiesFilename, logTe_logLLo);
        plotJointPDF(jointPDF, 0.0, 14.0, -1.5, 0.5, specifiedlogTe, specifiedlogLLo, specifiedGbpGrp, specifiedG, logTe_logLLo, true); //true for logscale
    //} else {
    //    jointPDF = binning(isochrones, weights, specifiedBV, specifiedMv, probabilitiesFilename, BV_Mv);
    //    plotJointPDF(jointPDF, 0.0, 14.0, -2.5, 0.5, specifiedBV, specifiedMv, BV_Mv, false);
    }
    saveProbabilitiesToTextIso(jointPDF, "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/outputdata/isogridjointpdfs/" + param1Stream.str() + "_" + param2Stream.str() + "_" + param3Stream.str() + "_" + param4Stream.str() + ".txt");
    return 0;
}
#endif // TEST_GAIAISOGRID_MAIN -DTEST_GAIAISOGRID_MAIN to compile with main function