#ifndef ELIASHBERGCHAIN_FILEHANDLING_H
#define ELIASHBERGCHAIN_FILEHANDLING_H

#include <vector>
#include <complex>
#include <string>
#include "parameters.h"

/**
 *
 * @param data - 1D data array to be written to output
 * @param filename - data will be given a default name
 */
void writeReal1DArrayToHdf5(const std::vector<double> data, const std::string filename);

void writeTwoReal1DArrayToHdf5(const std::vector<double> dat1,
                               const std::vector<double> dat2,
                               const std::string filename,
                               const std::string datName1,
                               const std::string datName2);

void writeBetasAndDsToHdf5(const std::vector<double> betas,
                           const std::vector<int> ds,
                           const std::vector<double> gaps,
                           const std::string filename);

void writeBetasToHdf5(const std::vector<double> betas,
                      const std::vector<double> gaps,
                      const std::string filename,
                      const parameters &prms);

void writeCandBToHdf5(const std::vector<double> couplings,
                      const std::vector<double> betaCs,
                      const std::string filename,
                      const parameters &prms);

double readInDouble(const std::string &fileName, const std::string &dataName);

/**
 *
 * @param readInArray array which is to be filled with data
 * @param fileName name of file to read
 */
void readInComplex1DArray(std::vector<std::complex<double>> &readInArray, const std::string &fileName);


/**
 *
 * @param readInArray array which is filled with to be read data
 * @param fileName of to be read data
 */
void readInComplex2DArray(std::vector<std::complex<double>> &readInArray, const std::string &fileName);


/**
 *
 * @param readInArray array which is filled with to be read data
 * @param fileName of to be read data
 */
void readInComplex3DArray(std::vector<std::complex<double>> &readInArray, const std::string &fileName);

void readInComplex4DArray(std::vector<std::complex<double>> &readInArray, const std::string &fileName);

void readInGFObject(std::vector<std::complex<double>> &readInArray, const std::string &fileName);

/**
 * Write object that has dimensions of a GF to a hdf5 file
 * @param data
 * @param filename
 */
void writeGFObjectToHdf5(const std::vector<std::complex<double>> &data,
                         const std::string &filename,
                         const parameters &prms);


/**
 * Read in parameters from a given file
 * @param fileName
 * @return the parameters in a specified class
 */
parameters readInParameters(const std::string& fileName);


#endif //ELIASHBERGCHAIN_FILEHANDLING_H
