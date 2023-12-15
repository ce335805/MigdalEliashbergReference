#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <string>
//#include <H5public.h>
#include <algorithm>
#include <cassert>
#include <complex>
#include <fstream>

#include "fileHandling.h"
#include "globals.h"
#include "parameters.h"


#include "H5Cpp.h"


/**
 *
 * @param data - 1D data array to be written to output
 * @param filename - data will be given a default name
 */
void writeReal1DArrayToHdf5(const std::vector<double> data, const std::string filename) {

  H5::H5File file(filename, H5F_ACC_TRUNC);

  const hsize_t dataSize = data.size();
  H5::DataSpace dataSpace(1, &dataSize);
  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);

  H5::DataSet dataset = file.createDataSet("Data", datatype, dataSpace);
  dataset.write(&data[0], datatype);
}

void writeTwoReal1DArrayToHdf5(const std::vector<double> dat1,
                               const std::vector<double> dat2,
                               const std::string filename,
                               const std::string datName1,
                               const std::string datName2) {

  H5::H5File file(filename, H5F_ACC_TRUNC);

  const hsize_t data1Size = dat1.size();
  const hsize_t data2Size = dat2.size();
  H5::DataSpace dataSpace1(1, &data1Size);
  H5::DataSpace dataSpace2(1, &data2Size);
  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);

  H5::DataSet dataset1 = file.createDataSet(datName1, datatype, dataSpace1);
  H5::DataSet dataset2 = file.createDataSet(datName2, datatype, dataSpace2);
  dataset1.write(&dat1[0], datatype);
  dataset2.write(&dat2[0], datatype);
}

void writeBetasAndDsToHdf5(const std::vector<double> betas,
                               const std::vector<int> ds,
                               const std::vector<double> gaps,
                               const std::string filename) {

  std::vector<double> dsDouble (ds.size(), 0.);
  for(ulong dInd = 0ul; dInd < dsDouble.size(); ++dInd){
    dsDouble[dInd] = double(ds[dInd]);
  }

  H5::H5File file(filename, H5F_ACC_TRUNC);

  const hsize_t betasSize = betas.size();
  const hsize_t dsSize = ds.size();
  H5::DataSpace dataSpace1(1, &betasSize);
  H5::DataSpace dataSpace2(1, &dsSize);
  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);

  H5::DataSet datasetBetas = file.createDataSet("Betas", datatype, dataSpace1);
  H5::DataSet datasetDs = file.createDataSet("Ds", datatype, dataSpace2);
  datasetBetas.write(&betas[0], datatype);
  datasetDs.write(&dsDouble[0], datatype);

  const hsize_t dataShapeGaps[3] = {dsSize, betasSize, 2ul};
  H5::DataSpace dataSpaceGaps(3, dataShapeGaps);
  datatype.setOrder(H5T_ORDER_LE);

  H5::DataSet datasetGaps = file.createDataSet("Gaps", datatype, dataSpaceGaps);
  datasetGaps.write(&gaps[0], datatype);

}

void writeBetasToHdf5(const std::vector<double> betas,
                           const std::vector<double> gaps,
                           const std::string filename,
                           const parameters &prms) {

  H5::H5File file(filename, H5F_ACC_TRUNC);

  const hsize_t betasSize = betas.size();
  H5::DataSpace dataSpace1(1, &betasSize);
  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);

  H5::DataSet datasetBetas = file.createDataSet("Betas", datatype, dataSpace1);
  datasetBetas.write(&betas[0], datatype);

  const hsize_t dataShapeGaps[2] = {betasSize, prms.getNBands()};
  H5::DataSpace dataSpaceGaps(2, dataShapeGaps);
  datatype.setOrder(H5T_ORDER_LE);

  H5::DataSet datasetGaps = file.createDataSet("Gaps", datatype, dataSpaceGaps);
  datasetGaps.write(&gaps[0], datatype);

}

void writeCandBToHdf5(const std::vector<double> couplings,
                      const std::vector<double> betaCs,
                      const std::string filename,
                      const parameters &prms) {

  H5::H5File file(filename, H5F_ACC_TRUNC);

  const hsize_t couplingsSize = couplings.size();
  H5::DataSpace dataSpace1(1, &couplingsSize);
  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);

  H5::DataSet datasetCouplings = file.createDataSet("Couplings", datatype, dataSpace1);
  datasetCouplings.write(&couplings[0], datatype);
  
  const hsize_t betaCsSize = betaCs.size();
  const hsize_t dataShapeBetaCs[1] = {betaCsSize};
  H5::DataSpace dataSpaceBetaCs(1, dataShapeBetaCs);
  datatype.setOrder(H5T_ORDER_LE);

  H5::DataSet datasetBetaCs = file.createDataSet("BetaCs", datatype, dataSpaceBetaCs);
  datasetBetaCs.write(&betaCs[0], datatype);

}


void readInComplex1DArray(std::vector<std::complex<double>> &readInArray, const std::string &fileName) {
  H5::H5File file(fileName, H5F_ACC_RDONLY);

  H5::DataSet realDataset = file.openDataSet("Real");
  H5::DataSet imagDataset = file.openDataSet("Imag");

  H5T_class_t typeClass = realDataset.getTypeClass();
  assert(typeClass == H5T_FLOAT);
  typeClass = imagDataset.getTypeClass();
  assert(typeClass == H5T_FLOAT);

  H5::DataSpace realDataSpace = realDataset.getSpace();
  H5::DataSpace imagDataSpace = imagDataset.getSpace();

  int rank = realDataSpace.getSimpleExtentNdims();
  assert(rank == 1);
  rank = imagDataSpace.getSimpleExtentNdims();
  assert(rank == 1);

  hsize_t dimsOut[1];
  realDataSpace.getSimpleExtentDims(dimsOut, nullptr);
  assert(dimsOut[0] == readInArray.size());
  imagDataSpace.getSimpleExtentDims(dimsOut, nullptr);
  assert(dimsOut[0] == readInArray.size());

  std::vector<double> realInput(readInArray.size(), 0.0);
  std::vector<double> imagInput(readInArray.size(), 0.0);

  const hsize_t inDimension[1] = {readInArray.size()};
  H5::DataSpace memspace(1, inDimension);

  realDataset.read(&realInput[0], H5::PredType::NATIVE_DOUBLE, memspace, realDataSpace);
  imagDataset.read(&imagInput[0], H5::PredType::NATIVE_DOUBLE, memspace, imagDataSpace);

  for (auto ind = 0ul; ind < readInArray.size(); ++ind) {
    readInArray[ind] = std::complex<double>(realInput[ind], imagInput[ind]);
  }
}

double readInDouble(const std::string &fileName, const std::string &dataName) {
  H5::H5File file(fileName, H5F_ACC_RDONLY);

  std::vector<double> readInArray(1ul, 0.);

  H5::DataSet realDataset = file.openDataSet(dataName);

  H5T_class_t typeClass = realDataset.getTypeClass();
  assert(typeClass == H5T_FLOAT);
  assert(typeClass == H5T_FLOAT);

  H5::DataSpace realDataSpace = realDataset.getSpace();

  int rank = realDataSpace.getSimpleExtentNdims();
  assert(rank == 1);

  hsize_t dimsOut[1];
  realDataSpace.getSimpleExtentDims(dimsOut, nullptr);
  assert(dimsOut[0] == readInArray.size());

  const hsize_t inDimension[1] = {readInArray.size()};
  H5::DataSpace memspace(1, inDimension);

  realDataset.read(&readInArray[0], H5::PredType::NATIVE_DOUBLE, memspace, realDataSpace);

  return readInArray[0];
}



/**
 *
 * @param readInArray array which is filled with to be read data
 * @param fileName of to be read data
 */
void readInComplex2DArray(std::vector<std::complex<double>> &readInArray, const std::string &fileName) {
  std::cout << "Filename='" << fileName << std::endl;
  H5::H5File file(fileName, H5F_ACC_RDONLY);

  H5::DataSet realDataset = file.openDataSet("Real");
  H5::DataSet imagDataset = file.openDataSet("Imag");

  H5T_class_t typeClass = realDataset.getTypeClass();
  assert(typeClass == H5T_FLOAT);
  typeClass = imagDataset.getTypeClass();
  assert(typeClass == H5T_FLOAT);

  H5::DataSpace realDataSpace = realDataset.getSpace();
  H5::DataSpace imagDataSpace = imagDataset.getSpace();

  int rank = realDataSpace.getSimpleExtentNdims();
  assert(rank == 2);
  rank = imagDataSpace.getSimpleExtentNdims();
  assert(rank == 2);

  hsize_t dimsOut[2];
  realDataSpace.getSimpleExtentDims(dimsOut, nullptr);
  assert(dimsOut[0] * dimsOut[1] == readInArray.size());
  imagDataSpace.getSimpleExtentDims(dimsOut, nullptr);
  assert(dimsOut[0] * dimsOut[1] == readInArray.size());

  std::vector<double> realInput(readInArray.size(), 0.0);
  std::vector<double> imagInput(readInArray.size(), 0.0);

  const hsize_t inDimension[1] = {readInArray.size()};
  H5::DataSpace memspace(1, inDimension);

  realDataset.read(&realInput[0], H5::PredType::NATIVE_DOUBLE, memspace, realDataSpace);
  imagDataset.read(&imagInput[0], H5::PredType::NATIVE_DOUBLE, memspace, imagDataSpace);

  for (auto ind = 0ul; ind < readInArray.size(); ++ind) {
    readInArray[ind] = std::complex<double>(realInput[ind], imagInput[ind]);
  }
}

/**
 *
 * @param readInArray array which is filled with to be read data
 * @param fileName of to be read data
 */
void readInComplex3DArray(std::vector<std::complex<double>> &readInArray, const std::string &fileName) {
  H5::H5File file(fileName, H5F_ACC_RDONLY);

  H5::DataSet realDataset = file.openDataSet("Real");
  H5::DataSet imagDataset = file.openDataSet("Imag");

  H5T_class_t typeClass = realDataset.getTypeClass();
  assert(typeClass == H5T_FLOAT);
  typeClass = imagDataset.getTypeClass();
  assert(typeClass == H5T_FLOAT);

  H5::DataSpace realDataSpace = realDataset.getSpace();
  H5::DataSpace imagDataSpace = imagDataset.getSpace();

  int rank = realDataSpace.getSimpleExtentNdims();
  assert(rank == 3);
  rank = imagDataSpace.getSimpleExtentNdims();
  assert(rank == 3);

  hsize_t dimsOut[3];
  realDataSpace.getSimpleExtentDims(dimsOut, nullptr);
  assert(dimsOut[0] * dimsOut[1] * dimsOut[2] == readInArray.size());
  imagDataSpace.getSimpleExtentDims(dimsOut, nullptr);
  assert(dimsOut[0] * dimsOut[1] * dimsOut[2] == readInArray.size());

  std::vector<double> realInput(readInArray.size(), 0.0);
  std::vector<double> imagInput(readInArray.size(), 0.0);

  const hsize_t inDimension[1] = {readInArray.size()};
  H5::DataSpace memspace(1, inDimension);

  realDataset.read(&realInput[0], H5::PredType::NATIVE_DOUBLE, memspace, realDataSpace);
  imagDataset.read(&imagInput[0], H5::PredType::NATIVE_DOUBLE, memspace, imagDataSpace);

  for (auto ind = 0ul; ind < readInArray.size(); ++ind) {
    readInArray[ind] = std::complex<double>(realInput[ind], imagInput[ind]);
  }
}

void readInComplex4DArray(std::vector<std::complex<double>> &readInArray, const std::string &fileName) {
    H5::H5File file(fileName, H5F_ACC_RDONLY);

    H5::DataSet realDataset = file.openDataSet("Real");
    H5::DataSet imagDataset = file.openDataSet("Imag");

    H5T_class_t typeClass = realDataset.getTypeClass();
    assert(typeClass == H5T_FLOAT);
    typeClass = imagDataset.getTypeClass();
    assert(typeClass == H5T_FLOAT);

    H5::DataSpace realDataSpace = realDataset.getSpace();
    H5::DataSpace imagDataSpace = imagDataset.getSpace();

    int rank = realDataSpace.getSimpleExtentNdims();
    assert(rank == 4);
    rank = imagDataSpace.getSimpleExtentNdims();
    assert(rank == 4);

    hsize_t dimsOut[4];
    realDataSpace.getSimpleExtentDims(dimsOut, nullptr);
    assert(dimsOut[0] * dimsOut[1] * dimsOut[2] * dimsOut[3] == readInArray.size());
    imagDataSpace.getSimpleExtentDims(dimsOut, nullptr);
    assert(dimsOut[0] * dimsOut[1] * dimsOut[2] * dimsOut[3] == readInArray.size());

    std::vector<double> realInput(readInArray.size(), 0.0);
    std::vector<double> imagInput(readInArray.size(), 0.0);

    const hsize_t inDimension[1] = {readInArray.size()};
    H5::DataSpace memspace(1, inDimension);

    realDataset.read(&realInput[0], H5::PredType::NATIVE_DOUBLE, memspace, realDataSpace);
    imagDataset.read(&imagInput[0], H5::PredType::NATIVE_DOUBLE, memspace, imagDataSpace);

    for (auto ind = 0ul; ind < readInArray.size(); ++ind) {
        readInArray[ind] = std::complex<double>(realInput[ind], imagInput[ind]);
    }
}

void readInGFObject(std::vector<std::complex<double>> &readInArray, const std::string &fileName) {

  H5::H5File file(fileName, H5F_ACC_RDONLY);

  H5::DataSet realDataset = file.openDataSet("Real");
  H5::DataSet imagDataset = file.openDataSet("Imag");

  H5T_class_t typeClass = realDataset.getTypeClass();
  assert(typeClass == H5T_FLOAT);
  typeClass = imagDataset.getTypeClass();
  assert(typeClass == H5T_FLOAT);

  H5::DataSpace realDataSpace = realDataset.getSpace();
  H5::DataSpace imagDataSpace = imagDataset.getSpace();

  int rank = realDataSpace.getSimpleExtentNdims();
  assert(rank == 6);
  rank = imagDataSpace.getSimpleExtentNdims();
  assert(rank == 6);

  hsize_t dimsOut[6];
  realDataSpace.getSimpleExtentDims(dimsOut, nullptr);
  assert(dimsOut[0] * dimsOut[1] * dimsOut[2] * dimsOut[3] * dimsOut[4] * dimsOut[5] == readInArray.size());
  imagDataSpace.getSimpleExtentDims(dimsOut, nullptr);
  assert(dimsOut[0] * dimsOut[1] * dimsOut[2] * dimsOut[3] * dimsOut[4] * dimsOut[5] == readInArray.size());

  std::vector<double> realInput(readInArray.size(), 0.0);
  std::vector<double> imagInput(readInArray.size(), 0.0);

  const hsize_t inDimension[1] = {readInArray.size()};
  H5::DataSpace memspace(1, inDimension);

  realDataset.read(&realInput[0], H5::PredType::NATIVE_DOUBLE, memspace, realDataSpace);
  imagDataset.read(&imagInput[0], H5::PredType::NATIVE_DOUBLE, memspace, imagDataSpace);

  for (auto ind = 0ul; ind < readInArray.size(); ++ind) {
    readInArray[ind] = std::complex<double>(realInput[ind], imagInput[ind]);
  }
}

void writeGFObjectToHdf5(const std::vector<std::complex<double>> &data,
                         const std::string &filename,
                         const parameters &prms) {

  assert(data.size() == prms.getGFSize());

  std::vector<double> realData(prms.getGFSize(), 0.0);
  std::vector<double> imagData(prms.getGFSize(), 0.0);

  std::transform(data.begin(),
                 data.end(),
                 realData.begin(),
                 [](const std::complex<double> entry) -> double {
                     return entry.real();
                 });
  std::transform(data.begin(),
                 data.end(),
                 imagData.begin(),
                 [](const std::complex<double> entry) -> double {
                     return entry.imag();
                 });

  H5::H5File file(filename, H5F_ACC_TRUNC);

  const hsize_t dataShape[6] = {prms.getNK(), prms.getNFreq(), prms.getNBands(), 2, prms.getNBands(), 2};
  H5::DataSpace dataSpace(6, dataShape);
  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);

  H5::DataSet datasetReal = file.createDataSet("Real", datatype, dataSpace);
  H5::DataSet datasetImag = file.createDataSet("Imag", datatype, dataSpace);
  datasetReal.write(&realData[0], datatype);
  datasetImag.write(&imagData[0], datatype);

}

parameters readInParameters(const std::string &fileName) {

  std::string line;

  bool phonMediatedRead(false);
  ulong NFreqRead(0ul);
  ulong NKRead(0ul);
  ulong NxRead(0ul);
  ulong NyRead(0ul);
  ulong NBandsRead(0ul);
  double BetaRead(0.);
  double hubbardURead(0.);
  std::complex<double> alphaSeedRead(0., 0.);
  double speedOfLightRead(0.);
  ulong maxIterHFRead(0ul);
  std::string bandstructureFileRead;
  std::string lmcFileRead;
  std::string photDispFileRead;
  std::string couplingFileRead;

  std::ifstream parameterFile(fileName);
  if (parameterFile.is_open()) {
    getline(parameterFile, line);
    phonMediatedRead = bool(std::stoi(line));

    getline(parameterFile, line);
    NFreqRead = ulong(std::stoi(line));

    getline(parameterFile, line);
    NKRead = ulong(std::stoi(line));

    getline(parameterFile, line);
    NxRead = ulong(std::stoi(line));

    getline(parameterFile, line);
    NyRead = ulong(std::stoi(line));

    getline(parameterFile, line);
    NBandsRead = ulong(std::stoi(line));

    getline(parameterFile, line);
    BetaRead = std::stod(line);

    getline(parameterFile, line);
    hubbardURead = std::stod(line);

    getline(parameterFile, line);
    alphaSeedRead = std::complex<double>(std::stod(line), 0.);

    getline(parameterFile, line);
    maxIterHFRead = ulong(std::stoi(line));

    getline(parameterFile, line);
    bandstructureFileRead = line;

    getline(parameterFile, line);
    lmcFileRead = line;

    getline(parameterFile, line);
    photDispFileRead = line;

    getline(parameterFile, line);
    couplingFileRead = line;

    parameterFile.close();
  } else std::cout << "Unable to open parameter file";

  ulong dimensionRead(0ul);
  if(NxRead == 1 || NyRead == 1){
      dimensionRead = 1ul;
  } else {
      dimensionRead = 2ul;
  }
  parameters readInParameters(phonMediatedRead, NFreqRead, dimensionRead, NKRead, NxRead, NyRead, NBandsRead, BetaRead, hubbardURead, alphaSeedRead, maxIterHFRead,
                              bandstructureFileRead, lmcFileRead, photDispFileRead, couplingFileRead, true);

  return readInParameters;

}


