#ifndef ELIASHBERGCHAIN_PARAMETERS_H
#define ELIASHBERGCHAIN_PARAMETERS_H

#include <iostream>

#include <cstdlib>
#include <complex>

#include "globals.h"
#include <cassert>

class parameters {
    bool phonMediated;
    ulong NFreq;
    const ulong dimension;
    const ulong NK;
    const ulong Nx;
    const ulong Ny;
    const ulong NBands;
    double Beta;
    const double hubbardU;
    const std::complex<double> alphaseed;
    const ulong maxIterHF;
    std::string bandstructureFile;
    std::string lmcFile;
    const std::string photDispFile;
    const std::string couplingFile;
    double intrabandCoupling;
    double interbandCoupling;

    //const ulong NFREQBOSONIC = 2ul * NFreq - 1;

public:
    parameters(const bool phonMediatedIni,
               const ulong NFreqIni,
               const ulong dimensionIni,
               const ulong NKIni,
               const ulong NxIni,
               const ulong NyIni,
               const ulong NBandsIni,
               const double BetaIni,
               const double hubbardUIni,
               const std::complex<double> alphaseedIni,
               const ulong maxIterHFIni,
               const std::string bandstructureFileIni,
               const std::string lmcFileIni,
               const std::string photDispFileIni,
               const std::string couplingFileIni)
            : phonMediated(phonMediatedIni),
              NFreq(NFreqIni),
              dimension(dimensionIni),
              NK(NKIni),
              Nx(NxIni),
              Ny(NyIni),
              NBands(NBandsIni),
              Beta(BetaIni),
              hubbardU(hubbardUIni),
              alphaseed(alphaseedIni),
              maxIterHF(maxIterHFIni),
              bandstructureFile(bandstructureFileIni),
              lmcFile(lmcFileIni),
              photDispFile(photDispFileIni),
              couplingFile(couplingFileIni) {
        intrabandCoupling = 0.;
        interbandCoupling = 1.;
    }

    parameters(const bool phonMediatedIni,
               const ulong NFreqIni,
               const ulong dimensionIni,
               const ulong NKIni,
               const ulong NxIni,
               const ulong NyIni,
               const ulong NBandsIni,
               const double BetaIni,
               const double hubbardUIni,
               const std::complex<double> alphaseedIni,
               const ulong maxIterHFIni,
               const std::string bandstructureFileIni,
               const std::string lmcFileIni,
               const std::string photDispFileIni,
               const std::string couplingFileIni,
               const bool verbouse)
            : phonMediated(phonMediatedIni),
              NFreq(NFreqIni),
              dimension(dimensionIni),
              NK(NKIni),
              Nx(NxIni),
              Ny(NyIni),
              NBands(NBandsIni),
              Beta(BetaIni),
              hubbardU(hubbardUIni),
              alphaseed(alphaseedIni),
              maxIterHF(maxIterHFIni),
              bandstructureFile(bandstructureFileIni),
              lmcFile(lmcFileIni),
              photDispFile(photDispFileIni),
              couplingFile(couplingFileIni) {
        intrabandCoupling = 0.;
        interbandCoupling = 1.;
        if (verbouse) {
            std::cout << "Initialized parameters" << '\n';
            std::cout << "phonMediated = " << phonMediated << '\n';
            std::cout << "NFreq = " << NFreq << '\n';
            std::cout << "dimension = " << dimension << '\n';
            std::cout << "NK = " << NK << '\n';
            std::cout << "Nx = " << Nx << '\n';
            std::cout << "Ny = " << Ny << '\n';
            std::cout << "NBands = " << NBands << '\n';
            std::cout << "Beta = " << Beta << '\n';
            std::cout << "HubbardU = " << hubbardU << '\n';
            std::cout << "alphaseed = " << alphaseed << '\n';
            std::cout << "maxIterHF = " << maxIterHF << '\n';
            std::cout << "bandstructureFile = " << bandstructureFile << '\n';
            std::cout << "lmcFile = " << lmcFile << '\n';
            std::cout << "photDispFile = " << photDispFile << '\n';
            std::cout << "couplingFile = " << couplingFile << '\n';
        }

    }

    inline bool isPhonMediated() const { return phonMediated; }

    inline ulong getNFreq() const { return NFreq; }

    inline ulong getDimension() const { return dimension; }

    inline ulong getNK() const { return NK; }

    inline ulong getNx() const { return Nx; }

    inline ulong getNy() const { return Ny; }

    inline ulong getNBands() const { return NBands; }

    inline double getBeta() const { return Beta; }

    inline double getHubbardU() const { return hubbardU; }

    inline double getIntrabandCoupling() const { return intrabandCoupling; }

    inline double getInterbandCoupling() const { return interbandCoupling; }

    inline std::complex<double> getAlphaSeed() const { return alphaseed; }

    inline ulong getMaxIterHF() const { return maxIterHF; }

    inline ulong getNFreqBosonic() const { return 2ul * NFreq - 1ul; }

    inline std::string getBandstructureFile() const { return bandstructureFile; }

    inline std::string getLMCFile() const { return lmcFile; }

    inline std::string getPhotDispFile() const { return photDispFile; }

    inline std::string getCouplingFile() const { return couplingFile; }

    inline ulong getGFSize() const { return NK * NFreq * NBands * NBands * 4ul; }

    inline unsigned long greenArgumentToIndex(const unsigned long kInd,
                                              const unsigned long wInd,
                                              const unsigned long n1,
                                              const unsigned long i1,
                                              const unsigned long n2,
                                              const unsigned long i2) const {
        return kInd * 4ul * NBands * NBands * NFreq + wInd * 4ul * NBands * NBands + n1 * 4ul * NBands +
               i1 * 2ul * NBands + n2 * 2ul + i2;
    }

    inline unsigned long greenArgumentToIndex(const unsigned long kYInd,
                                              const unsigned long kXInd,
                                              const unsigned long wInd,
                                              const unsigned long n1,
                                              const unsigned long i1,
                                              const unsigned long n2,
                                              const unsigned long i2) const {
        return kYInd * 4ul * NBands * NBands * NFreq * Nx + kXInd * 4ul * NBands * NBands * NFreq +
               wInd * 4ul * NBands * NBands + n1 * 4ul * NBands +
               i1 * 2ul * NBands + n2 * 2ul + i2;
    }

    inline ulong greenArgumentToIndexPadded(const unsigned long kInd,
                                            const unsigned long wInd,
                                            const unsigned long n1,
                                            const unsigned long i1,
                                            const unsigned long n2,
                                            const unsigned long i2) const {
        return kInd * 4ul * NBands * NBands * 2ul * NFreq + wInd * 4ul * NBands * NBands + n1 * 4ul * NBands +
               i1 * 2ul * NBands + n2 * 2ul + i2;
    }

    inline ulong greenArgumentToIndexPadded(const unsigned long kYInd,
                                            const unsigned long kXInd,
                                            const unsigned long wInd,
                                            const unsigned long n1,
                                            const unsigned long i1,
                                            const unsigned long n2,
                                            const unsigned long i2) const {
        return kYInd * 4ul * NBands * NBands * 2ul * NFreq * Nx + kXInd * 4ul * NBands * NBands * 2ul * NFreq +
               wInd * 4ul * NBands * NBands +
               n1 * 4ul * NBands + i1 * 2ul * NBands + n2 * 2ul + i2;
    }


    inline unsigned long lmcArgumentToIndex(const ulong kInd,
                                            const ulong n1,
                                            const ulong n2,
                                            const ulong polarization) const {
        return kInd * NBands * NBands * dimension + n1 * NBands * dimension + n2 * dimension + polarization;
    }

    inline unsigned long lmcArgumentToIndex(const ulong kYInd,
                                            const ulong kXInd,
                                            const ulong n1,
                                            const ulong n2,
                                            const ulong polarization) const {
        return kYInd * 2ul * Nx * NBands * NBands * dimension + kXInd * NBands * NBands * dimension +
               n1 * NBands * dimension + n2 * dimension + polarization;
    }

    inline double wIndToW(const ulong wInd) const {
        return double(2. * double(wInd) - double(NFreq) + 1.) * PI / Beta;
    }

    inline double wIndToWExtended(const unsigned long wInd, const ulong extent) const {
        return double(2. * double(wInd) - (double(NFreq) * double(extent)) + 1) * PI / Beta;
    }

    inline double omegaIndToW(const unsigned long omegaInd) const {
        return 2. * (double(omegaInd) - double(NFreq) + 1.) * PI / Beta;
    }

    inline double omegaIndToWSmall(const ulong omegaInd) const {
        return (2. * double(omegaInd) - double(NFreq)) * PI / Beta;
    }

    inline ulong w1Mw2Small(const ulong w1, const ulong w2) const {
        assert(w1 + NFreq / 2ul >= w2);
        assert(w1 + NFreq / 2ul < NFreq + w2);
        return NFreq / 2ul + w1 - w2;
    }


    inline ulong minusKInd(const ulong kInd) const {
        return (NK - kInd) % NK;
    }

    inline ulong minusWInd(const ulong wInd) const {
        return NFreq - wInd - 1;
    }

    inline ulong minusWIndExtended(const ulong wInd, const ulong extent) const {
        return (NFreq * extent) - wInd - 1;
    }

    void setBeta(const double betaSet) {
        std::cout << "Setting Beta to new value: " << betaSet << '\n';
        Beta = betaSet;
    }

    void setNFreq(const ulong nFreqSet) {
        std::cout << "Setting NFreq to new value: " << nFreqSet << '\n';
        NFreq = nFreqSet;
    }

    void setIntrabandCoupling(const double intrabandCouplingSet) {
        std::cout << "Setting intrabandCoupling to new value: " << intrabandCouplingSet << '\n';
        intrabandCoupling = intrabandCouplingSet;
    }

    void setInterbandCoupling(const double interbandCouplingSet) {
        std::cout << "Setting interbandCoupling to new value: " << interbandCouplingSet << '\n';
        interbandCoupling = interbandCouplingSet;
    }

    void setBandstructureFile(const std::string bandstructureFileSet) {
        bandstructureFile = bandstructureFileSet;
    }

    void setLMCFile(const std::string lmcFileSet) {
        lmcFile = lmcFileSet;
    }

};


#endif //ELIASHBERGCHAIN_PARAMETERS_H
