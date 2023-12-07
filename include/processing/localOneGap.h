#ifndef ELIASHBERGCHAIN_LOCALONEGAP_H
#define ELIASHBERGCHAIN_LOCALONEGAP_H

#include "globals.h"


inline double nF(const double E, const double beta) { return 1. / (std::exp(beta * E) + 1); }

inline double nFE2(const double E, const double beta) { return 1. / (std::exp(beta * E) + 1) + drivenElectrons; }


inline double nB(const double E, const double beta) {
  return 1. / (std::exp(beta * E) - 1) + drivenPhonons;
}

double gapOneDeltaLocalBig(const std::vector<std::complex<double>> &bands,
                           const double Om,
                           const double mu,
                           const double betaC,
                           const parameters &prms);

double gapOneDeltaQ(const std::vector<std::complex<double>> &bands,
                    const double Om,
                    const double mu,
                    const double betaC,
                    const parameters &prms);

double betaCOneGapLocal(const double Om, const double mu, const parameters &prms);

double betaCOneGapApproxLocal(const double Om, const double mu, const parameters &prms);

double betaCOneGapLocalCoupling(const double Om, const double mu, const parameters &prms);

void tcAsOfMuOneGapLocal(const std::vector<double> &muArr, const double Om, parameters &prms);

void tcAsOfOmOneGapLocal(const double mu, const std::vector<double> &omArr, parameters &prms);

void tcAsOfNPHOneGapLocal(const double mu, const double om, const std::vector<double> &NphArr, parameters &prms);

#endif //ELIASHBERGCHAIN_LOCALONEGAP_H
