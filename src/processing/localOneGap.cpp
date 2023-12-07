#include <vector>
#include <complex>

#include "parameters.h"
//#include "initialize.h"
//#include "dysonEq.h"
#include "fileHandling.h"
#include <math.h>
//#include "particleNumber.h"
#include <cmath>
#include "localOneGap.h"

//double drivenPhonons;

double gapOneDeltaLocalBig(const std::vector<std::complex<double>> &bands,
                           const double Om,
                           const double mu,
                           const double betaC,
                           const parameters &prms) {
  
  double w = PI / betaC;
  std::vector<std::complex<double>> sumKRes(prms.getNK(), std::complex<double>(0., 0.));
  
  double g4(prms.getInterbandCoupling() * prms.getInterbandCoupling() * prms.getInterbandCoupling() *
            prms.getInterbandCoupling());
  //double prefac(g4 * Om * Om / prms.getNK() / prms.getNK());
  double prefac(g4 * Om * Om / prms.getNK());
  
  const std::complex<double> del(0., 1e-8);
  
  //add attractive Hubbard U
  double UAttractive = 0.0;

#pragma omp parallel for
  for (long kInd = 0ul; kInd < long(prms.getNK()); ++kInd) {
    double E1 = abs(bands[kInd * prms.getNBands() + 0ul] - mu);
    //for (long kDInd = 0ul; kDInd < long(prms.getNK()); ++kDInd) {
    
    //double E2 = abs(bands[kDInd * prms.getNBands() + 1ul] - mu);
    double E2 = abs(bands[0ul * prms.getNBands() + 1ul] - mu);
    
    std::complex<double> nFE1Terms = ((Om * Om - (E1 - II * w) * (E1 - II * w) + 2 * (E1 * E1 - Om * Om - w * w) * nF(E1, betaC)) *
                                      (-((E1 + E2 - Om) * (E1 - E2 + Om) * (E2 + Om)) - 2 * E2 * (E1 * E1 - E2 * E2 + Om * Om) * nB(Om, betaC) + 2 * Om * (E1 * E1 + E2 * E2 - Om * Om) * nFE2(E2, betaC))) /
                                     ((E1 - E2 - Om) * (E1 + E2 - Om + del) * Om * (E1 - E2 + Om + del) * (E1 + E2 + Om) * (E1 * E1 * E1 * E1 + 2 * E1 * E1 * (-Om * Om + w * w) + (Om * Om + w * w) * (Om * Om + w * w)));
    
    std::complex<double> nFE2MOmTerms = (E1 * (nB(Om, betaC) + nFE2(E2, betaC)) * ((E2 + II * w) * (E2 - 2 * Om + II * w) + 2 * (-E2 * E2 + 2 * E2 * Om + w * w) * nF(-E2 + Om, betaC))) /
                                        ((E1 + E2 - Om + del) * Om * (E1 - E2 + Om + del) * (E2 * E2 + w * w) * ((E2 - 2 * Om) * (E2 - 2 * Om) + w * w));
    
    
    std::complex<double> nFE2POmTerms = (E1 * (1 + nB(Om, betaC) - nFE2(E2, betaC)) * (-(E2 * (E2 + 2 * Om)) + 2. * II * (E2 + Om) * w + w * w + 2 * (E2 * E2 + 2 * E2 * Om - w * w) * nF(E2 + Om, betaC))) /
                                        ((E1 - E2 - Om) * Om * (E1 + E2 + Om) * (E2 * E2 + w * w) * ((E2 + 2 * Om) * (E2 + 2 * Om) + w * w));
    
    
    std::complex<double> nB2OmTerms = (E1 * (4 * E2 * ((E1 - Om) * (E1 + Om) * (E2 * E2 * E2 * E2 - 6 * E2 * E2 * Om * Om + 8 * Om * Om * Om * Om) +
                                                       (E2 * E2 * E2 * E2 - 12 * E2 * E2 * Om * Om + 6 * Om * Om * Om * Om + 2 * E1 * E1 * (E2 * E2 + Om * Om)) * w * w +
                                                       (E1 * E1 + 2 * E2 * E2 - 3 * Om * Om) * w * w * w * w + w * w * w * w * w * w) * nB(Om, betaC) * nB(Om, betaC) +
                                             (E1 * E1 - (Om + II * w) * (Om + II * w)) * (E2 * E2 - (2 * Om + II * w) * (2 * Om + II * w)) *
                                             ((E2 + Om) * (E2 - II * w) * (E2 - 2 * Om + II * w) + 2 * Om * (E2 * E2 - 2. * II * Om * w - w * w) * nFE2(E2, betaC)) +
                                             2 * nB(Om, betaC) * (E2 * (E2 - 2 * Om) * (E1 - Om) * (E1 + Om) * (E2 + 2 * Om) * (2 * E2 * E2 - E2 * Om - 4 * Om * Om) -
                                                                  2. * II * E2 * Om * (E1 * E1 * E2 * E2 + E2 * E2 * E2 * E2 - 7 * E2 * E2 * Om * Om + 8 * Om * Om * Om * Om) * w +
                                                                  (2 * E2 * E2 * E2 * E2 * E2 - E2 * E2 * E2 * E2 * Om - 24 * E2 * E2 * E2 * Om * Om + 16 * E2 * E2 * Om * Om * Om + 12 * E2 * Om * Om * Om * Om -
                                                                   20 * Om * Om * Om * Om * Om +
                                                                   4 * E1 * E1 * (E2 * E2 * E2 + E2 * Om * Om + Om * Om * Om)) * w * w - 2. * II * E2 * Om * (E1 * E1 + 3 * E2 * E2 + Om * Om) * w * w * w +
                                                                  (4 * E2 * E2 * E2 - 6 * E2 * Om * Om - Om * Om * Om + E1 * E1 * (2 * E2 + Om)) * w * w * w * w - 4. * II * E2 * Om * w * w * w * w * w +
                                                                  (2 * E2 + Om) * w * w * w * w * w * w +
                                                                  2 * Om * (E2 * E2 * (E1 - Om) * (E1 + Om) * (E2 * E2 - 4 * Om * Om) + (E2 * E2 * E2 * E2 - 4 * (E1 * E1 + 4 * E2 * E2) * Om * Om + 20 * Om * Om * Om * Om) * w * w +
                                                                            (-E1 + Om) * (E1 + Om) * w * w * w * w - w * w * w * w * w * w) * nFE2(E2, betaC)))) /
                                      (Om * Om * (E2 * E2 + w * w) * (E1 * E1 * E1 * E1 + 2 * E1 * E1 * (-Om * Om + w * w) + (Om * Om + w * w) * (Om * Om + w * w)) *
                                       (E2 * E2 * E2 * E2 + 2 * E2 * E2 * (-4 * Om * Om + w * w) + (4 * Om * Om + w * w) * (4 * Om * Om + w * w)));
    
    sumKRes[kInd] += prefac / (E1 * E2) * (nFE1Terms + nFE2MOmTerms + nFE2POmTerms + nB2OmTerms);
    //}
    sumKRes[kInd] += UAttractive / double(prms.getNK()) * std::tanh(E1 * betaC / 2.) / (2. * E1);
  }
  
  std::complex<double> sum(0., 0.);
  for (ulong kInd = 0ul; kInd < prms.getNK(); ++kInd) {
    sum += sumKRes[kInd];
  }
  //std::cout << sum.imag() << std::endl;
  return sum.real();
}

double gapOneDeltaLocalBigCoupling(
    const std::vector<std::complex<double>> &bands,
    const std::vector<std::complex<double>> &coupling,
    const double Om,
    const double mu,
    const double betaC,
    const parameters &prms) {
  
  double w = PI / betaC;
  std::vector<std::complex<double>> sumKRes(prms.getNK(), std::complex<double>(0., 0.));
  
  double g4(prms.getInterbandCoupling() * prms.getInterbandCoupling() * prms.getInterbandCoupling() *
            prms.getInterbandCoupling());
  double prefac(g4 * Om * Om / prms.getNK() / prms.getNK());
  //double prefac(g4 * Om * Om / prms.getNK());
  
  const std::complex<double> del(0., 1e-8);
  
  //add attractive Hubbard U -- negative for repulsive
  const double d = 5.;//nm
  const double fractionFromD = 0.5 * 1. / 2. * d / 0.142 * 4. * PI / 3.;
  //double UAttractive = -9.3 * 0.2 / fractionFromD;
  double UAttractive = -1.1 / fractionFromD;

#pragma omp parallel for
  for (long kxInd = 0ul; kxInd < long(prms.getNx()); ++kxInd) {
    for (long kyInd = 0ul; kyInd < long(prms.getNy()); ++kyInd) {
      ulong kInd = kxInd * prms.getNx() + kyInd;
      double E1 = abs(bands[kInd * prms.getNBands() + 0ul] - mu);
      for (long kDxInd = 0ul; kDxInd < long(prms.getNx()); ++kDxInd) {
        for (long kDyInd = 0ul; kDyInd < long(prms.getNy()); ++kDyInd) {
          ulong kDInd = kDxInd * prms.getNx() + kDyInd;
          double E2 = abs(bands[kDInd * prms.getNBands() + 1ul] - mu);
          
          
          std::complex<double> nFE1Terms = ((Om * Om - (E1 - II * w) * (E1 - II * w) + 2 * (E1 * E1 - Om * Om - w * w) * nF(E1, betaC)) *
                                            (-((E1 + E2 - Om) * (E1 - E2 + Om) * (E2 + Om)) - 2 * E2 * (E1 * E1 - E2 * E2 + Om * Om) * nB(Om, betaC) + 2 * Om * (E1 * E1 + E2 * E2 - Om * Om) * nFE2(E2, betaC))) /
                                           ((E1 - E2 - Om) * (E1 + E2 - Om + del) * Om * (E1 - E2 + Om + del) * (E1 + E2 + Om) * (E1 * E1 * E1 * E1 + 2 * E1 * E1 * (-Om * Om + w * w) + (Om * Om + w * w) * (Om * Om + w * w)));
          
          std::complex<double> nFE2MOmTerms = (E1 * (nB(Om, betaC) + nFE2(E2, betaC)) * ((E2 + II * w) * (E2 - 2 * Om + II * w) + 2 * (-E2 * E2 + 2 * E2 * Om + w * w) * nF(-E2 + Om, betaC))) /
                                              ((E1 + E2 - Om + del) * Om * (E1 - E2 + Om + del) * (E2 * E2 + w * w) * ((E2 - 2 * Om) * (E2 - 2 * Om) + w * w));
          
          
          std::complex<double> nFE2POmTerms = (E1 * (1 + nB(Om, betaC) - nFE2(E2, betaC)) * (-(E2 * (E2 + 2 * Om)) + 2. * II * (E2 + Om) * w + w * w + 2 * (E2 * E2 + 2 * E2 * Om - w * w) * nF(E2 + Om, betaC))) /
                                              ((E1 - E2 - Om) * Om * (E1 + E2 + Om) * (E2 * E2 + w * w) * ((E2 + 2 * Om) * (E2 + 2 * Om) + w * w));
          
          
          std::complex<double> nB2OmTerms = (E1 * (4 * E2 * ((E1 - Om) * (E1 + Om) * (E2 * E2 * E2 * E2 - 6 * E2 * E2 * Om * Om + 8 * Om * Om * Om * Om) +
                                                             (E2 * E2 * E2 * E2 - 12 * E2 * E2 * Om * Om + 6 * Om * Om * Om * Om + 2 * E1 * E1 * (E2 * E2 + Om * Om)) * w * w +
                                                             (E1 * E1 + 2 * E2 * E2 - 3 * Om * Om) * w * w * w * w + w * w * w * w * w * w) * nB(Om, betaC) * nB(Om, betaC) +
                                                   (E1 * E1 - (Om + II * w) * (Om + II * w)) * (E2 * E2 - (2 * Om + II * w) * (2 * Om + II * w)) *
                                                   ((E2 + Om) * (E2 - II * w) * (E2 - 2 * Om + II * w) + 2 * Om * (E2 * E2 - 2. * II * Om * w - w * w) * nFE2(E2, betaC)) +
                                                   2 * nB(Om, betaC) * (E2 * (E2 - 2 * Om) * (E1 - Om) * (E1 + Om) * (E2 + 2 * Om) * (2 * E2 * E2 - E2 * Om - 4 * Om * Om) -
                                                                        2. * II * E2 * Om * (E1 * E1 * E2 * E2 + E2 * E2 * E2 * E2 - 7 * E2 * E2 * Om * Om + 8 * Om * Om * Om * Om) * w +
                                                                        (2 * E2 * E2 * E2 * E2 * E2 - E2 * E2 * E2 * E2 * Om - 24 * E2 * E2 * E2 * Om * Om + 16 * E2 * E2 * Om * Om * Om + 12 * E2 * Om * Om * Om * Om -
                                                                         20 * Om * Om * Om * Om * Om +
                                                                         4 * E1 * E1 * (E2 * E2 * E2 + E2 * Om * Om + Om * Om * Om)) * w * w - 2. * II * E2 * Om * (E1 * E1 + 3 * E2 * E2 + Om * Om) * w * w * w +
                                                                        (4 * E2 * E2 * E2 - 6 * E2 * Om * Om - Om * Om * Om + E1 * E1 * (2 * E2 + Om)) * w * w * w * w - 4. * II * E2 * Om * w * w * w * w * w +
                                                                        (2 * E2 + Om) * w * w * w * w * w * w +
                                                                        2 * Om * (E2 * E2 * (E1 - Om) * (E1 + Om) * (E2 * E2 - 4 * Om * Om) + (E2 * E2 * E2 * E2 - 4 * (E1 * E1 + 4 * E2 * E2) * Om * Om + 20 * Om * Om * Om * Om) * w * w +
                                                                                  (-E1 + Om) * (E1 + Om) * w * w * w * w - w * w * w * w * w * w) * nFE2(E2, betaC)))) /
                                            (Om * Om * (E2 * E2 + w * w) * (E1 * E1 * E1 * E1 + 2 * E1 * E1 * (-Om * Om + w * w) + (Om * Om + w * w) * (Om * Om + w * w)) *
                                             (E2 * E2 * E2 * E2 + 2 * E2 * E2 * (-4 * Om * Om + w * w) + (4 * Om * Om + w * w) * (4 * Om * Om + w * w)));
          
          
          const double b1x = 0.;
          const double b1y = 4. * PI / 3.;
          const double b2x = 4. * PI / (3. * std::sqrt(3.)) * 3. / 2.;
          const double b2y = 4. * PI / (3. * std::sqrt(3.)) * (-std::sqrt(3.) / 2.);
          
          const double k1x = 0.;
          const double k1y = b1y / prms.getNx() / fractionFromD * kxInd + b1y / 3. - b1y / fractionFromD / 2.;
          const double k2x = b2x / prms.getNx() / fractionFromD * kyInd + 2. * b2x / 3. - b2x / fractionFromD / 2.;
          const double k2y = b2y / prms.getNx() / fractionFromD * kyInd + 2. * b2y / 3. - b2y / fractionFromD / 2.;
          
          const double Kx = b1x / 3. + 2. * b2x / 3.;
          const double Ky = b1y / 3. + 2. * b2y / 3.;
          //KPoint = b1 / 3. + 2. * b2 / 3.
          
          const double dist = (k1x + k2x - Kx) * (k1x + k2x - Kx) + (k1y + k2y - Ky) * (k1y + k2y - Ky);
          
          const double couplingCalcK = std::sqrt(dist) < 0.142 / d ? 1. : 0.;
          
          const double kD1x = 0.;
          const double kD1y = b1y / prms.getNx() / fractionFromD * kDxInd + b1y / 3. - b1y / fractionFromD / 2.;
          const double kD2x = b2x / prms.getNx() / fractionFromD * kDyInd + 2. * b2x / 3. - b2x / fractionFromD / 2.;
          const double kD2y = b2y / prms.getNx() / fractionFromD * kDyInd + 2. * b2y / 3. - b2y / fractionFromD / 2.;
          
          const double distQ = (k1x + k2x - kD1x - kD2x) * (k1x + k2x - kD1x - kD2x) + (k1y + k2y - kD1y - kD2y) * (k1y + k2y - kD1y - kD2y);
          
          const double couplingCalcQ = std::sqrt(distQ) < 0.142 / d ? 1. : 0.;
          
          //std::complex<double> couplingK = coupling[kInd];
          std::complex<double> couplingKD = coupling[kDInd];
          
          sumKRes[kInd] += couplingCalcQ * couplingKD * couplingCalcK * prefac / (E1 * E2) * (nFE1Terms + nFE2MOmTerms + nFE2POmTerms + nB2OmTerms);
        }
      }
      std::complex<double> couplingK = coupling[kInd];
      sumKRes[kInd] += couplingK * UAttractive / double(prms.getNK()) * std::tanh(E1 * betaC / 2.) / (2. * E1);
    }
  }
  
  std::complex<double> sum(0., 0.);
  for (ulong kInd = 0ul; kInd < prms.getNK(); ++kInd) {
    sum += sumKRes[kInd];
  }
  //std::cout << sum.imag() << std::endl;
  return sum.real();
}

double gapOneDeltaLocalShort(const std::vector<std::complex<double>> &bands,
                             const double Om,
                             const double mu,
                             const double betaC,
                             const parameters &prms) {
  
  std::vector<std::complex<double>> sumKRes(prms.getNK(), std::complex<double>(0., 0.));
  
  double g4(prms.getInterbandCoupling() * prms.getInterbandCoupling() * prms.getInterbandCoupling() *
            prms.getInterbandCoupling());
  //double prefac(g4 / prms.getNK() / prms.getNK());
  double prefac(g4 / prms.getNK());
  
  const std::complex<double> del(0., 1e-6);

#pragma omp parallel for
  for (long kInd = 0ul; kInd < long(prms.getNK()); ++kInd) {
    double E1 = abs(bands[kInd * prms.getNBands() + 0ul] - mu);
    //for (long kDInd = 0ul; kDInd < long(prms.getNK()); ++kDInd) {
    //double E2 = abs(bands[kDInd * prms.getNBands() + 1ul] - mu);
    double E2 = abs(bands[0ul * prms.getNBands() + 1ul] - mu);
    
    std::complex<double> bcsTerms =
        (1. + nB(Om, betaC) - nFE2(E2, betaC)) * (1. - 2. * nF(E1, betaC))
        / (2. * E1 * Om * Om * Om);
    
    std::complex<double> resonantTerms =
        -(nB(Om, betaC) + nFE2(E2, betaC)) * (E1 * (1. - 2. * nF(E2 - Om, betaC)) - (E2 - Om) * (1. - 2. * nF(E1, betaC)))
        / (E1 * ((E2 - Om) * (E2 - Om) - E1 * E1 + del) * Om * Om);
    
    sumKRes[kInd] += prefac * (bcsTerms + resonantTerms);
    //}
  }
  
  std::complex<double> sum(0., 0.);
  for (ulong kInd = 0ul; kInd < prms.getNK(); ++kInd) {
    sum += sumKRes[kInd];
  }
  //std::cout << sum.imag() << std::endl;
  return sum.real();
}

double betaCOneGapLocal(const double Om, const double mu, const parameters &prms) {
  std::vector<std::complex<double>> bandStructure(prms.getNK() * prms.getNBands(), std::complex<double>(0., 0.));
  readInComplex2DArray(bandStructure, prms.getBandstructureFile());
  assert(bandStructure.size() == prms.getNK() * prms.getNBands());
  
  double betaCMax(1e12);
  double betaCMin(0.);
  ulong maxIter(100ul);
  ulong iter(0ul);
  
  double betaC(1.);
  double resultGapEquation(0.);
  
  while (abs(resultGapEquation - 1.) > 1e-7 && iter < maxIter) {
    
    resultGapEquation = gapOneDeltaLocalBig(bandStructure, Om, mu, betaC, prms);
    //resultGapEquation = gapOneDeltaLocalShort(bandStructure, Om, mu, betaC, prms);
    
    if (resultGapEquation < 1) {
      if (iter < 10) {
        betaC += 20.;
        ++iter;
        continue;
      }
      betaCMin = betaC;
      betaC = (betaC + betaCMax) / 2.;
    } else {
      betaCMax = betaC;
      betaC = (betaC + betaCMin) / 2.;
    }
    if (abs(betaC - betaCMax) < 1e-6)
      break;
    ++iter;
  }
  std::cout << "result Gap Eq. = " << resultGapEquation << '\n';
  return betaC;
}

double betaCOneGapLocalApprox(const double Om, const double mu, const parameters &prms) {
  std::vector<std::complex<double>> bandStructure(prms.getNK() * prms.getNBands(), std::complex<double>(0., 0.));
  readInComplex2DArray(bandStructure, prms.getBandstructureFile());
  assert(bandStructure.size() == prms.getNK() * prms.getNBands());
  
  double betaCMax(1e12);
  double betaCMin(0.);
  ulong maxIter(100ul);
  ulong iter(0ul);
  
  double betaC(1.);
  double resultGapEquation(0.);
  
  while (abs(resultGapEquation - 1.) > 1e-7 && iter < maxIter) {
    
    resultGapEquation = gapOneDeltaLocalShort(bandStructure, Om, mu, betaC, prms);
    
    if (resultGapEquation < 1) {
      if (iter < 10) {
        betaC += 20.;
        ++iter;
        continue;
      }
      betaCMin = betaC;
      betaC = (betaC + betaCMax) / 2.;
    } else {
      betaCMax = betaC;
      betaC = (betaC + betaCMin) / 2.;
    }
    if (abs(betaC - betaCMax) < 1e-6)
      break;
    ++iter;
  }
  std::cout << "result Gap Eq. = " << resultGapEquation << '\n';
  return betaC;
}

double betaCOneGapLocalCoupling(const double Om, const double mu, const parameters &prms) {
  std::vector<std::complex<double>> bandStructure(prms.getNK() * prms.getNBands(), std::complex<double>(0., 0.));
  std::vector<std::complex<double>> coupling(prms.getNK(), std::complex<double>(0., 0.));
  readInComplex2DArray(bandStructure, prms.getBandstructureFile());
  readInComplex1DArray(coupling, prms.getLMCFile());
  assert(bandStructure.size() == prms.getNK() * prms.getNBands());
  
  double betaCMax(1e6);
  double betaCMin(0.);
  ulong maxIter(100ul);
  ulong iter(0ul);
  
  double betaC(1.);
  double resultGapEquation(0.);
  
  while (abs(resultGapEquation - 1.) > 1e-7 && iter < maxIter) {
    
    //std::cout << "res Gap = " << resultGapEquation << '\n';
    //std::cout << "betaC = " << betaC << '\n';
    
    resultGapEquation = gapOneDeltaLocalBigCoupling(bandStructure, coupling, Om, mu, betaC, prms);
    
    if (resultGapEquation < 1) {
      if (iter < 10) {
        betaC += 100.;
        ++iter;
        continue;
      }
      betaCMin = betaC;
      betaC = (betaC + betaCMax) / 2.;
    } else {
      betaCMax = betaC;
      betaC = (betaC + betaCMin) / 2.;
    }
    if (abs(betaC - betaCMax) < 1e-6)
      break;
    ++iter;
  }
  std::cout << "result Gap Eq. = " << resultGapEquation << '\n';
  return betaC;
}

double betaCOneGapApproxLocal(const double Om, const double mu, const parameters &prms) {
  std::vector<std::complex<double>> bandStructure(prms.getNK() * prms.getNBands(), std::complex<double>(0., 0.));
  readInComplex2DArray(bandStructure, prms.getBandstructureFile());
  assert(bandStructure.size() == prms.getNK() * prms.getNBands());
  
  double betaCMax(1e12);
  double betaCMin(0.);
  ulong maxIter(100ul);
  ulong iter(0ul);
  
  double betaC(1.);
  double resultGapEquation(0.);
  
  while (abs(resultGapEquation - 1.) > 1e-7 && iter < maxIter) {
    
    //resultGapEquation = gapOneDeltaLocalBig(bandStructure, Om, mu, betaC, prms);
    resultGapEquation = gapOneDeltaLocalShort(bandStructure, Om, mu, betaC, prms);
    
    if (resultGapEquation < 1) {
      if (iter < 10) {
        betaC += 20.;
        ++iter;
        continue;
      }
      betaCMin = betaC;
      betaC = (betaC + betaCMax) / 2.;
    } else {
      betaCMax = betaC;
      betaC = (betaC + betaCMin) / 2.;
    }
    if (abs(betaC - betaCMax) < 1e-6)
      break;
    ++iter;
  }
  std::cout << "result Gap Eq. = " << resultGapEquation << '\n';
  return betaC;
}


void tcAsOfMuOneGapLocal(const std::vector<double> &muArr, const double Om, parameters &prms) {
  
  prms.setIntrabandCoupling(0.);
  prms.setInterbandCoupling(1.);
  std::vector<double> betaCArr(muArr.size(), 0.);
  
  //double coupling = std::sqrt(1. / 137. * 3 * 1e8 / (1e12 * 1e-9 * 24));
  //std::cout << "coupling = " << coupling << '\n';
  
  for (ulong muInd = 0ul; muInd < muArr.size(); ++muInd) {
    std::cout << "Starting search at mu = " << muArr[muInd] << '\n';
    betaCArr[muInd] = betaCOneGapLocal(Om, muArr[muInd], prms);
    std::cout << "betaC = " << betaCArr[muInd] << '\n';
  }
  
  std::string BosStr("Bos" + std::to_string(int(100 * drivenPhonons)));
  std::string ElStr("El" + std::to_string(int(1000 * drivenElectrons)));
  
  std::string filename("data/results/tcAsOfMuOneGapHubU");
  filename = filename + ElStr;
  filename = filename + BosStr + ".hdf5";
  std::cout << "filename = " << filename << std::endl;
  writeCandBToHdf5(muArr, betaCArr, filename, prms);
  
}


void tcAsOfOmOneGapLocal(const double mu, const std::vector<double> &omArr, parameters &prms) {
  
  prms.setIntrabandCoupling(0.);
  prms.setInterbandCoupling(0.6);
  std::vector<double> betaCArr(omArr.size(), 0.);
  std::vector<double> betaCApproxArr(omArr.size(), 0.);
  
  for (ulong omInd = 0ul; omInd < omArr.size(); ++omInd) {
    std::cout << "Starting search at Om = " << omArr[omInd] << '\n';
    betaCArr[omInd] = betaCOneGapLocal(omArr[omInd], mu, prms);
    betaCApproxArr[omInd] = betaCOneGapLocalApprox(omArr[omInd], mu, prms);
    std::cout << "betaC = " << betaCArr[omInd] << '\n';
  }
  
  std::string BosStr("Bos" + std::to_string(int(100 * drivenPhonons)));
  
  std::string filename("data/results/tcAsOfOmOneGap");
  filename = filename + BosStr + ".hdf5";
  std::cout << "filename = " << filename << std::endl;
  writeCandBToHdf5(omArr, betaCArr, filename, prms);
  
  std::string filenameA("data/results/tcAsOfOmOneGapShort");
  filenameA = filenameA + BosStr + ".hdf5";
  std::cout << "filename Approx = " << filenameA << std::endl;
  writeCandBToHdf5(omArr, betaCApproxArr, filenameA, prms);
  
}

void tcAsOfNPHOneGapLocal(const double mu, const double om, const std::vector<double> &NphArr, parameters &prms) {
  
  prms.setIntrabandCoupling(0.);
  prms.setInterbandCoupling(0.6);
  std::vector<double> betaCArr(NphArr.size(), 0.);
  
  for (ulong nPhInd = 0ul; nPhInd < NphArr.size(); ++nPhInd) {
    ////// ---- set numbre of phonons ------ ////////
    //drivenPhonons = NphArr[nPhInd];
    std::cout << "Starting search at nPh = " << NphArr[nPhInd] << '\n';
    //betaCArr[nPhInd] = betaCOneGapLocalCoupling(om, mu, prms);
    betaCArr[nPhInd] = betaCOneGapLocal(om, mu, prms);
    std::cout << "betaC = " << betaCArr[nPhInd] << '\n';
  }
  
  std::string OmStr("Bos" + std::to_string(int(100 * om)));
  
  std::string filename("data/results/tcAsOfNphOneGap");
  filename = filename + OmStr + ".hdf5";
  std::cout << "filename = " << filename << std::endl;
  writeCandBToHdf5(NphArr, betaCArr, filename, prms);
  
}