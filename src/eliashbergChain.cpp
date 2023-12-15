#include <vector>
#include <complex>
#include <fileHandling.h>
#include "iostream"

#include "localOneGap.h"

#include "parameters.h"

int main() {
  
  ///////////////////// Computation of estimated critical temperature for graphene on STO ///////////////////////
  parameters prms = readInParameters("parameters/graphenePhonLoc.txt");
  const double d = 5.;//distance from interface in nm
  //The fraction of the entire BZ that is considered in the simulation
  const double fractionFromD = 0.5 * 1. / 2. * d / 0.142 * 4. * PI / 3.;
  //Strength of effective light-matter coupling excluding the division by sqrt(d)
  const double gTimesSqrtD = 0.0002144001430638806;
  const double couplingG = 1. / fractionFromD * gTimesSqrtD / std::sqrt(d * 1e-9);
  
  //Intra-band coupling is neglected
  prms.setIntrabandCoupling(0.);
  //Inter-band coupling is set to the previously calculated value
  prms.setInterbandCoupling(couplingG);
  std::cout << "g = " << couplingG << std::endl;
  //Frequency of surface mode of STO
  const double Om = 0.016125725310656633;
  //const double Om = 0.001125725310656633;
  //chemical potential in eV
  const double mu = -0.015;
  
  double betaC = betaCOneGapLocalCoupling(Om, mu, prms);
  std::cout << "mu = " << mu << std::endl;
  std::cout << "betaC = " << betaC << std::endl;
  
  ///////////// Computation for Sawtooth Chain ////////////////////
  //parameters prms = readInParameters("parameters/sawtoothPhonLoc.txt");
  //const ulong nOm (400ul);
  //const double mu = -1.5857864376269046;
  //std::vector<double> omArr (nOm, 0.);
  //for(ulong omInd = 0ul; omInd < omArr.size(); ++omInd){
  //  omArr[omInd] = 2.5 + 0.7 * double(omInd) / double(nOm);
  //}

  //tcAsOfOmOneGapLocal(mu, omArr, prms);
  //const ulong nNPH(200ul);
  //const double om = 2.97;
  //std::vector<double> nPhArr (nNPH, 0.);
  //for(ulong nPhInd = 0ul; nPhInd < nPhArr.size(); ++nPhInd){
  //  nPhArr[nPhInd] = 0. + 5. * double(nPhInd) / double(nPhArr.size());
  //}
  //tcAsOfNPHOneGapLocal(mu, om, nPhArr, prms);
  
  return 0;
}