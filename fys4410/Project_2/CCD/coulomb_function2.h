#ifndef COULOMB_FUNCTIONS2_H
#define COULOMB_FUNCTIONS2_H

#include <armadillo>
#include "abstract_coulomb.h"

using namespace arma;

// struct channelBundle{
//   int chanNx;
//   int chanNy;
//   int chanNz;
//   int chanSz;
//   int chanTz;
// };

class coulomb_function2: public abstract_coulomb{
public:

  double hbarOmega;
  int shellMax;
  // int Nspstates;
  // int Nparticles;
  // int Nchannels;
  // int ** indexMap;
  // double * spEnergy;
  // channelBundle * chanValue;
  // channelBundle * chanModValue;

  coulomb_function2(double hbarOmega_In, int shellMax_In, int Nparticles_In);// : abstractSPbasis(densityIn, tzMaxIn, shellMaxIn, NparticlesIn) {}
  void generateIndexMap();
  void generateBasis();
  int checkSympqrs(int p, int q, int r, int s);
  int checkModSympqrs(int p, int q, int r, int s);
  int checkChanSym(int p, int q, int ichan);
  int checkChanSym2(int p, int q, int chanMl, int chanSz);
int checkChanModSym2(int p, int q, int chanMl, int chanSz);
  int checkChanModSym(int p, int q, int ichan);
  void setUpTwoStateChannels();
  void printBasis();
  void deallocate();

  double calc_TBME(int p, int q, int r, int s , mat mapping);
  int spinExchangeMtxEle(int i, int j, int k, int l);
  int kron_del(int i, int j);
  double specialGamma(double x);
  double bin_coef(int n, int k);
  double fac_over_fac(int a1, int a2);
  int factorial(int n);
};

#endif /* QDOTS2D_HPP */
