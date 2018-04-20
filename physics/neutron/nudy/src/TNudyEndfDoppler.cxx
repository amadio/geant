#include <math.h>
#include <iostream>
#include "Geant/TNudyEndfDoppler.h"
#include <iomanip>
#include <algorithm>

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfDoppler)
#endif

    TNudyEndfDoppler::TNudyEndfDoppler()
{
}

TNudyEndfDoppler::TNudyEndfDoppler(double isigDiff, double aw, double t1, double t2, std::vector<double> &x1,
                                   std::vector<double> &x2)
    : fF0K2P(0), fF1K2P(0), fF2K2P(0), fF3K2P(0), fF4K2P(0)
{
#define XNEPSM(S) fabs(fONE - fTW3 * sqrt(fONE / 3.0) * pow((fONE + S * (fONE + S)), fTHH) / (S * (fONE + S)))
#define FTAIL(X, fY)                                                                 \
  fOVSQPI *((1.0 + 2.0 * fY * fY) * sqrt(PI) * fHALF * (erf(X - fY) - erf(X + fY)) - \
            (X + fY) * exp(-(X - fY) * (X - fY)) + (X - fY) * exp(-(X + fY) * (X + fY)))
  fSigDiff = isigDiff;
  fAWRI    = aw;
  fTk      = t2 - t1;
  if (t1 == t2) {
    for (int j = 0, x1Size = x1.size(); j != x1Size; ++j) {
      fSigma.push_back(x2[j]);
    }
  } else {
    fALPHA = fAWRI / (fBoltz * fTk);
    fNcrs  = x1.size();
    while (fRATHIG - fRATLOW > 1E-7) {
      fRATHLF = fHALF * (fRATLOW + fRATHIG);
      fHTEST  = XNEPSM(fRATHLF);
      if (fHTEST < isigDiff) {
        fRATLOW = fRATHLF;
      } else {
        fRATHIG = fRATHLF;
      }
    } // end of while loop
    fRATMAX = fRATLOW * fRATLOW;
    if (x1[0] == x1[1]) {
      x1.erase(x1.begin() + 1);
      x2.erase(x2.begin() + 1);
    }
    fIPP               = 0;
    fSize              = x1.size();
    fJLoop             = 0;
    fXss               = 0.5 * (2 * x2[0] + (x2[1] - x2[0]) * (-x1[0]) / (x1[1] - x1[0]));
    if (fXss < 0) fXss = 0;
    for (int k = 0; k < fSize; k++) {
      //         std::cout <<"Loop begins "<< x1[fJLoop] <<"  "<< x2[fJLoop] << std::endl;
      fY2  = x1[fJLoop] * fALPHA;
      fY   = sqrt(fY2);
      fZKT = sqrt(x1[fIPP + 1] * fALPHA) - fY;
      while (fZKT < -fZLIMI) {
        fIPP += 1;
        fZKT = sqrt(x1[fIPP + 1] * fALPHA) - fY;
        if (fZKT + fZLIMI > 0.05) {
          fIPP -= 2;
          fZKT = sqrt(x1[fIPP + 1] * fALPHA) - fY;
          break;
        }
      } // end of while loop
      fMipp   = fIPP;
      fXSUM   = 0.0;
      fFTAIL1 = 0.0;
      fKPP    = fIPP;
      fE2     = x1[fKPP];
      fS2     = x2[fKPP];
      fZK2    = sqrt(fE2 * fALPHA) - fY;
      fZK22   = fZK2 * fZK2;
      fEXPA   = exp(-fZK22);
      fF0K2   = erf(fZK2);
      fF1K2   = fOVSQPI * (fONE - fEXPA);
      fF2K2   = fHALF * fF0K2 - fOVSQPI * fZK2 * fEXPA;
      fF3K2   = fOVSQPI * (fONE - (1 + fZK22) * fEXPA);
      fF4K2   = fHALF * fTHH * fF0K2 - fOVSQPI * fZK2 * (fTHH + fZK22) * fEXPA;

      if (fY < fZLIMI) {
        fZK2P  = fZK2 + 2 * fY;
        fZK22P = fZK2P * fZK2P;
        fEXPAP = exp(-fZK22P);
        fF0K2P = erf(fZK2P);
        fF1K2P = fOVSQPI * (fONE - fEXPAP);
        fF2K2P = fHALF * fF0K2P - fOVSQPI * fZK2P * fEXPAP;
        fF3K2P = fOVSQPI * (fONE - (1 + fZK22P) * fEXPAP);
        fF4K2P = fHALF * fTHH * fF0K2P - fOVSQPI * fZK2P * (fTHH + fZK22P) * fEXPAP;
      }

      while (fZK2 < fZLIMI && fKPP < (int)x1.size() - 1) {
        fE1  = fE2;
        fS1  = fS2;
        fKPP = fKPP + 1;
        fE2  = x1[fKPP];
        fS2  = x2[fKPP];
        while (fE2 == fE1) {
          fKPP = fKPP + 1;
          fE2  = x1[fKPP];
          fS2  = x2[fKPP];
        }
        if (fE2 - fE1 == 0.0 || fE2 == 0 || fALPHA == 0.0) {
          std::cout << "Doppler fails between " << fE1 << "  " << fE2 << "  " << fY2 << "  " << fALPHA << "  " << fAWRI
                    << std::endl;
          continue;
        }
        fF0K1 = fF0K2;
        fF1K1 = fF1K2;
        fF2K1 = fF2K2;
        fF3K1 = fF3K2;
        fF4K1 = fF4K2;
        fZK2  = sqrt(fE2 * fALPHA) - fY;
        fZK22 = fZK2 * fZK2;
        fEXPA = exp(-fZK22);
        fF0K2 = erf(fZK2);
        fF1K2 = fOVSQPI * (fONE - fEXPA);
        fF2K2 = fHALF * fF0K2 - fOVSQPI * fZK2 * fEXPA;
        fF3K2 = fOVSQPI * (fONE - (1 + fZK22) * fEXPA);
        fF4K2 = fHALF * fTHH * fF0K2 - fOVSQPI * fZK2 * (fTHH + fZK22) * fEXPA;
        fFACT = fONE / (fE2 - fE1);
        fAK   = (fE2 * fS1 - fE1 * fS2) * fFACT;
        fCK   = (fS2 - fS1) * fFACT / fALPHA;
        fCKY  = fCK * fY;
        fCKY2 = fCK * fY2;
        fXSUM = fXSUM + fCK * (fF4K2 - fF4K1) + 4 * fCKY * (fF3K2 - fF3K1) + (fAK + 6 * fCKY2) * (fF2K2 - fF2K1) +
                2 * fY * (fAK + 2 * fCKY2) * (fF1K2 - fF1K1) + fY2 * (fAK + fCKY2) * (fF0K2 - fF0K1);

        if (fY < fZLIMI) {
          fZK1P  = fZK2P;
          fF0K1P = fF0K2P;
          fF1K1P = fF1K2P;
          fF2K1P = fF2K2P;
          fF3K1P = fF3K2P;
          fF4K1P = fF4K2P;
          fZK2P  = fZK2 + 2 * fY;
          fZK22P = fZK2P * fZK2P;
          fEXPAP = exp(-fZK22P);
          fF0K2P = erf(fZK2P);
          fF1K2P = fOVSQPI * (fONE - fEXPAP);
          fF2K2P = fHALF * fF0K2P - fOVSQPI * fZK2P * fEXPAP;
          fF3K2P = fOVSQPI * (fONE - (1 + fZK22P) * fEXPAP);
          fF4K2P = fHALF * fTHH * fF0K2P - fOVSQPI * fZK2P * (fTHH + fZK22P) * fEXPAP;
          fXSUM =
              fXSUM - (fCK * (fF4K2P - fF4K1P) - 4 * fCKY * (fF3K2P - fF3K1P) + (fAK + 6 * fCKY2) * (fF2K2P - fF2K1P) -
                       2 * fY * (fAK + 2 * fCKY2) * (fF1K2P - fF1K1P) + fY2 * (fAK + fCKY2) * (fF0K2P - fF0K1P));
        } // end of if
      }   // end of if
      // while loop
      if (fXSUM < 0) fXSUM = 0;
      //      std::cout << fHALF * fXSUM / fY2 <<" xsum1 \t"<< fXss  <<" fFTAIL1 \t"<< 2*x2[0] <<" fFTAIL2 \t"<< x2[
      //      fJLoop ] <<std::endl;
      fFTAIL1 = fXss * (FTAIL(sqrt(x1[0] * fALPHA), fY) - FTAIL(fZERO, fY));
      fXSUM   = fXSUM + fFTAIL1;
      fFTAIL2 = x2[x1.size() - 1] * (FTAIL(sqrt((x1[x1.size() - 1] + 0.1 * x1[x1.size() - 1]) * fALPHA), fY) -
                                     FTAIL(sqrt(x1[x1.size() - 1] * fALPHA), fY));
      fXSUM = fXSUM + fFTAIL2;
      fSigma.push_back(fHALF * fXSUM / fY2);
      //        std::cout << fHALF * fXSUM / fY2 <<" xsum2 \t"<< fHALF * fFTAIL1 / fY2 <<" fFTAIL1 \t"<< fHALF * fFTAIL2
      //        / fY2 <<" fFTAIL2 \t"<< x2[ fJLoop ] <<std::endl;
      if (fJLoop > 0 && k < fSize - 1) {
        fMLoop = 0;
        // 	std::cout <<" before RecursionLinear1 "<< x1[ fJLoop - 1 ] <<"  "<< x1[ fJLoop ] << std::endl;
        RecursionLinear1(x1, x2, x1[fJLoop - 1], x2[fJLoop - 1], fSigma[fJLoop - 1], x1[fJLoop], x2[fJLoop],
                         fSigma[fJLoop]);
        fJLoop += fMLoop;
      }
      fJLoop++;
    } // end of for loop
  }
}

double TNudyEndfDoppler::RecursionLinear1(std::vector<double> &x1, std::vector<double> &x2, double x, double y,
                                          double sig, double xd, double yd, double sigd)
{
#define XNEPSM(S) fabs(fONE - fTW3 * sqrt(fONE / 3.0) * pow((fONE + S * (fONE + S)), fTHH) / (S * (fONE + S)))
#define FTAILX(X, fY)                                                                \
  fOVSQPI *((1.0 + 2.0 * fY * fY) * sqrt(PI) * fHALF * (erf(X - fY) - erf(X + fY)) - \
            (X + fY) * exp(-(X - fY) * (X - fY)) + (X - fY) * exp(-(X + fY) * (X + fY)))
  if (y <= 0.0 && yd <= 0.0) return 0;
  if (fMLoop > 500) return 0;
  double mid     = 0.5 * (x + xd);
  double sigmid1 = y + (yd - y) * (mid - x) / (xd - x);
  double sigmid2 = sig + (sigd - sig) * (mid - x) / (xd - x);
  //   std::cout << sig <<"  "<< sigd <<"  "<< sigmid2 <<"  "<< fMLoop << std::endl;
  //   std::cout << x <<"  "<< xd <<"  "<< mid << std::endl;
  //   std::cout << y <<"  "<< yd <<"  "<< sigmid1 << std::endl;

  std::vector<double>::iterator itx;
  itx        = std::find(x1.begin(), x1.end(), x);
  int xindex = itx - x1.begin();
  x1.insert(x1.begin() + xindex + 1, 1, mid);
  x2.insert(x2.begin() + xindex + 1, 1, sigmid1);
  double sigmid3 = BroadMore(x1, x2, mid);
  double errmid  = fabs(sigmid2 - sigmid3) / sigmid3;
  //  if ( errmid >= fSigDiff ) {
  if (errmid >= fSigDiff) {
    //  if (fabs(sig/sigd -1) > 1E-2 || errmid >= fSigDiff ) {
    //  if (m1 > 0 && m3 > 0 && mdiff > 1E-3  && mdiff != 1) {
    fSigma.insert(fSigma.begin() + xindex + 1, 1, sigmid3);
    fMLoop++;
    //       std::cout << x <<"  \t"<< xd <<" error \t"<< errmid <<"  \t"<<sigmid3<<"  \t"<<sigmid2 << std::endl;
    //       std::cout << sig <<"  \t"<< sigd <<"  \t"<< mid << std::endl;
    //       std::cout << x <<"  \t"<< y <<"  \t"<< xd <<"  \t"<< yd << std::endl;
    RecursionLinear1(x1, x2, x, y, sig, mid, sigmid1, sigmid3);
    RecursionLinear1(x1, x2, mid, sigmid1, sigmid3, xd, yd, sigd);
    return 0;
  } else {
    x1.erase(x1.begin() + xindex + 1);
    x2.erase(x2.begin() + xindex + 1);
    return 0;
  }
  return 0;
}
double TNudyEndfDoppler::BroadMore(std::vector<double> &x1, std::vector<double> &x2, double ixp)
{

  fIPP  = fMipp;
  fKPP  = fIPP;
  fNcrs = x1.size();
  fY2   = ixp * fALPHA;
  fY    = sqrt(fY2);
  fZKT  = sqrt(x1[fIPP + 1] * fALPHA) - fY;
  while (fZKT > -fZLIMI) {
    fIPP -= 1;
    fZKT = sqrt(x1[fIPP - 1] * fALPHA) - fY;
    if (fIPP <= 0) {
      fIPP = 0;
      break;
    }
  } // end of while loop
  fXSUM   = 0.0;
  fFTAIL1 = 0.0;
  fKPP    = fIPP;
  fE2     = x1[fKPP];
  fS2     = x2[fKPP];
  fZK2    = sqrt(fE2 * fALPHA) - fY;
  fZK22   = fZK2 * fZK2;
  fEXPA   = exp(-fZK22);
  fF0K2   = erf(fZK2);
  fF1K2   = fOVSQPI * (fONE - fEXPA);
  fF2K2   = fHALF * fF0K2 - fOVSQPI * fZK2 * fEXPA;
  fF3K2   = fOVSQPI * (fONE - (1 + fZK22) * fEXPA);
  fF4K2   = fHALF * fTHH * fF0K2 - fOVSQPI * fZK2 * (fTHH + fZK22) * fEXPA;
  if (fY < fZLIMI) {
    fZK2P  = fZK2 + 2 * fY;
    fZK22P = fZK2P * fZK2P;
    fEXPAP = exp(-fZK22P);
    fF0K2P = erf(fZK2P);
    fF1K2P = fOVSQPI * (fONE - fEXPAP);
    fF2K2P = fHALF * fF0K2P - fOVSQPI * fZK2P * fEXPAP;
    fF3K2P = fOVSQPI * (fONE - (1 + fZK22P) * fEXPAP);
    fF4K2P = fHALF * fTHH * fF0K2P - fOVSQPI * fZK2P * (fTHH + fZK22P) * fEXPAP;
  }
  while (fZK2 < fZLIMI && fKPP < (int)x1.size() - 1) {
    fE1  = fE2;
    fS1  = fS2;
    fKPP = fKPP + 1;
    fE2  = x1[fKPP];
    fS2  = x2[fKPP];
    while (fE2 == fE1) {
      fKPP = fKPP + 1;
      fE2  = x1[fKPP];
      fS2  = x2[fKPP];
    }
    if (fE2 - fE1 == 0.0 || fY2 == 0.0 || fALPHA == 0.0) {
      std::cout << "Doppler fails " << fE1 << "  " << fE2 << "  " << fY2 << "  " << fALPHA << "  " << fAWRI
                << std::endl;
      continue;
    }
    fF0K1 = fF0K2;
    fF1K1 = fF1K2;
    fF2K1 = fF2K2;
    fF3K1 = fF3K2;
    fF4K1 = fF4K2;
    fZK2  = sqrt(fE2 * fALPHA) - fY;
    fZK22 = fZK2 * fZK2;
    fEXPA = exp(-fZK22);
    fF0K2 = erf(fZK2);
    fF1K2 = fOVSQPI * (fONE - fEXPA);
    fF2K2 = fHALF * fF0K2 - fOVSQPI * fZK2 * fEXPA;
    fF3K2 = fOVSQPI * (fONE - (1 + fZK22) * fEXPA);
    fF4K2 = fHALF * fTHH * fF0K2 - fOVSQPI * fZK2 * (fTHH + fZK22) * fEXPA;
    fFACT = fONE / (fE2 - fE1);
    fAK   = (fE2 * fS1 - fE1 * fS2) * fFACT;
    fCK   = (fS2 - fS1) * fFACT / fALPHA;
    fCKY  = fCK * fY;
    fCKY2 = fCK * fY2;
    fXSUM = fXSUM + fCK * (fF4K2 - fF4K1) + 4 * fCKY * (fF3K2 - fF3K1) + (fAK + 6 * fCKY2) * (fF2K2 - fF2K1) +
            2 * fY * (fAK + 2 * fCKY2) * (fF1K2 - fF1K1) + fY2 * (fAK + fCKY2) * (fF0K2 - fF0K1);
    if (fY < fZLIMI) {
      fZK1P  = fZK2P;
      fF0K1P = fF0K2P;
      fF1K1P = fF1K2P;
      fF2K1P = fF2K2P;
      fF3K1P = fF3K2P;
      fF4K1P = fF4K2P;
      fZK2P  = fZK2 + 2 * fY;
      fZK22P = fZK2P * fZK2P;
      fEXPAP = exp(-fZK22P);
      fF0K2P = erf(fZK2P);
      fF1K2P = fOVSQPI * (fONE - fEXPAP);
      fF2K2P = fHALF * fF0K2P - fOVSQPI * fZK2P * fEXPAP;
      fF3K2P = fOVSQPI * (fONE - (1 + fZK22P) * fEXPAP);
      fF4K2P = fHALF * fTHH * fF0K2P - fOVSQPI * fZK2P * (fTHH + fZK22P) * fEXPAP;
      fXSUM  = fXSUM - (fCK * (fF4K2P - fF4K1P) - 4 * fCKY * (fF3K2P - fF3K1P) + (fAK + 6 * fCKY2) * (fF2K2P - fF2K1P) -
                       2 * fY * (fAK + 2 * fCKY2) * (fF1K2P - fF1K1P) + fY2 * (fAK + fCKY2) * (fF0K2P - fF0K1P));
    } // end of if
  }   // end of if
  // while loop
  if (fXSUM < 0) fXSUM = 0;
  fFTAIL1              = fXss * (FTAILX(sqrt(x1[0] * fALPHA), fY) - FTAILX(fZERO, fY));
  fXSUM                = fXSUM + fFTAIL1;
  fFTAIL2              = x2[x1.size() - 1] * (FTAILX(sqrt((x1[x1.size() - 1] + 0.1 * x1[x1.size() - 1]) * fALPHA), fY) -
                                 FTAILX(sqrt(x1[x1.size() - 1] * fALPHA), fY));
  fXSUM = fXSUM + fFTAIL2;
  return fHALF * fXSUM / fY2;
}
