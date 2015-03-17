#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "GUMollerBhabha.h"
#include <iostream>

#include "backend/Backend.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_BOTH 
GUMollerBhabha::GUMollerBhabha(Random_t* states, int threadId) 
  :
  fRandomState(states), fThreadId(threadId),
  fMinX(1.), fMaxX(1001.), fDeltaX(0.1), 
  fMinY(0.), fMaxY(1001.), fDeltaY(0.1),
  fNrow(100), fNcol(100) 
{
  //replace hard coded numbers by default constants

  fAliasSampler = new GUAliasSampler(fRandomState,fThreadId,
                                     10,fMinX,fMaxX,fNrow,fNcol);

  //eventually arguments of BuildTable should be replaced by members of *this
  //and dropped from the function signature. Same for BuildPdfTable
  BuildTable(10,fMinX,fMaxX,fNrow,fNcol);   
}

VECPHYS_CUDA_HEADER_BOTH 
GUMollerBhabha::GUMollerBhabha(Random_t* states, int threadId,
                                             GUAliasSampler* sampler) 
  :
  fRandomState(states), fThreadId(threadId),
  fMinX(1.), fMaxX(1001.), fDeltaX(0.1), 
  fMinY(0.), fMaxY(1001.), fDeltaY(0.1),
  fNrow(100), fNcol(100) 
{
  //replace hard coded numbers by default constants
  fAliasSampler = sampler;
}

//need another Ctor with setable parameters

VECPHYS_CUDA_HEADER_BOTH 
GUMollerBhabha::~GUMollerBhabha() 
{
  if(fAliasSampler) delete fAliasSampler;
}


VECPHYS_CUDA_HEADER_BOTH void 
GUMollerBhabha::BuildTable( int Z, 
                                   const double xmin, 
                                   const double xmax,
                                   const int nrow,
                                   const int ncol)
{
  //for now, the model does not own pdf.  Otherwise, pdf should be the 
  //data member of *this and set its point to the fpdf of fAliasSampler 
  double *pdf = new double [nrow*ncol];

  BuildPdfTable(Z,xmin,xmax,nrow,ncol,pdf); 
  fAliasSampler->BuildAliasTables(nrow,ncol,pdf);

  delete [] pdf;
}

VECPHYS_CUDA_HEADER_BOTH void 
GUMollerBhabha::BuildPdfTable(int Z, 
                                     const double xmin, 
                                     const double xmax,
                                     const int nrow,
                                     const int ncol,
                                     double *p
                                     )
{
  // Build the probability density function (MollerBhabha pdf) 
  // in the energy randge [xmin,xmax] with an equal bin size
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  //           xmin (miminum energy)
  //           xmax (maxinum energy)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  //build pdf  
  double dx = (xmax - xmin)/nrow;
  double xo =  xmin + 0.5*dx;


  for(int i = 0; i < nrow ; ++i) {
    //for each input energy bin
    double x = xo + dx*i;

    double ymin = x/(1+2.0*x*inv_electron_mass_c2);
    double dy = (x - ymin)/(ncol-1);
    double yo = ymin + 0.5*dy;
  
    double sum = 0.;

    for(int j = 0; j < ncol ; ++j) {
      //for each output energy bin
      double y = yo + dy*j;
      double xsec = CalculateDiffCrossSection(0,x,y);
      p[i*ncol+j] = xsec;
      sum += xsec;
    }

    //normalization
    sum = 1.0/sum;

    for(int j = 0; j < ncol ; ++j) {
      p[i*ncol+j] *= sum;
    }
  }
}

// function implementing the cross section for MollerBhabha

VECPHYS_CUDA_HEADER_BOTH double 
GUMollerBhabha::CalculateDiffCrossSection( int Zelement, //dummy for now
                                           double kineticEnergy, 
					   double deltaRayEnergy) const
{
  // based on Geant3 : Simulation of the delta-ray production (PHY331-1)
  // input  : kineticEnergy (incomming photon energy)
  //          deltaRayEnergy (scattered photon energy)
  // output : dcross (differential cross section) 

  double dcross = 0.0;

  double tau = kineticEnergy/electron_mass_c2;
  double gam = tau + 1.0;
  double gamma2 = gam*gam;
  double epsil = deltaRayEnergy/kineticEnergy;

  //Moller (e-e-) scattering only
  //Bhabha (e+e-) scattering not implemented for now

  double fgam = (2*gam-1.0)/gamma2;
  double x = 1/epsil;
  double y = 1/(1.0-epsil);

  dcross = (1-fgam + x*(x-fgam) + y*y*fgam);

  return dcross;
}

} // end namespace impl
} // end namespace vecphys
