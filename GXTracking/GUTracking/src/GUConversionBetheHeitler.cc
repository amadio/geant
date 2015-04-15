#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "GUConversionBetheHeitler.h"

#include "backend/Backend.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_BOTH GUConversionBetheHeitler::
GUConversionBetheHeitler(Random_t* states, int threadId) 
  :
  fRandomState(states), fThreadId(threadId),
  fMinX(1.e-8),  fMaxX(1000.), 
  fMaxZelement(100),       // Elements up to Z=100
  fNrow(100), fNcol(100) 
{
  //replace hard coded numbers by default constants
  fAliasSampler = new GUAliasSampler(fRandomState, fThreadId,
                                     fMaxZelement,
                                     fMinX, fMaxX,
                                     fNrow, fNcol);

  for( int z= 1; z < fMaxZelement; ++z)
  {
    //eventually arguments of BuildTable should be replaced by members of *this
    //  and dropped from the function signature. Same for BuildPdfTable
    BuildOneTable(z, fMinX, fMaxX, fNrow, fNcol);
  }
}

VECPHYS_CUDA_HEADER_BOTH GUConversionBetheHeitler::
GUConversionBetheHeitler(Random_t* states, int threadId,
                         GUAliasSampler* sampler) 
  :
  fRandomState(states), fThreadId(threadId),
  fMinX(1.e-8),  fMaxX(1000.), 
  fNrow(100), fNcol(100) 
{
  //replace hard coded numbers by default constants
  fAliasSampler = sampler;
}

//need another Ctor with setable parameters

VECPHYS_CUDA_HEADER_BOTH 
GUConversionBetheHeitler::~GUConversionBetheHeitler() 
{
  if(fAliasSampler) delete fAliasSampler;
}

VECPHYS_CUDA_HEADER_BOTH void 
GUConversionBetheHeitler::BuildOneTable( int Z, 
                                         const double xmin, 
                                         const double xmax,
                                         const int nrow,
                                         const int ncol)
{
  //for now, the model does not own pdf.  Otherwise, pdf should be the 
  //data member of *this and set its point to the fpdf of fAliasSampler 
  double *pdf = new double [nrow*ncol];

  BuildPdfTable(Z,xmin,xmax,nrow,ncol,pdf); 
  fAliasSampler->BuildAliasTable(Z,nrow,ncol,pdf);

  delete [] pdf;
}

VECPHYS_CUDA_HEADER_BOTH void 
GUConversionBetheHeitler::BuildPdfTable(int Z, 
                                        const double xmin, 
                                        const double xmax,
                                        const int nrow,
                                        const int ncol,
                                        double *p)
{
  // Build the probability density function (BetheHeitler pdf) 
  // in the energy randge [xmin,xmax] with an equal bin size
  //
  // input  :  Z    (atomic number) 
  //           xmin (miminum photon energy)
  //           xmax (maxinum photon energy)
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

    double ymin = electron_mass_c2;
    double ymax = x - electron_mass_c2;

    double dy = (ymax - ymin)/(ncol-1); //!!!this should be passed to sampler?
    double yo = ymin + 0.5*dy;
  
    double sum = 0.;

    for(int j = 0; j < ncol ; ++j) {
      //for each output energy bin
      double y = yo + dy*j;
      double xsec = CalculateDiffCrossSection(Z,x,y);
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

// function implementing the cross section for KleinNishina
// TODO: need to get electron properties from somewhere

VECPHYS_CUDA_HEADER_BOTH double 
GUConversionBetheHeitler::CalculateDiffCrossSection( int Zelement, 
                                                     double gammaEnergy, 
                                                     double electEnergy) const
{ 
  // based on Geant4 : G4BetheHeitlerModel
  // input  : gammaEnergy (incomming photon energy)
  //          electEnergy (converted electron/positron energy)
  // output : dsigma  (differential cross section) 

  double epsil  = electEnergy/gammaEnergy;

  double epsil0 = electron_mass_c2/gammaEnergy ;
  if(epsil0 > 1.0) { return 0; }

  // Extract Coulomb factor for this Element
  //F(Z)
  int logZ3 = log(1.0*int(Zelement + 0.5))/3.0; 
  double FZ = 8.*logZ3 ; //8.*(anElement->GetIonisation()->GetlogZ3());
  if (gammaEnergy > 50. /* *MeV */) { 
    FZ += 8.*ComputeCoulombFactor(Zelement); 
  }
  
  //delta -> screenvar
  int Z3 = pow(1.0*int(Zelement + 0.5),1/3.0);
  double screenfac = 136.*epsil0/Z3; //(anElement->GetIonisation()->GetZ3());
  double screenvar = screenfac/(epsil*(1-epsil));
  
  double dsigma = ScreenFunction1(screenvar)*(epsil*epsil+(1.-epsil)*(1.-epsil))                + ScreenFunction2(screenvar)*(2.0/3)*epsil*(1.0-epsil);
  
  return dsigma;
}

VECPHYS_CUDA_HEADER_BOTH double 
GUConversionBetheHeitler::ComputeCoulombFactor(double fZeff) const
{
  // Compute Coulomb correction factor (Phys Rev. D50 3-1 (1994) page 1254)

  const double k1 = 0.0083 , k2 = 0.20206 ,k3 = 0.0020 , k4 = 0.0369 ;
  const double fine_structure_const = (1.0/137); //check unit

  double az1 = fine_structure_const*fZeff;
  double az2 = az1 * az1;
  double az4 = az2 * az2;

  double fCoulomb = (k1*az4 + k2 + 1./(1.+az2))*az2 - (k3*az4 + k4)*az4;
  return fCoulomb;
}


VECPHYS_CUDA_HEADER_BOTH double 
GUConversionBetheHeitler::ScreenFunction1(double screenVariable) const
{
  // compute the value of the screening function 3*PHI1 - PHI2
  double screenVal;
  
  if (screenVariable > 1.)
    screenVal = 42.24 - 8.368*log(screenVariable+0.952);
  else
    screenVal = 42.392 - screenVariable*(7.796 - 1.961*screenVariable);
  
  return screenVal;
}

VECPHYS_CUDA_HEADER_BOTH double 
GUConversionBetheHeitler::ScreenFunction2(double screenVariable) const
{
  // compute the value of the screening function 1.5*PHI1 - 0.5*PHI2
  double screenVal;
  
  if (screenVariable > 1.)
    screenVal = 42.24 - 8.368*log(screenVariable+0.952);
  else
    screenVal = 41.405 - screenVariable*(5.828 - 0.8945*screenVariable);

  return screenVal;
}

} // end namespace impl
} // end namespace vecphys
