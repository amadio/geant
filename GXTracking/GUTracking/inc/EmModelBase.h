#ifndef EmModelBase_H
#define EmModelBase_H 1

#include "backend/Backend.h"
#include "base/SystemOfUnits.h"

#include "GUConstants.h"
#include "GUTrack.h"
#include "GUAliasSampler.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

template <class EmModel>
class EmModelBase {

public:

  VECPHYS_CUDA_HEADER_HOST
  EmModelBase(Random_t* states, int tid);

  VECPHYS_CUDA_HEADER_BOTH 
  EmModelBase(Random_t* states, int tid, GUAliasSampler* sampler);

  VECPHYS_CUDA_HEADER_BOTH 
  ~EmModelBase();

  VECPHYS_CUDA_HEADER_HOST
  void BuildCrossSectionTable();

  VECPHYS_CUDA_HEADER_HOST
  void BuildAliasTable();

  //scalar
  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH 
  void AtomicCrossSection(GUTrack&  projectile,   
                          const int targetElement,
                          double&   sigma);

  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH 
  void Interact(GUTrack&  projectile,   
                const int targetElement,
                GUTrack&  secondary );

  //vector
#ifndef VECPHYS_NVCC
  template <typename Backend>
  void AtomicCrossSection(GUTrack_v& inProjectile,  
                          const int* targetElements,
                          double*    sigma);     

  template <typename Backend>
  void Interact(GUTrack_v& inProjectile,  
                const int* targetElements,
                GUTrack_v& outSecondaryV);     
#endif

  //validation 
  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void AtomicCrossSectionG4(GUTrack&  inProjectile,
                            const int targetElement,
                            double&   sigma);

  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void InteractG4(GUTrack&  inProjectile,
                  const int targetElement,
                  GUTrack&  outSecondary);

  VECPHYS_CUDA_HEADER_BOTH
  GUAliasSampler* GetSampler() {return fAliasSampler;}

  VECPHYS_CUDA_HEADER_BOTH
  void SetSampler(GUAliasSampler* sampler) { fAliasSampler = sampler ;}


protected:
  // Auxiliary methods
  VECPHYS_CUDA_HEADER_BOTH
  void SetLowEnergyLimit(double lowLimit) { fLowEnergyLimit = lowLimit; }

  VECPHYS_CUDA_HEADER_BOTH
  void SetHighEnergyLimit(double highLimit) { fHighEnergyLimit = highLimit; }

  VECPHYS_CUDA_HEADER_BOTH double 
  ComputeCoulombFactor(double fZeff);

private:
  // Implementation methods
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void RotateAngle(typename Backend::double sinTheta,
                   typename Backend::double xhat,
                   typename Backend::double yhat,
                   typename Backend::double zhat,
                   typename Backend::double &xr,
                   typename Backend::double &yr,
                   typename Backend::double &zr);

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void ConvertXtoFinalState(double energyIn, 
                            double energyOut, 
                            double sinTheta, 
                            GUTrack& primary, 
                            GUTrack& secondary);

  //data members
protected:
  Random_t* fRandomState;
  int       fThreadId;

  double    fLowEnergyLimit;  
  double    fHighEnergyLimit;  

  //Sampling Tables
  GUAliasSampler* fAliasSampler; 
  Precision fMinX; // lower bound of the sampling variable
  Precision fMaxX; // upper bound of the sampling variable
  int       fNrow; // number of input energy bins of the alias table
  int       fNcol; // number of output energy bins of the alias table
  Precision fMaxZelement; 
};

//Implementation
template <class EmModel>
VECPHYS_CUDA_HEADER_HOST 
EmModelBase<EmModel>::EmModelBase(Random_t* states, int tid) 
  : fRandomState(states), fThreadId(tid), 
    fLowEnergyLimit(0.1*keV), fHighEnergyLimit(10.*TeV),
    fMinX(1.e-8),  fMaxX(1000.),  
    fNrow(100), fNcol(100),
    fMaxZelement(maximumZ)
{
  fAliasSampler = 0;
}

template <class EmModel>
VECPHYS_CUDA_HEADER_BOTH 
EmModelBase<EmModel>::EmModelBase(Random_t* states, int tid, GUAliasSampler* sampler) 
  : fRandomState(states), fThreadId(tid), 
    fLowEnergyLimit(0.1*keV), fHighEnergyLimit(10.*TeV),
    fMinX(1.e-8),  fMaxX(1000.), 
    fNrow(100), fNcol(100),
    fMaxZelement(maximumZ)
{
  fAliasSampler = sampler; 
}

template <class EmModel>
VECPHYS_CUDA_HEADER_BOTH 
EmModelBase<EmModel>::~EmModelBase()
{
  if(fAliasSampler) delete fAliasSampler;
}

template <class EmModel>
VECPHYS_CUDA_HEADER_HOST 
void EmModelBase<EmModel>::BuildCrossSectionTable() 
{ 
  //dummy interface for now
  for( int z= 1 ; z < fMaxZelement; ++z)
  {
    static_cast<EmModel*>(this)->BuildCrossSectionTablePerAtom(z); 
  }
}

template <class EmModel>
VECPHYS_CUDA_HEADER_HOST 
void EmModelBase<EmModel>::BuildAliasTable() 
{ 
  //replace hard coded numbers by default constants
  fAliasSampler = new GUAliasSampler(fRandomState, fThreadId, fMaxZelement,
                                     fMinX, fMaxX, fNrow, fNcol);

  //temporary array
  double *pdf = new double [(fNrow+1)*fNcol];
  for( int z= 1 ; z < fMaxZelement; ++z)
  {
    static_cast<EmModel*>(this)->BuildPdfTable(z,pdf); 
    fAliasSampler->BuildAliasTable(z,fNrow,fNcol,pdf);
  }
  delete [] pdf;
}

template <class EmModel>
template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void EmModelBase<EmModel>::AtomicCrossSection(GUTrack&  inProjectile,
                                              const int targetElement,
                                              double&   sigma ) 
{
  sigma = 0.;
  double energyIn= inProjectile.E;
  if(energyIn > fLowEnergyLimit) {
    sigma = static_cast<EmModel*>(this)-> template CrossSectionKernel<Backend>(energyIn,targetElement);
  }
}

template <class EmModel>
template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void EmModelBase<EmModel>::Interact(GUTrack&  inProjectile,
                                    const int targetElement,
                                    GUTrack&  outSecondary ) 
{
  double energyIn= inProjectile.E;
  double energyOut, sinTheta;

  static_cast<EmModel*>(this)-> template InteractKernel<Backend>(energyIn,targetElement,energyOut,sinTheta);

  //update final states of the primary and store the secondary
  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
}
  
#ifndef VECPHYS_NVCC
template <class EmModel>
template <typename Backend>
void EmModelBase<EmModel>::AtomicCrossSection(GUTrack_v& inProjectile,
                                              const int* targetElements,
                                              double*    sigma) 
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::double double;

  for(int j = 0; j < inProjectile.numTracks  ; ++j) {
    assert( (targetElements[j] > 0)  && (targetElements[j] <= maximumZ ) );
  }

  int ibase= 0;
  int numChunks= (inProjectile.numTracks/double::Size);

  for(int i=0; i < numChunks ; ++i) {
    double energyIn(&inProjectile.E[ibase]);
    Index_t  zElement(targetElements[ibase]);

    double sigmaOut = static_cast<EmModel*>(this)-> template CrossSectionKernel<Backend>(energyIn,zElement);

    sigmaOut.store(&sigma[ibase]);
    ibase+= double::Size;
  }

  //leftover - do scalar
  for(int i = numChunks*double::Size ; i < inProjectile.numTracks ; ++i) {
    sigma[i] = static_cast<EmModel*>(this)-> template CrossSectionKernel<kScalar>(inProjectile.E[i],targetElements[i]);
  }
}

template <class EmModel>
template <typename Backend>
void EmModelBase<EmModel>::Interact(GUTrack_v& inProjectile,  
                                    const int* targetElements,
                                    GUTrack_v& outSecondary) 
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::double double;

  for(int j = 0; j < inProjectile.numTracks  ; ++j) {
    assert( (targetElements[j] > 0)  && (targetElements[j] <= maximumZ ) );
  }

  int ibase= 0;
  int numChunks= (inProjectile.numTracks/double::Size);

  for(int i= 0; i < numChunks ; ++i) {
    double energyIn(&inProjectile.E[ibase]);
    double px(&inProjectile.px[ibase]);
    double py(&inProjectile.py[ibase]);
    double pz(&inProjectile.pz[ibase]);
    double sinTheta;
    double energyOut;
    Index_t  zElement(targetElements[ibase]);

    static_cast<EmModel*>(this)-> template InteractKernel<Backend>(energyIn,zElement,energyOut,sinTheta);

    //need to rotate the angle with respect to the line of flight
    double invp = 1./energyIn;
    double xhat = px*invp;
    double yhat = py*invp;
    double zhat = pz*invp;

    double uhat = 0.;
    double vhat = 0.;
    double what = 0.;

    RotateAngle<Backend>(sinTheta,xhat,yhat,zhat,uhat,vhat,what);

    // Update primary
    energyOut.store(&inProjectile.E[ibase]);
    double pxFinal, pyFinal, pzFinal;
     
    pxFinal= energyOut*uhat;
    pyFinal= energyOut*vhat;
    pzFinal= energyOut*what;
    pxFinal.store(&inProjectile.px[ibase]);
    pyFinal.store(&inProjectile.py[ibase]);
    pzFinal.store(&inProjectile.pz[ibase]);

    // create Secondary
    double secE = energyIn - energyOut; 
    double pxSec= secE*(xhat-uhat);
    double pySec= secE*(yhat-vhat);
    double pzSec= secE*(zhat-what);

    secE.store(&outSecondary.E[ibase]);
    pxSec.store(&outSecondary.px[ibase]);
    pySec.store(&outSecondary.py[ibase]);
    pzSec.store(&outSecondary.pz[ibase]);

    ibase+= double::Size;
  }

  //leftover - do scalar (temporary)
  for(int i = numChunks*double::Size ; i < inProjectile.numTracks ; ++i) {

    double senergyIn= inProjectile.E[i];
    double senergyOut, ssinTheta;

    static_cast<EmModel*>(this)-> template InteractKernel<kScalar>(senergyIn,targetElements[i],senergyOut,ssinTheta);

    //need to rotate the angle with respect to the line of flight
    double sinvp = 1./senergyIn;
    double sxhat = inProjectile.px[i]*sinvp;
    double syhat = inProjectile.py[i]*sinvp;
    double szhat = inProjectile.pz[i]*sinvp;

    double suhat = 0.;
    double svhat = 0.;
    double swhat = 0.;

    RotateAngle<kScalar>(ssinTheta,sxhat,syhat,szhat,suhat,svhat,swhat);

    //update primary
    inProjectile.E[i]  = senergyOut;
    inProjectile.px[i] = senergyOut*suhat;
    inProjectile.py[i] = senergyOut*svhat;
    inProjectile.pz[i] = senergyOut*swhat;

    //create secondary
    outSecondary.E[i]  = (senergyIn-senergyOut); 
    outSecondary.px[i] = outSecondary.E[i]*(sxhat-suhat);
    outSecondary.py[i] = outSecondary.E[i]*(syhat-svhat);
    outSecondary.pz[i] = outSecondary.E[i]*(szhat-swhat);
    //fill other information
  }
}
#endif

template <class EmModel>
template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void EmModelBase<EmModel>::AtomicCrossSectionG4(GUTrack&  inProjectile,
                                                const int targetElement,
                                                double&   sigma)
{
  sigma = 0.;
  double energyIn = inProjectile.E;

  if(energyIn > fLowEnergyLimit) {
    sigma = static_cast<EmModel*>(this)->GetG4CrossSection(energyIn,targetElement);
  }
}

template <class EmModel>
template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void EmModelBase<EmModel>::InteractG4(GUTrack&  inProjectile,
                                      const int targetElement,
                                      GUTrack&  outSecondary)
{

  Precision energyIn = inProjectile.E;
  Precision energyOut;
  Precision sinTheta;

  static_cast<EmModel*>(this)->SampleByCompositionRejection(targetElement,energyIn,energyOut,sinTheta);

  //update final states of the primary and store the secondary
  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
}

template <class EmModel>
template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void
EmModelBase<EmModel>::RotateAngle(typename Backend::double sinTheta,
                                  typename Backend::double xhat,
                                  typename Backend::double yhat,
                                  typename Backend::double zhat,
                                  typename Backend::double &xr,
                                  typename Backend::double &yr,
                                  typename Backend::double &zr)
{
  typedef typename Backend::int    int;
  typedef typename Backend::double double;
  typedef typename Backend::bool   bool;

  double phi = UniformRandom<Backend>(fRandomState,int(fThreadId));
  double pt = xhat*xhat + yhat*yhat;

  double cosphi, sinphi;
  sincos(phi, &sinphi, &cosphi);

  double uhat = sinTheta*cosphi; // cos(phi);
  double vhat = sinTheta*sinphi; // sin(phi);
  double what = Sqrt((1.-sinTheta)*(1.+sinTheta));

  bool positive = ( pt > 0. );
  bool negativeZ = ( zhat < 0. );

  //mask operation???
  if(positive) {
    double phat = Sqrt(pt);
    xr = (xhat*zhat*uhat - yhat*vhat)/phat + xhat*what;
    yr = (yhat*zhat*uhat - xhat*vhat)/phat + yhat*what;
    zr = -phat*uhat + zhat*what;
  }
  else if(negativeZ) {
    xr = -xhat;
    yr =  yhat;
    zr = -zhat;
  }  
  else {
    xr = xhat;
    yr = yhat;
    zr = zhat;
  }
}

template <class EmModel>
template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void EmModelBase<EmModel>::ConvertXtoFinalState(double energyIn, 
                                                double energyOut, 
                                                double sinTheta, 
                                                GUTrack& inProjectile,
                                                GUTrack& outSecondary )
{
  //need to rotate the angle with respect to the line of flight
  double invp = 1./energyIn;
  double xhat = inProjectile.px*invp;
  double yhat = inProjectile.py*invp;
  double zhat = inProjectile.pz*invp;

  double uhat = 0.;
  double vhat = 0.;
  double what = 0.;

  RotateAngle<Backend>(sinTheta,xhat,yhat,zhat,uhat,vhat,what);

  //update primary
  inProjectile.E  = energyOut;
  inProjectile.px = energyOut*uhat;
  inProjectile.py = energyOut*vhat;
  inProjectile.pz = energyOut*what;

  //create secondary
  outSecondary.E  = (energyIn-energyOut); 
  outSecondary.px = outSecondary.E*(xhat-uhat);
  outSecondary.py = outSecondary.E*(yhat-vhat);
  outSecondary.pz = outSecondary.E*(zhat-what);
  //fill other information
}

template <class EmModel>
VECPHYS_CUDA_HEADER_BOTH double 
EmModelBase<EmModel>::ComputeCoulombFactor(double Zeff)
{
  // Compute Coulomb correction factor (Phys Rev. D50 3-1 (1994) page 1254)

  const double k1 = 0.0083 , k2 = 0.20206 ,k3 = 0.0020 , k4 = 0.0369 ;
  const double fine_structure_const = (1.0/137); //check unit

  double az1 = fine_structure_const*Zeff;
  double az2 = az1 * az1;
  double az4 = az2 * az2;

  double coulombFactor = (k1*az4 + k2 + 1./(1.+az2))*az2 - (k3*az4 + k4)*az4;
  return coulombFactor;
}

} // end namespace impl
} // end namespace vecphys

#endif
