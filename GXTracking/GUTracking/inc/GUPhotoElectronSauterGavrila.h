#ifndef GUPhotoElectronSauterGavrila_H
#define GUPhotoElectronSauterGavrila_H 1

#include "backend/Backend.h"

#include "GUConstants.h"
#include "GUTrack.h"
#include "PhysicalConstants.h"
#include "StaticSandiaData.h"

// add the sincos function on MAC because sincos is not part of math.h
#ifdef __APPLE__ // possibly other conditions
inline void sincos(double x, double *s, double *c){
  __sincos(x,s,c);
}
#endif

namespace vecphys {

  //VECPHYS_DEVICE_DECLARE_CONV( GUPhotoElectronSauterGavrila )
VECPHYS_DEVICE_FORWARD_DECLARE( class GUAliasSampler; )

inline namespace VECPHYS_IMPL_NAMESPACE {

class GUPhotoElectronSauterGavrila
{
public:

  VECPHYS_CUDA_HEADER_HOST
  GUPhotoElectronSauterGavrila(Random_t* states, int threadId = -1); 

  VECPHYS_CUDA_HEADER_BOTH
  GUPhotoElectronSauterGavrila(Random_t* states, int threadId, 
                        GUAliasSampler* sampler); 

  VECPHYS_CUDA_HEADER_BOTH
  ~GUPhotoElectronSauterGavrila();

  // VECPHYS_CUDA_HEADER_BOTH
  // void GetSampleParameters(double x, int &irow, int &icol,double &t);
  //  --> Apparent method above is neither defined, nor used ... 
  
  // Core Interface methods
  // -------------------------------------------
  
  // Generate secondaries 
  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void Interact(GUTrack& projectile,    // In/Out: Updated to new state - choice
                int      targetElement, // Q: Need Material index instead ? 
                GUTrack& secondary ) const;

  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void InteractG4(GUTrack& inProjectile,
                  int      targetElement,
                  GUTrack& outSecondary );

  // Vector version - stage 2 - Need the definition of these vector types
#ifndef VECPHYS_NVCC
  template <typename Backend>
  void Interact(GUTrack_v& inProjectile,    // In/Out
                const int *targetElements,  // Number equal to num of tracks
                GUTrack_v& outSecondaryV    // Empty vector for secondaries
                )const;     
#endif

  // Core method(s)
  // -------------------------------------------  
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH void 
  InteractKernel(typename Backend::Double_t energyIn, 
                 typename Backend::Index_t   zElement,
                 typename Backend::Double_t& energyOut,
                 typename Backend::Double_t& sinTheta) const;


  // Alternative Implementation method(s) - for reference/comparisons
  // ----------------------------------------------------------------
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void SampleByCompositionRejection(typename Backend::Double_t energyIn,
                                    typename Backend::Double_t& energyOut,
                                    typename Backend::Double_t& sinTheta);

  //  Initialisation methods
  // -------------------------------------------

  // Initializes this class and its sampler 
  VECPHYS_CUDA_HEADER_HOST
  void BuildOneTable( int Z,
                      const double xmin,
                      const double xmax,
                      const int nrow,
                      const int ncol);

  VECPHYS_CUDA_HEADER_HOST
  void BuildPdfTable(int Z,
                     const double xmin,
                     const double xmax,
                     const int nrow,
                     const int ncol,
                     double *p);

  VECPHYS_CUDA_HEADER_HOST
  void BuildLogPdfTable(int Z,
                        const double xmin,
                        const double xmax,
                        const int nrow,
                        const int ncol,
                        double *p);
  
public:
  // Auxiliary methods
  VECPHYS_CUDA_HEADER_BOTH
  GUAliasSampler* GetSampler() {return fAliasSampler;}

  VECPHYS_CUDA_HEADER_BOTH
  void SetSampler(GUAliasSampler* sampler) { fAliasSampler = sampler ;}

public:
  // Implementation methods - used to implement Interact
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void
  RotateAngle(typename Backend::Double_t sinTheta,
              typename Backend::Double_t xhat,
              typename Backend::Double_t yhat,
              typename Backend::Double_t zhat,
              typename Backend::Double_t &xr,
              typename Backend::Double_t &yr,
              typename Backend::Double_t &zr) const;

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  SampleSinTheta(typename Backend::Double_t energyIn,
                 typename Backend::Double_t energyOut) const; 

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t 
  TotalCrossSection(typename Backend::Double_t energy,
                    typename Backend::Int_t zElement) const;

  template<class Backend>
  inline
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t 
  GetPhotoElectronEnergy(typename Backend::Double_t energyIn,
                         typename Backend::Int_t zElement) const;

private: 
  // Implementation methods 

  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSectionK( int Zelement, double Ein, double outEphoton ) const;

  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton ) const;

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void ConvertXtoFinalState(double energyIn, 
                            double energyOut, 
                            double sinTheta, 
                            GUTrack& primary, 
                            GUTrack& secondary) const;

private:
  GUAliasSampler* fAliasSampler; 

  // Helper data members for GPU random -- to be replaced by use of a GPU manager class
  Random_t* fRandomState;
  int       fThreadId;

  Precision fMinX;   // E Minimum - lowest energy for projectile
  Precision fMaxX;
  // Precision fDeltaX;

  // Precision fMinY, fMaxY, fDeltaY;  // Energy limits for outgoing particles ? Not used

  Precision fMaxZelement; // 
  
  //Sampling Tables
  int fNrow;
  int fNcol;
};

//Implementation

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUPhotoElectronSauterGavrila::
InteractKernel(typename Backend::Double_t energyIn, 
               typename Backend::Index_t   zElement,
               typename Backend::Double_t& energyOut,
               typename Backend::Double_t& sinTheta) const
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::Double_t Double_t;

  //energy of photo-electron: Sandia parameterization
  energyOut = GetPhotoElecticEnergy(energyIn,zElement) ;

  //sample angular distribution of photo-electron

  Index_t   index;
  Index_t   icol;
  Double_t  fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn,index,icol,fraction);

  Double_t probNA;
  Double_t aliasInd;

  //this did not used to work - Fixed SW
  fAliasSampler->GatherAlias<Backend>(index,zElement,probNA,aliasInd);
  
  Double_t mininum = -1.0;
  Double_t deltaE = 2.0;

  Double_t  costTheta = mininum + fAliasSampler->SampleX<Backend>(deltaE,probNA,
					        aliasInd,icol,fraction);
  sinTheta = sqrt((1+costTheta)*(1-costTheta));
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
void
GUPhotoElectronSauterGavrila::RotateAngle(typename Backend::Double_t sinTheta,
                                          typename Backend::Double_t xhat,
                                          typename Backend::Double_t yhat,
                                          typename Backend::Double_t zhat,
                                          typename Backend::Double_t &xr,
                                          typename Backend::Double_t &yr,
                                          typename Backend::Double_t &zr) const
{
  typedef typename Backend::Double_t Double_t;
  typedef typename Backend::Bool_t   Bool_t;

  Double_t phi = UniformRandom<Backend>(fRandomState,fThreadId);

  Double_t pt = xhat*xhat + yhat*yhat;

  Double_t cosphi, sinphi;
  sincos(phi, &sinphi, &cosphi);

  Double_t uhat = sinTheta*cosphi; // cos(phi);
  Double_t vhat = sinTheta*sinphi; // sin(phi);
  Double_t what = Sqrt((1.-sinTheta)*(1.+sinTheta));

  Bool_t positive = ( pt > 0. );
  Bool_t negativeZ = ( zhat < 0. );

  //mask operation???
  if(positive) {
    Double_t phat = Sqrt(pt);
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

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t 
GUPhotoElectronSauterGavrila::
SampleSinTheta(typename Backend::Double_t energyIn,
               typename Backend::Double_t energyOut) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  //angle of the scatterred photon

  Double_t epsilon = energyOut/energyIn;

  Bool_t condition = epsilon > 1.0;

  MaskedAssign( condition, 1.0 , &epsilon );

  Double_t E0_m    = inv_electron_mass_c2*energyIn;
  Double_t onecost = (1.0 - epsilon)/(epsilon*E0_m);
  Double_t sint2   = onecost*(2.-onecost);

  Double_t sinTheta = 0.5;
  Bool_t condition2 = sint2 < 0.0;

  MaskedAssign(  condition2, 0.0, &sinTheta );   // Set sinTheta = 0
  MaskedAssign( !condition2, Sqrt(sint2), &sinTheta );   
  
  return sinTheta;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUPhotoElectronSauterGavrila::
SampleByCompositionRejection(typename Backend::Double_t  energyIn,
                             typename Backend::Double_t& energyOut,
                             typename Backend::Double_t& cosTheta)
{
  typedef typename Backend::Double_t Double_t;

  //use the scalar implementation which is equivalent to Geant4
  energyOut = GetPhotoElectronEnergy<kScalar>(energyIn,10);

  //sample angular direction according to SauterGavrilaAngularDistribution

  Double_t tau = energyIn/electron_mass_c2;

  Double_t gamma     = tau + 1.0;
  Double_t beta      = sqrt(tau*(tau + 2.0))/gamma;

  Double_t A = (1-beta)/beta;
  Double_t Ap2 = A + 2;
  Double_t B   = 0.5*beta*gamma*(gamma-1.)*(gamma-2.);
  Double_t grej = 2*(1+A*B)/A;
  
  Double_t z;
  Double_t g;

  do { 
    Double_t q = UniformRandom<Backend>(fRandomState,fThreadId);
    z = 2*A*(2*q + Ap2*sqrt(q))/(Ap2*Ap2 - 4*q);
    g = (2 - z)*(1.0/(A + z) + B);
  } while(g < UniformRandom<Backend>(fRandomState,fThreadId)*grej);
  
  cosTheta = 1 - z;
}

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUPhotoElectronSauterGavrila::Interact(GUTrack& inProjectile,
                                     int      targetElement,
                                     GUTrack& outSecondary ) const
{
  double energyIn= inProjectile.E;
  double energyOut, sinTheta;
#ifdef CHECK
  if( (energyIn <= fMinX) || (energyIn > fMaxX) )
  {
    printf(" Illegal input Energy = %f min = %f max = %f\n",
	   energyIn,fMinX,fMaxX);
  }
#endif 
  //  assert( (energyIn >= fMinX)  && (energyIn <= fMaxX) );
  InteractKernel<Backend>(energyIn, targetElement, energyOut, sinTheta);
  
  //update final states of the primary and store the secondary
  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
}
  
#ifndef VECPHYS_NVCC
template <typename Backend>
void GUPhotoElectronSauterGavrila::Interact( GUTrack_v& inProjectile,    // In/Out
                                      const int *targetElements,  // Number equal to num of tracks
                                      GUTrack_v& outSecondary    // Empty vector for secondaries
                                      ) const
{
  typedef typename Backend::Double_t Double_t;
  typedef typename Backend::Index_t  Index_t;

  for(int j = 0; j < inProjectile.numTracks  ; ++j) {
     assert( (targetElements[j] > 0)  && (targetElements[j] <= fMaxZelement) );
  }
  
  int ibase= 0;
  int numChunks= (inProjectile.numTracks/Double_t::Size);

  for(int i=0; i < numChunks ; ++i) {
    Double_t energyIn(inProjectile.E[ibase]);
    Double_t px(inProjectile.px[ibase]);
    Double_t py(inProjectile.py[ibase]);
    Double_t pz(inProjectile.pz[ibase]);
    Double_t sinTheta;
    Double_t energyOut;
    Index_t  zElement(targetElements[ibase]);

    InteractKernel<Backend>(energyIn, zElement, energyOut, sinTheta);

    //need to rotate the angle with respect to the line of flight
    Double_t invp = 1./energyIn;
    Double_t xhat = px*invp;
    Double_t yhat = py*invp;
    Double_t zhat = pz*invp;

    Double_t uhat = 0.;
    Double_t vhat = 0.;
    Double_t what = 0.;

    RotateAngle<Backend>(sinTheta,xhat,yhat,zhat,uhat,vhat,what);

    // Update primary
    energyOut.store(&inProjectile.E[ibase]);
    Double_t pxFinal, pyFinal, pzFinal;
     
    pxFinal= energyOut*uhat;
    pyFinal= energyOut*vhat;
    pzFinal= energyOut*what;
    pxFinal.store(&inProjectile.px[ibase]);
    pyFinal.store(&inProjectile.py[ibase]);
    pzFinal.store(&inProjectile.pz[ibase]);

    // create Secondary
    Double_t secE = energyIn - energyOut; 
    Double_t pxSec= secE*(xhat-uhat);
    Double_t pySec= secE*(yhat-vhat);
    Double_t pzSec= secE*(zhat-what);

    secE.store(&outSecondary.E[ibase]);
    pxSec.store(&outSecondary.px[ibase]);
    pySec.store(&outSecondary.py[ibase]);
    pzSec.store(&outSecondary.pz[ibase]);

    ibase+= Double_t::Size;
  }
}    
#endif

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUPhotoElectronSauterGavrila::ConvertXtoFinalState(double energyIn, 
                                                 double energyOut, 
                                                 double sinTheta, 
                                                 GUTrack& inProjectile,
                                                 GUTrack& outSecondary ) const
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

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUPhotoElectronSauterGavrila::InteractG4(GUTrack& inProjectile,
                                       int      targetElement,
                                       GUTrack& outSecondary )
{
  Precision energyIn;

  energyIn = inProjectile.E;

  Precision energyOut;
  Precision sinTheta;
  SampleByCompositionRejection<Backend>(energyIn,energyOut,sinTheta);

  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
  
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t 
GUPhotoElectronSauterGavrila::
TotalCrossSection(typename Backend::Double_t energy,
                  typename Backend::Int_t  zElement) const
{
  typedef typename Backend::Double_t Double_t;

  Double_t sigma = 0.;

  //Sandia parameterization for Z < 100
  int Z = zElement;

  int    fCumulInterval[101]  = {0};
  double fSandiaCof[4]        = {0.0};

  fCumulInterval[0] = 1;

  //scan
  for (int iz = 1; iz < 101; ++iz) {     
    fCumulInterval[iz] = fCumulInterval[iz-1] + fNbOfIntervals[iz];
  }

  double Emin  = fSandiaTable[fCumulInterval[Z-1]][0]*keV;

  int interval = fNbOfIntervals[Z] - 1;
  int row = fCumulInterval[Z-1] + interval;

  while ((interval>0) && (energy<fSandiaTable[row][0]*keV)) {
    --interval;
    row = fCumulInterval[Z-1] + interval;
  }

  if (energy >= Emin) {        
    double AoverAvo = Z*amu/fZtoAratio[Z];
    fSandiaCof[0]=AoverAvo*funitc[1]*fSandiaTable[row][1];     
    fSandiaCof[1]=AoverAvo*funitc[2]*fSandiaTable[row][2];     
    fSandiaCof[2]=AoverAvo*funitc[3]*fSandiaTable[row][3];     
    fSandiaCof[3]=AoverAvo*funitc[4]*fSandiaTable[row][4];
  }
  else {
    fSandiaCof[0] = fSandiaCof[1] = fSandiaCof[2] = fSandiaCof[3] = 0.;
  }     

  Double_t energy2 = energy*energy;
  Double_t energy3 = energy*energy2;
  Double_t energy4 = energy2*energy2;

  Double_t sgima = fSandiaCof[0]/energy  + fSandiaCof[1]/energy2 +
    fSandiaCof[2]/energy3 + fSandiaCof[3]/energy4;

  return sigma;
}

template<class Backend>
inline
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
GUPhotoElectronSauterGavrila::
GetPhotoElectronEnergy(typename Backend::Double_t energy,
                       typename Backend::Int_t  zElement) const
{
  // this method is not vectorizable and only for the scalar backend

  typedef typename Backend::Int_t Int_t;
  typedef typename Backend::Double_t Double_t;

  // Photo electron energy
  Double_t energyOut = 0.;

  // Select atomic shell
  assert (zElement>0 && zElement <101);
  Int_t nShells = fNumberOfShells[zElement];

  Int_t i = 0;  
  Double_t bindingEnergy =0;

  for( ; i < nShells ; ++i) {
    bindingEnergy = fBindingEnergies[fIndexOfShells[zElement] + i]*eV;
    if(energy >= bindingEnergy ) { break; }
  }

  // Normally one shell is available 
  if (i < nShells) { 
    bindingEnergy = fBindingEnergies[fIndexOfShells[zElement] + i]*eV;

    // update by deexcitation goes here

    energyOut = energy - bindingEnergy;
  }
  return energyOut;
}

#ifndef VECPHYS_NVCC
template<>
inline
VECPHYS_CUDA_HEADER_BOTH
typename kVc::Double_t
GUPhotoElectronSauterGavrila::
GetPhotoElectronEnergy<kVc>(typename kVc::Double_t energy,
                            typename kVc::Int_t  zElement) const
{
  kVc::Double_t energyOut;

  for(int i = 0; i < kVc::kSize ; ++i) {
    energyOut[i] = GetPhotoElectronEnergy<kScalar>(energy[i],zElement[i]);
  }

  return energyOut;
}
#endif

} // end namespace impl
} // end namespace vecphys

#endif
