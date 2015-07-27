#ifndef GUComptonKleinNishina_H
#define GUComptonKleinNishina_H 1

#include "backend/Backend.h"

#include "PhysicalConstants.h"
#include "GUTrack.h"
//#include "SystemOfUnits.h"

#include "GUAliasSampler.h"

// add the sincos function on MAC because sincos is not part of math.h
#ifdef __APPLE__ // possibly other conditions
inline void sincos(double x, double *s, double *c){
  __sincos(x,s,c);
}
#endif

namespace vecphys {

  //VECPHYS_DEVICE_DECLARE_CONV( GUComptonKleinNishina )
  //VECPHYS_DEVICE_FORWARD_DECLARE( class GUAliasSampler; )

inline namespace VECPHYS_IMPL_NAMESPACE {

class GUComptonKleinNishina
{
public:

  VECPHYS_CUDA_HEADER_HOST
  GUComptonKleinNishina(Random_t* states, int threadId = -1); 

  VECPHYS_CUDA_HEADER_BOTH
  GUComptonKleinNishina(Random_t* states, int threadId, 
                        GUAliasSampler* sampler); 

  VECPHYS_CUDA_HEADER_BOTH
  ~GUComptonKleinNishina();

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
  InteractKernel(typename Backend::double energyIn, 
                 typename Backend::Index_t   zElement,
                 typename Backend::double& energyOut,
                 typename Backend::double& sinTheta) const;


  // Alternative Implementation method(s) - for reference/comparisons
  // ----------------------------------------------------------------
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void SampleByCompositionRejection(typename Backend::double energyIn,
                                    typename Backend::double& energyOut,
                                    typename Backend::double& sinTheta);

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
  RotateAngle(typename Backend::double sinTheta,
              typename Backend::double xhat,
              typename Backend::double yhat,
              typename Backend::double zhat,
              typename Backend::double &xr,
              typename Backend::double &yr,
              typename Backend::double &zr) const;

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::double
  SampleSinTheta(typename Backend::double energyIn,
                 typename Backend::double energyOut) const; 

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::double 
  TotalCrossSection(typename Backend::double energy,
                    typename Backend::double Z) const;

private: 
  // Implementation methods 

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
VECPHYS_CUDA_HEADER_BOTH void 
GUComptonKleinNishina::InteractKernel(typename Backend::double  energyIn, 
                                      typename Backend::Index_t   zElement,
                                      typename Backend::double& energyOut,
                                          typename Backend::double& sinTheta)
                                      const
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::double double;

  Index_t   index;
  Index_t   icol;
  double  fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn,index,icol,fraction);

  double probNA;
  double aliasInd;

  //this did not used to work - Fixed SW
  fAliasSampler->GatherAlias<Backend>(index,zElement,probNA,aliasInd);
  
  double mininumE = energyIn/(1+2.0*energyIn*inv_electron_mass_c2);
  double deltaE = energyIn - mininumE;

  energyOut = mininumE + fAliasSampler->SampleX<Backend>(deltaE,probNA,
					        aliasInd,icol,fraction);
  sinTheta = SampleSinTheta<Backend>(energyIn,energyOut);
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
void
GUComptonKleinNishina::RotateAngle(typename Backend::double sinTheta,
                                   typename Backend::double xhat,
                                   typename Backend::double yhat,
                                   typename Backend::double zhat,
                                   typename Backend::double &xr,
                                   typename Backend::double &yr,
                                   typename Backend::double &zr) const
{
  typedef typename Backend::double double;
  typedef typename Backend::Bool_t   Bool_t;

  double phi = UniformRandom<Backend>(fRandomState,fThreadId);

  double pt = xhat*xhat + yhat*yhat;

  double cosphi, sinphi;
  sincos(phi, &sinphi, &cosphi);

  double uhat = sinTheta*cosphi; // cos(phi);
  double vhat = sinTheta*sinphi; // sin(phi);
  double what = Sqrt((1.-sinTheta)*(1.+sinTheta));

  Bool_t positive = ( pt > 0. );
  Bool_t negativeZ = ( zhat < 0. );

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

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::double 
GUComptonKleinNishina::
SampleSinTheta(typename Backend::double energyIn,
               typename Backend::double energyOut) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::double double;

  //angle of the scatterred photon

  double epsilon = energyOut/energyIn;

  Bool_t condition = epsilon > 1.0;

  MaskedAssign( condition, 1.0 , &epsilon );

  double E0_m    = inv_electron_mass_c2*energyIn;
  double onecost = (1.0 - epsilon)/(epsilon*E0_m);
  double sint2   = onecost*(2.-onecost);

  double sinTheta = 0.5;
  Bool_t condition2 = sint2 < 0.0;

  MaskedAssign(  condition2, 0.0, &sinTheta );   // Set sinTheta = 0
  MaskedAssign( !condition2, Sqrt(sint2), &sinTheta );   
  
  return sinTheta;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUComptonKleinNishina::
SampleByCompositionRejection(typename Backend::double  energyIn,
                             typename Backend::double& energyOut,
                             typename Backend::double& sinTheta)
{
  typedef typename Backend::double double;

  double epsilon, epsilonsq, onecost, sint2, greject ;
  
  double E0_m = energyIn*inv_electron_mass_c2 ;
  double eps0       = 1./(1. + 2.*E0_m);
  double epsilon0sq = eps0*eps0;
  double alpha1     = - log(eps0);
  double alpha2     = 0.5*(1.- epsilon0sq);
  
  do {
    if( alpha1/(alpha1+alpha2) > UniformRandom<Backend>(fRandomState,fThreadId) ) {
      epsilon   = exp(-alpha1*UniformRandom<Backend>(fRandomState,fThreadId));
      epsilonsq = epsilon*epsilon; 
    } 
    else {
      epsilonsq = epsilon0sq+(1.- epsilon0sq)*UniformRandom<Backend>(fRandomState,fThreadId);
      epsilon   = sqrt(epsilonsq);
    }
    
    onecost = (1.- epsilon)/(epsilon*E0_m);
    sint2   = onecost*(2.-onecost);
    greject = 1. - epsilon*sint2/(1.+ epsilonsq);
    
  } while (greject < UniformRandom<Backend>(fRandomState,fThreadId));
  
  energyOut = epsilon*energyIn;
  sinTheta = (sint2 < 0.0) ? 0.0 : sqrt(sint2);
}

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUComptonKleinNishina::Interact(GUTrack& inProjectile,
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
void GUComptonKleinNishina::Interact( GUTrack_v& inProjectile,    // In/Out
                                      const int *targetElements,  // Number equal to num of tracks
                                      GUTrack_v& outSecondary    // Empty vector for secondaries
                                      ) const
{
  typedef typename Backend::double double;
  typedef typename Backend::Index_t  Index_t;

  for(int j = 0; j < inProjectile.numTracks  ; ++j) {
     assert( (targetElements[j] > 0)  && (targetElements[j] <= fMaxZelement) );
  }
  
  int ibase= 0;
  int numChunks= (inProjectile.numTracks/double::Size);

  for(int i=0; i < numChunks ; ++i) {
    double energyIn(inProjectile.E[ibase]);
    double px(inProjectile.px[ibase]);
    double py(inProjectile.py[ibase]);
    double pz(inProjectile.pz[ibase]);
    double sinTheta;
    double energyOut;
    Index_t  zElement(targetElements[ibase]);

    InteractKernel<Backend>(energyIn, zElement, energyOut, sinTheta);

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
}    
#endif

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUComptonKleinNishina::ConvertXtoFinalState(double energyIn, 
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
void GUComptonKleinNishina::InteractG4(GUTrack& inProjectile,
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
typename Backend::double 
GUComptonKleinNishina::
TotalCrossSection(typename Backend::double energy,
                  typename Backend::double Z) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::double double;

  double sigma = 0.;

  //low energy limit
  Bool_t condition = energy > 10*keV;
  MaskedAssign( !condition, 0.0 , &sigma );

  if(Any(condition)) {
    double Z2 =  Z*Z;
    //put coff's to a constant header
    double p1 =  2.7965e-1 +  1.9756e-5*Z + -3.9178e-7*Z2;
    double p2 = -1.8300e-1 + -1.0205e-2*Z +  6.8241e-5*Z2;
    double p3 =  6.7527    + -7.3913e-2*Z +  6.0480e-5*Z2;
    double p4 = -1.9798e+1 +  2.7079e-2*Z +  3.0274e-4*Z2;
 
    double X = energy/electron_mass_c2;
    double X2 = X*Z;

    double tmpsigma = p1*log(1.+2.*X)/X
          + (p2 + p3*X + p4*X2)/(1. + 20.*X + 230.*X2 + 440.*X2*X);
    tmpsigma *= Z*barn;

    MaskedAssign( condition, tmpsigma , &sigma );
  }

  // 5% level improvements for low energy at 10-40keV for Hydrogen
  // and at 10-15keV for all Z    

  return sigma;
}

} // end namespace impl
} // end namespace vecphys

#endif
