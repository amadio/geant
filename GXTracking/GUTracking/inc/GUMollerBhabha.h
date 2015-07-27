#ifndef GUMollerBhabha_H
#define GUMollerBhabha_H 1

#include "backend/Backend.h"

#include "GUConstants.h"
#include "GUTrack.h"
#include "PhysicalConstants.h"

namespace vecphys {

  //VECPHYS_DEVICE_DECLARE_CONV( GUMollerBhabha )
VECPHYS_DEVICE_FORWARD_DECLARE( class GUAliasSampler; )

inline namespace VECPHYS_IMPL_NAMESPACE {

class GUMollerBhabha
{
public:

  VECPHYS_CUDA_HEADER_HOST
  GUMollerBhabha(Random_t* states, int threadId = -1); 

  VECPHYS_CUDA_HEADER_BOTH
  GUMollerBhabha(Random_t* states, int threadId, 
                        GUAliasSampler* sampler); 

  VECPHYS_CUDA_HEADER_BOTH
  ~GUMollerBhabha();

  // Generate secondaries 
  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void Interact(GUTrack& projectile,    //In/Out: Updated to new state - choice
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


private: 
  // Implementation methods 

  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSection( int Zelement, double Ein, double Eout) const;

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

  Precision fMinX;
  Precision fMaxX;
  //  Precision fDeltaX;

  //  Precision fMinY;
  //  Precision fMaxY;
  //  Precision fDeltaY;

  Precision fMaxZelement; 

  //Sampling Tables
  int fNrow;
  int fNcol;
};

//Implementation

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void 
GUMollerBhabha::InteractKernel(typename Backend::double  energyIn, 
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
  
  double mininumE = 0.1*keV;
  double deltaE = energyIn/2.0 - mininumE;

  energyOut = mininumE + fAliasSampler->SampleX<Backend>(deltaE,probNA,
                                                aliasInd,icol,fraction);
  sinTheta = SampleSinTheta<Backend>(energyIn,energyOut);
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
void
GUMollerBhabha::RotateAngle(typename Backend::double sinTheta,
                            typename Backend::double xhat,
                            typename Backend::double yhat,
                            typename Backend::double zhat,
                            typename Backend::double &xr,
                            typename Backend::double &yr,
                            typename Backend::double &zr) const
{
  typedef typename Backend::double double;
  typedef typename Backend::bool   bool;

  double phi = UniformRandom<Backend>(fRandomState,fThreadId);

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

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::double 
GUMollerBhabha::
SampleSinTheta(typename Backend::double energyIn,
               typename Backend::double energyOut) const
{
  typedef typename Backend::bool   bool;
  typedef typename Backend::double double;

  //angle of the scatterred electron

  double energy = energyIn + electron_mass_c2;
  double totalMomentum = sqrt(energyIn*(energy + electron_mass_c2));

  double deltaMomentum = sqrt(energyOut * (energyOut + 2.0*electron_mass_c2));
  double cost =  energyOut * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum);
  double sint2 = (1.0 - cost)*(1. + cost);

  double sinTheta = 0.5;
  bool condition2 = sint2 < 0.0;

  MaskedAssign(  condition2, 0.0, &sinTheta );   // Set sinTheta = 0
  MaskedAssign( !condition2, Sqrt(sint2), &sinTheta );   
  
  return sinTheta;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUMollerBhabha::
SampleByCompositionRejection(typename Backend::double  energyIn,
			     typename Backend::double& energyOut,
			     typename Backend::double& sinTheta)
{
  typedef typename Backend::double double;

  //cut energy
  double tmin = 0.1*keV; // minimum delta-ray energy which should be setable  
  //  double tmax = fLimit*energyIn; //fLimit=0.5 for e- and 1 for e+
  double tmax = 0.5*energyIn; //fLimit=0.5 for e- and 1 for e+

  if(fMaxX < tmax) { tmax = fMaxX; }
  if(tmin >= tmax) { return; }

  double energy = energyIn + electron_mass_c2;
  double totalMomentum = sqrt(energyIn*(energy + electron_mass_c2));

  double xmin   = tmin/energyIn;
  double xmax   = tmax/energyIn;
  double gam    = energy/electron_mass_c2;
  double gamma2 = gam*gam;
  double beta2  = 1.0 - 1.0/gamma2;
  double x, z, q, grej;

  //Moller (e-e-) scattering 
  //Bhabha (e+e-) scattering if positron (not implemented for now)

  double gg = (2.0*gam - 1.0)/gamma2;
  double y = 1.0 - xmax;
  grej = 1.0 - gg*xmax + xmax*xmax*(1.0 - gg + (1.0 - gg*y)/(y*y));
  
  do {
    q =  UniformRandom<Backend>(fRandomState,fThreadId);
    x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
    y = 1.0 - x;
    z = 1.0 - gg*x + x*x*(1.0 - gg + (1.0 - gg*y)/(y*y));
  } while (grej * UniformRandom<Backend>(fRandomState,fThreadId) > z);
  
  energyOut = x*energyIn;

  double deltaMomentum = sqrt(energyOut * (energyOut + 2.*electron_mass_c2));
  double cost = energyOut * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum);
  double sint2 = (1.0 - cost)*(1. + cost);

  sinTheta = (sint2 < 0.0) ? 0.0 : sqrt(sint2);
}

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUMollerBhabha::Interact(GUTrack& inProjectile,
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
void GUMollerBhabha::Interact( GUTrack_v& inProjectile,    // In/Out
                               const int *targetElements,  // Number equal to num of tracks
                               GUTrack_v& outSecondary    // Empty vector for secondaries
                                      ) const
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::double double;

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
void GUMollerBhabha::ConvertXtoFinalState(double energyIn, 
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
void GUMollerBhabha::InteractG4(GUTrack& inProjectile,
                                int      targetElement,
                                GUTrack& outSecondary )
{
  Precision energyIn;

  energyIn = inProjectile.E;

  Precision energyOut = -999.;
  Precision sinTheta = -999.;
  SampleByCompositionRejection<Backend>(energyIn,energyOut,sinTheta);

  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
  
}

} // end namespace impl
} // end namespace vecphys

#endif
