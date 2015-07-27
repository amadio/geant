#ifndef GUSeltzerBerger_H
#define GUSeltzerBerger_H 1

#include "backend/Backend.h"

#include "GUConstants.h"
#include "GUTrack.h"
#include "PhysicalConstants.h"
#include "Physics2DVector.h"
#include <fstream>

namespace vecphys {

VECPHYS_DEVICE_FORWARD_DECLARE( class GUAliasSampler; )

inline namespace VECPHYS_IMPL_NAMESPACE {

class GUSeltzerBerger
{
public:

  VECPHYS_CUDA_HEADER_HOST
  GUSeltzerBerger(Random_t* states, int threadId = -1); 

  VECPHYS_CUDA_HEADER_BOTH
  GUSeltzerBerger(Random_t* states, int threadId, 
		  GUAliasSampler* sampler,
		  Physics2DVector* sbData); 

  VECPHYS_CUDA_HEADER_BOTH
  ~GUSeltzerBerger();

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
  void SampleByCompositionRejection(typename Backend::int     elementZ,
                                    typename Backend::double  energyIn,
                                    typename Backend::double& energyOut,
                                    typename Backend::double& sinTheta);

  //  Initialisation methods
  // -------------------------------------------

  // Initializes this class and its sampler 
  VECPHYS_CUDA_HEADER_HOST
  void BuildOneTable(int Z,
                     const double xmin,
                     const double xmax,
                     const int nrow,
                     const int ncol);

  VECPHYS_CUDA_HEADER_HOST
  void BuildLogPdfTable(int Z,
                        const double xmin,
                        const double xmax,
                        const int nrow,
                        const int ncol,
                        double *p);

  VECPHYS_CUDA_HEADER_HOST bool
  RetrieveSeltzerBergerData(std::ifstream& in, Physics2DVector *vec2D);

public:
  // Auxiliary methods
  VECPHYS_CUDA_HEADER_BOTH
  GUAliasSampler* GetSampler() {return fAliasSampler;}

  VECPHYS_CUDA_HEADER_BOTH
  void SetSampler(GUAliasSampler* sampler) { fAliasSampler = sampler ;}

  VECPHYS_CUDA_HEADER_BOTH
  Physics2DVector* GetSBData() {return fDataSB;}

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
  SampleSinTheta(typename Backend::double energyIn) const;

private: 
  // Implementation methods 

  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSection(int Zelement, double Ein, double Eout) const;

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void ConvertXtoFinalState(double energyIn, 
			    double energyOut, 
			    double sinTheta, 
			    GUTrack& primary, 
			    GUTrack& secondary) const;

private:
  GUAliasSampler* fAliasSampler; 
  Physics2DVector* fDataSB;

  // Helper data members for GPU random
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
GUSeltzerBerger::InteractKernel(typename Backend::double  energyIn, 
                                typename Backend::Index_t   zElement,
                                typename Backend::double& energyOut,
                                typename Backend::double& sinTheta) const
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

  //Seltzer-Berger specific
  
  // To-do: apply densityFactor (dummy for now) and calculate deltaY correctly
  // densityFactor = (Migdal constant)x(electron density of the material); 
  double densityFactor = 1.0;
  
  double emin = Min(fMinX, energyIn);
  double emax = Min(fMaxX, energyIn);

  double totalEnergy = energyIn + electron_mass_c2;
  double densityCorr = densityFactor*totalEnergy*totalEnergy;
  double minY = Log(emin*emin + densityCorr);
  double maxY = Log(emax*emax + densityCorr);
  double deltaY = maxY - minY;

  double yhat = fAliasSampler->SampleX<Backend>(deltaY,probNA,
                                                  aliasInd,icol,fraction);

  energyOut =  Sqrt(Max(Exp(minY + yhat)- densityCorr,0.0));
  sinTheta = SampleSinTheta<Backend>(energyOut);
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
void
GUSeltzerBerger::RotateAngle(typename Backend::double sinTheta,
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
GUSeltzerBerger::
SampleSinTheta(typename Backend::double energyIn) const
{
  typedef typename Backend::bool   bool;
  typedef typename Backend::double double;

  //angle of the radiated photon 
  //based on G4DipBustGenerator::SampleDirection
  
  double c = 4. - 8.*UniformRandom<Backend>(fRandomState,fThreadId);
  double a;
  double signc; 
  bool condition = c > 0.;
  MaskedAssign(  condition,  1. , &signc );
  MaskedAssign( !condition, -1. , &signc );
  MaskedAssign(  condition,  c , &a );
  MaskedAssign( !condition, -c , &a );

  double delta  = Sqrt(a*a+4.);
  delta += a;
  delta *= 0.5; 

  //To-do:  Vc does not support pow 
  //  double cofA = -signc*Pow(delta, 1./3.);
  double cofA = -signc*Sqrt(delta); //temporary replace Sqrt by pow

  double cosTheta = cofA - 1./cofA;

  double tau  = energyIn/electron_mass_c2;
  double beta = Sqrt(tau*(tau + 2.))/(tau + 1.);

  cosTheta = (cosTheta + beta)/(1 + cosTheta*beta);

  double sinTheta = Sqrt((1 - cosTheta)*(1 + cosTheta));

  return sinTheta;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUSeltzerBerger::
SampleByCompositionRejection(typename Backend::int     elementZ,
			     typename Backend::double  energyIn,
			     typename Backend::double& energyOut,
			     typename Backend::double& sinTheta)
{
  typedef typename Backend::double double;

  //check validity of elementZ

  //based on G4SeltzerBergerModel::SampleSecondaries

  double kineticEnergy = energyIn ;
  double emin = Min(fMinX, energyIn );
  double emax = Min(fMaxX, energyIn);

  if(emin >= emax) { return; }

  double densityFactor =1.0;
  
  double totalEnergy = energyIn + electron_mass_c2;
  double densityCorr = densityFactor*totalEnergy*totalEnergy;
  double totMomentum = Sqrt(kineticEnergy*(totalEnergy + electron_mass_c2));

  double xmin = Log(emin*emin + densityCorr);
  double xmax = Log(emax*emax + densityCorr);
  double    y = Log(energyIn/MeV);

  double gammaEnergy;
  double v; 

  // majoranta
  double x0 = emin/energyIn;
  double vmax = fDataSB[elementZ].Value(x0, y)*1.02;

  const double epeaklimit= 300*MeV; 
  const double elowlimit = 10*keV; 

  // majoranta corrected for e-
  bool isElectron = true;
  if(isElectron && x0 < 0.97 && 
     ((energyIn > epeaklimit) || (energyIn < elowlimit))) {
    double ylim = Min(fDataSB[elementZ].Value(0.97, 4*log(10.)),
                         1.1*fDataSB[elementZ].Value(0.97, y)); 
    if(ylim > vmax) { vmax = ylim; }
  }
  if(x0 < 0.05) { vmax *= 1.2; }

  do {
    double auxrand = UniformRandom<Backend>(fRandomState,fThreadId);
    double x = exp(xmin + auxrand*(xmax - xmin))-densityCorr;
    if(x < 0.0) { x = 0.0; }
    energyOut = sqrt(x);
    double x1 = energyOut/energyIn;
    v = fDataSB[elementZ].Value(x1, y);

    // correction for positrons - uncomment this if we add positron       
    /*    
    if(!isElectron) {
      double e1 = kineticEnergy - emin;
      double invbeta1 = (e1 + fMass)/sqrt(e1*(e1 + 2*fMass));
      double e2 = kineticEnergy - gammaEnergy;
      double invbeta2 = (e2 + fMass)/sqrt(e2*(e2 + 2*fMass));
      double xxx = twopi*fine_structure_const*currentZ*(invbeta1 - invbeta2); 
      
      if(xxx < -12. ) { v = 0.0; }
      else { v *= exp(xxx); }
    }
    */

  } while (v < vmax*UniformRandom<Backend>(fRandomState,fThreadId));

  //anagle
  sinTheta = SampleSinTheta<Backend>(energyOut);
}

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUSeltzerBerger::Interact(GUTrack& inProjectile,
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
void GUSeltzerBerger::Interact(GUTrack_v& inProjectile,   
                               const int *targetElements, 
                               GUTrack_v& outSecondary) const
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
void GUSeltzerBerger::ConvertXtoFinalState(double energyIn, 
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
void GUSeltzerBerger::InteractG4(GUTrack& inProjectile,
                                 int      targetElement,
                                 GUTrack& outSecondary )
{
  Precision energyIn;

  energyIn = inProjectile.E;

  Precision energyOut;
  Precision sinTheta;
  SampleByCompositionRejection<Backend>(targetElement, 
                                        energyIn, energyOut, sinTheta);

  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
  
}

} // end namespace impl
} // end namespace vecphys

#endif
