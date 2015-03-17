#ifndef GUMollerBhabha_H
#define GUMollerBhabha_H 1

#include "backend/Backend.h"

#include "GUConstants.h"
#include "GUTrack.h"

namespace vecphys {

  //VECPHYS_DEVICE_DECLARE_CONV( GUMollerBhabha )
VECPHYS_DEVICE_FORWARD_DECLARE( class GUAliasSampler; )

inline namespace VECPHYS_IMPL_NAMESPACE {

class GUMollerBhabha
{
public:

  VECPHYS_CUDA_HEADER_BOTH
  GUMollerBhabha(Random_t* states, int threadId = -1); 

  VECPHYS_CUDA_HEADER_BOTH
  GUMollerBhabha(Random_t* states, int threadId, 
                        GUAliasSampler* sampler); 

  VECPHYS_CUDA_HEADER_BOTH
  ~GUMollerBhabha();

  VECPHYS_CUDA_HEADER_BOTH
  GUAliasSampler* GetSampler() {return fAliasSampler;}

  VECPHYS_CUDA_HEADER_BOTH
  void SetSampler(GUAliasSampler* sampler) { fAliasSampler = sampler ;}
  
  VECPHYS_CUDA_HEADER_BOTH
  void GetSampleParameters(double x, int &irow, int &icol,double &t);

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

  // Initializes this class and its sampler 
  VECPHYS_CUDA_HEADER_BOTH
  void BuildTable( int Z,
                              const double xmin,
                              const double xmax,
                              const int nrow,
			      const int ncol);
  // QUESTION: This might depend on physics? So maybe we should place inside model? 

  VECPHYS_CUDA_HEADER_BOTH
  void BuildPdfTable(int Z,
                     const double xmin,
                     const double xmax,
                     const int nrow,
                     const int ncol,
                     double *p);

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH void 
  InteractKernel(typename Backend::Double_t energyIn, 
                 typename Backend::Double_t& energyOut,
		 typename Backend::Double_t& sinTheta) const;


  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  SampleSinTheta(typename Backend::Double_t energyIn,
                 typename Backend::Double_t energyOut) const; 


  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void SampleByCompositionRejection(typename Backend::Double_t energyIn,
				    typename Backend::Double_t& energyOut,
				    typename Backend::Double_t& sinTheta);

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

  Precision fMinX;
  Precision fMaxX;
  Precision fDeltaX;

  Precision fMinY;
  Precision fMaxY;
  Precision fDeltaY;

  //Sampling Tables
  int fNrow;
  int fNcol;
};

//Implementation

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
void
GUMollerBhabha::RotateAngle(typename Backend::Double_t sinTheta,
                                   typename Backend::Double_t xhat,
                                   typename Backend::Double_t yhat,
                                   typename Backend::Double_t zhat,
                                   typename Backend::Double_t &xr,
                                   typename Backend::Double_t &yr,
                                   typename Backend::Double_t &zr) const
{
  typedef typename Backend::Double_t Double_t;
  typedef typename Backend::Bool_t   Bool_t;

  Double_t phi = UniformRandom(fRandomState,fThreadId);

  Double_t pt = xhat*xhat + yhat*yhat;

  Double_t uhat = sinTheta*cos(phi);
  Double_t vhat = sinTheta*sin(phi);
  Double_t what = Sqrt((1.-sinTheta)*(1.+sinTheta));

  Bool_t positive = ( pt > 0. );
  Bool_t negative = ( zhat < 0. );

  //mask operation???
  if(positive) {
    Double_t phat = Sqrt(pt);
    xr = (xhat*zhat*uhat - yhat*vhat)/phat + xhat*what;
    yr = (yhat*zhat*uhat - xhat*vhat)/phat + yhat*what;
    zr = -phat*uhat + zhat*what;
  }
  else if(negative) {
    xr = -xhat;
    yr =  yhat;
    zr = -zhat;
  }  

}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t 
GUMollerBhabha::
SampleSinTheta(typename Backend::Double_t energyIn,
               typename Backend::Double_t energyOut) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  //angle of the scatterred electron

  Double_t energy = energyIn + electron_mass_c2;
  Double_t totalMomentum = sqrt(energyIn*(energy + electron_mass_c2));

  Double_t deltaMomentum = sqrt(energyOut * (energyOut + 2.0*electron_mass_c2));
  Double_t cost =  energyOut * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum);
  Double_t sint2 = (1.0 - cost)*(1. + cost);

  Double_t sinTheta = 0.5;
  Bool_t condition2 = sint2 < 0.0;

  MaskedAssign(  condition2, 0.0, &sinTheta );   // Set sinTheta = 0
  MaskedAssign( !condition2, Sqrt(sint2), &sinTheta );   
  
  return sinTheta;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUMollerBhabha::
SampleByCompositionRejection(typename Backend::Double_t  energyIn,
			     typename Backend::Double_t& energyOut,
			     typename Backend::Double_t& sinTheta)
{
  typedef typename Backend::Double_t Double_t;

  //cut energy
  Double_t tmin = 1.0;  
  //  Double_t tmax = fLimit*energyIn; //fLimit=0.5 for e- and 1 for e+
  Double_t tmax = 0.5*energyIn; //fLimit=0.5 for e- and 1 for e+

  if(fMaxX < tmax) { tmax = fMaxX; }
  if(tmin >= tmax) { return; }

  Double_t energy = energyIn + electron_mass_c2;
  Double_t totalMomentum = sqrt(energyIn*(energy + electron_mass_c2));

  Double_t xmin   = tmin/energyIn;
  Double_t xmax   = tmax/energyIn;
  Double_t gam    = energy/electron_mass_c2;
  Double_t gamma2 = gam*gam;
  Double_t beta2  = 1.0 - 1.0/gamma2;
  Double_t x, z, q, grej;

  //Moller (e-e-) scattering 
  //Bhabha (e+e-) scattering if positron (not implemented for now)

  Double_t gg = (2.0*gam - 1.0)/gamma2;
  Double_t y = 1.0 - xmax;
  grej = 1.0 - gg*xmax + xmax*xmax*(1.0 - gg + (1.0 - gg*y)/(y*y));
  
  do {
    q = UniformRandom(fRandomState, fThreadId);
    x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
    y = 1.0 - x;
    z = 1.0 - gg*x + x*x*(1.0 - gg + (1.0 - gg*y)/(y*y));
  } while (grej * UniformRandom(fRandomState,fThreadId) > z);
  
  energyOut = x*energyIn;

  Double_t deltaMomentum = sqrt(energyOut * (energyOut + 2.*electron_mass_c2));
  Double_t cost = energyOut * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum);
  Double_t sint2 = (1.0 - cost)*(1. + cost);

  sinTheta = (sint2 < 0.0) ? 0.0 : sqrt(sint2);
}

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUMollerBhabha::Interact(GUTrack& inProjectile,
                                     int      targetElement,
                                     GUTrack& outSecondary ) const
{
  double energyIn;
  double deltaE; //temporary - this should be dy in BuildPdfTable

  energyIn = inProjectile.E;
  deltaE =  energyIn - energyIn/(1+2.0*energyIn*inv_electron_mass_c2);

  int index;
  int icol;
  double fraction;

  fAliasSampler->SampleBin<Backend>(energyIn,index,icol,fraction);

  double probNA;   // Non-alias probability
  int aliasInd; 

  //  This is really an integer -- could be In  
  fAliasSampler->GetAlias(index,probNA,aliasInd);

  double energyOut = fAliasSampler->SampleX<Backend>(deltaE,probNA,aliasInd,
					       icol,fraction);

  //calcuate the scattered angle
  double sinTheta = SampleSinTheta<Backend>(energyIn,energyOut);

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
  typedef typename Backend::Double_t Double_t;

  Double_t energyIn;
  Double_t deltaE; //temporary - this should be dy in BuildPdfTable
  Index_t  index;
  Index_t  icol;
  Double_t fraction;

  Double_t px, py, pz;

  for(int i=0; i < inProjectile.numTracks/Double_t::Size ; ++i) {
    
    //gather
    // loads energies into a VC-"register" type called energyIn
    for(int j = 0; j < Double_t::Size ; ++j) {
      energyIn[j] = inProjectile.E[ i*Double_t::Size + j];
      deltaE[j] = energyIn[j] - energyIn[j]/(1+2.0*energyIn[j]*inv_electron_mass_c2);

      px[j] =  inProjectile.px[ i*Double_t::Size + j];
      py[j] =  inProjectile.py[ i*Double_t::Size + j];
      pz[j] =  inProjectile.pz[ i*Double_t::Size + j];
    }

    fAliasSampler->SampleBin<Backend>(energyIn,index,icol,fraction);

    Double_t probNA;   // Non-alias probability
    Double_t aliasInd; // This is really an integer -- could be Index_t !?

    //gather for alias table lookups
    fAliasSampler->GatherAlias<Backend>(index,probNA,aliasInd);

    Double_t energyOut = fAliasSampler->SampleX<Backend>(deltaE,probNA,aliasInd,
						         icol,fraction);

    //calcuate the scattered angle
    Double_t sinTheta = SampleSinTheta<Backend>(energyIn,energyOut);

    //need to rotate the angle with respect to the line of flight
    Double_t invp = 1./energyIn;
    Double_t xhat = px*invp;
    Double_t yhat = py*invp;
    Double_t zhat = pz*invp;

    Double_t uhat = 0.;
    Double_t vhat = 0.;
    Double_t what = 0.;

    RotateAngle<Backend>(sinTheta,xhat,yhat,zhat,uhat,vhat,what);

    //scatter 
    for(int j = 0; j < Double_t::Size ; ++j) {
      int it = i*Double_t::Size + j;

      //update primary
      inProjectile.E[it]  = energyOut[j];
      inProjectile.px[it] = energyOut[j]*uhat[j];
      inProjectile.py[it] = energyOut[j]*vhat[j];
      inProjectile.pz[it] = energyOut[j]*what[j];

      //create secondary
      double secE = energyIn[j] - energyOut[j]; 
      outSecondary.E[it]  = secE;
      outSecondary.px[it] = secE*(xhat[j]-uhat[j]);
      outSecondary.py[it] = secE*(yhat[j]-vhat[j]);
      outSecondary.pz[it] = secE*(zhat[j]-what[j]);
      //fill other information
    }
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

  Precision energyOut;
  Precision sinTheta;
  SampleByCompositionRejection<Backend>(energyIn,energyOut,sinTheta);

  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
  
}

} // end namespace impl
} // end namespace vecphys

#endif
