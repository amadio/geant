#ifndef GUConversionBetheHeitler_H
#define GUConversionBetheHeitler_H 1

#include "backend/Backend.h"

#include "GUConstants.h"
#include "GUTrack.h"

// add the sincos function on MAC because sincos is not part of math.h
#ifdef __APPLE__ // possibly other conditions
inline void sincos(double x, double *s, double *c){
  __sincos(x,s,c);
}
#endif

namespace vecphys {

VECPHYS_DEVICE_FORWARD_DECLARE( class GUAliasSampler; )

inline namespace VECPHYS_IMPL_NAMESPACE {

class GUConversionBetheHeitler
{
public:

  VECPHYS_CUDA_HEADER_HOST
  GUConversionBetheHeitler(Random_t* states, int threadId = -1); 

  VECPHYS_CUDA_HEADER_BOTH
  GUConversionBetheHeitler(Random_t* states, int threadId, 
                           GUAliasSampler* sampler); 

  VECPHYS_CUDA_HEADER_BOTH
  ~GUConversionBetheHeitler();

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
                 typename Backend::Double_t& energyElectron,
                 typename Backend::Double_t& energyPositron,
                 typename Backend::Double_t& sinThetaElectron,
                 typename Backend::Double_t& sinThetaPositron) const;


  // Alternative Implementation method(s) - for reference/comparisons
  // ----------------------------------------------------------------
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void SampleByCompositionRejection(typename Backend::Int_t     elementZ,
                                    typename Backend::Double_t  energyIn,
                                    typename Backend::Double_t& energyElectron,
                                    typename Backend::Double_t& energyPositron,
                                    typename Backend::Double_t& sinThetaElectron,
                                    typename Backend::Double_t& sinThetaPositron);

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
  void
  SampleSinTheta(typename Backend::Double_t energyElectron,
                 typename Backend::Double_t energyPositron,
		 typename Backend::Double_t sinThetaElectron,
		 typename Backend::Double_t sinThetaPositron) const; 

private: 
  // Implementation methods 

  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton ) const;

  //this should be a method of GUElement
  VECPHYS_CUDA_HEADER_BOTH 
  double ComputeCoulombFactor(double Zeff) const;

  VECPHYS_CUDA_HEADER_BOTH 
  double ScreenFunction1(double screenVariable) const;

  VECPHYS_CUDA_HEADER_BOTH 
  double ScreenFunction2(double screenVariable) const;

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void ConvertXtoFinalState(double energyIn, 
                            double energyElectron, 
                            double energyPositron, 
                            double sinThetaElectron, 
                            double sinThetaPositron, 
                            GUTrack& primary, 
                            GUTrack& secondary) const;

private:
  GUAliasSampler* fAliasSampler; 

  // Helper data members for GPU random -- to be replaced by use of a GPU manager class
  Random_t* fRandomState;
  int       fThreadId;

  Precision fMinX;   // E Minimum - lowest energy for projectile
  Precision fMaxX;

  Precision fMaxZelement; // 
  
  //Sampling Tables
  int fNrow;
  int fNcol;
};

//Implementation

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH 
void 
GUConversionBetheHeitler::InteractKernel(typename Backend::Double_t  energyIn, 
                                         typename Backend::Index_t   zElement,
                                         typename Backend::Double_t& energyElectron,
                                         typename Backend::Double_t& energyPositron,
                                         typename Backend::Double_t& sinThetaElectron,
                                         typename Backend::Double_t& sinThetaPositron)
                                         const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::Double_t Double_t;

  Index_t   index;
  Index_t   icol;
  Double_t  fraction;

  fAliasSampler->SampleBin<Backend>(energyIn,index,icol,fraction);

  Double_t probNA;
  Double_t aliasInd;

  //this did not used to work - Fixed SW
  fAliasSampler->GatherAlias<Backend>(index,zElement,probNA,aliasInd);
  
  Double_t mininumE = electron_mass_c2;
  Double_t deltaE = energyIn - mininumE;

  //electron energy
  Double_t energyOut = mininumE + fAliasSampler->SampleX<Backend>(deltaE,probNA,
					        aliasInd,icol,fraction);

  Bool_t condition = 0.5 > UniformRandom<Backend>(fRandomState,fThreadId);
  MaskedAssign( condition, energyOut, &energyElectron);     
  MaskedAssign( condition, energyIn - energyOut, &energyPositron);     

  MaskedAssign(!condition, energyOut, &energyPositron);     
  MaskedAssign(!condition, energyIn - energyOut, &energyElectron);     
  
  SampleSinTheta<Backend>(energyIn,energyElectron,energyPositron,
                          sinThetaElectron, sinThetaPositron);
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
void
GUConversionBetheHeitler::RotateAngle(typename Backend::Double_t sinTheta,
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
void
GUConversionBetheHeitler::
SampleSinTheta(typename Backend::Double_t energyElectron,
               typename Backend::Double_t energyPositron,
	       typename Backend::Double_t sinThetaElectron,
	       typename Backend::Double_t sinThetaPositron) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  //angles of the pair production (gamma -> e+e-)

  Double_t u;
  const double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  Bool_t condition = 9./(9. + d) > UniformRandom<Backend>(fRandomState,fThreadId);
  MaskedAssign( condition, -log( UniformRandom<Backend>(fRandomState,fThreadId)*
                       UniformRandom<Backend>(fRandomState,fThreadId))/a1, &u );
  MaskedAssign(!condition, -log( UniformRandom<Backend>(fRandomState,fThreadId)*
                       UniformRandom<Backend>(fRandomState,fThreadId))/a2, &u );

  Double_t TetEl = u*electron_mass_c2/energyElectron;
  Double_t TetPo = u*electron_mass_c2/energyPositron;

  //sinTheta - just return theta instead!
  sinThetaElectron =  sin(TetEl);
  sinThetaPositron = -sin(TetPo);
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUConversionBetheHeitler::
SampleByCompositionRejection(typename Backend::Int_t     elementZ,
                             typename Backend::Double_t  energyIn,
                             typename Backend::Double_t& energyElectron,
                             typename Backend::Double_t& energyPositron,
                             typename Backend::Double_t& sinThetaElectron,
                             typename Backend::Double_t& sinThetaPositron)
{
  typedef typename Backend::Double_t Double_t;

  Double_t epsil ;
  Double_t epsil0 = electron_mass_c2/energyIn ;
  if(epsil0 > 1.0) { return; }

  //  static 
  const double Egsmall=2.; //2.*MeV;

  // select randomly one element constituing the material - input

  if (energyIn < Egsmall) {
    epsil = epsil0+(0.5-epsil0)*UniformRandom<Backend>(fRandomState,fThreadId);
  } else {
    // Extract Coulomb factor for this element
    int logZ3 = log(1.0*int(elementZ + 0.5))/3.0; 
    double FZ = 8.*logZ3 ; //8.*(anElement->GetIonisation()->GetlogZ3());

    if (energyIn > 50. /*MeV*/ ) { FZ += 8.*ComputeCoulombFactor(elementZ);}
    // limits of the screening variable
    int Z3 = pow(1.0*int(elementZ + 0.5),1/3.0);
    double screenfac = 136.*epsil0/Z3; //(anElement->GetIonisation()->GetZ3());
    double screenmax = exp ((42.24 - FZ)/8.368) - 0.952 ;
    double screenmin = fmin(4.*screenfac,screenmax);

    // limits of the energy sampling
    double epsil1 = 0.5 - 0.5*sqrt(1. - screenmin/screenmax) ;
    double epsilmin = fmax(epsil0,epsil1) , epsilrange = 0.5 - epsilmin;

    // sample the energy rate of the created electron (or positron)

    double  screenvar, greject ;

    double F10 = ScreenFunction1(screenmin) - FZ;
    double F20 = ScreenFunction2(screenmin) - FZ;

    double NormF1 = fmax(F10*epsilrange*epsilrange,0.); 
    double NormF2 = fmax(1.5*F20,0.);
    do {
      if( NormF1/(NormF1+NormF2) > UniformRandom<Backend>(fRandomState,fThreadId)) {
        epsil = 0.5 - 
          epsilrange*pow(UniformRandom<Backend>(fRandomState,fThreadId), 0.333333);
        screenvar = screenfac/(epsil*(1-epsil));
        greject = (ScreenFunction1(screenvar) - FZ)/F10;
              
      } else { 
        epsil = epsilmin + epsilrange*UniformRandom<Backend>(fRandomState,fThreadId);
        screenvar = screenfac/(epsil*(1-epsil));
        greject = (ScreenFunction2(screenvar) - FZ)/F20;
      }
    } while( greject < UniformRandom<Backend>(fRandomState,fThreadId) );
  }   //  end of epsil sampling

  if (UniformRandom<Backend>(fRandomState,fThreadId) > 0.5) {
    energyElectron = (1.-epsil)*energyIn;
    energyPositron = epsil*energyIn;
  } else {
    energyPositron = (1.-epsil)*energyIn;
    energyElectron = epsil*energyIn;
  }

  //sample angle
  double u;
  const double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  //(9./(9.+d) = 0.25
  if (9./(9.+d) > UniformRandom<Backend>(fRandomState,fThreadId)) {
    u= - log( UniformRandom<Backend>(fRandomState,fThreadId)*UniformRandom<Backend>(fRandomState,fThreadId))/a1;
  }
  else {                            
    u= - log( UniformRandom<Backend>(fRandomState,fThreadId)*UniformRandom<Backend>(fRandomState,fThreadId))/a2;
  }

  //move to the constant file
  const double twopi = 2.*3.14159265358979323846;

  double TetEl = u*electron_mass_c2/energyElectron;
  double TetPo = u*electron_mass_c2/energyPositron;
  double Phi  = twopi * UniformRandom<Backend>(fRandomState,fThreadId);

  //return sinTheta
  sinThetaElectron = sin(TetEl);
  sinThetaPositron =-sin(TetPo);

}

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUConversionBetheHeitler::Interact(GUTrack& inProjectile,
                                        int      targetElement,
                                        GUTrack& outSecondary ) const
{
  double energyIn= inProjectile.E;
  double energyElectron, sinThetaElectron;
  double energyPositron, sinThetaPositron;
#ifdef CHECK
  if( (energyIn <= fMinX) || (energyIn > fMaxX) )
  {
    printf(" Illegal input Energy = %f min = %f max = %f\n",
	   energyIn,fMinX,fMaxX);
  }
#endif 
  //  assert( (energyIn >= fMinX)  && (energyIn <= fMaxX) );
  InteractKernel<Backend>(energyIn, targetElement,energyElectron, energyPositron, sinThetaElectron, sinThetaPositron);
  
  //update final states of the primary and store the secondary
  ConvertXtoFinalState<Backend>(energyIn,energyElectron, energyPositron,sinThetaElectron,sinThetaPositron,
                                inProjectile,outSecondary);
}
  
#ifndef VECPHYS_NVCC
template <typename Backend>
void GUConversionBetheHeitler::Interact(GUTrack_v& inProjectile,  
                                        const int *targetElements,
                                        GUTrack_v& outSecondary ) const
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
    Double_t sinThetaElectron;
    Double_t energyElectron;
    Double_t sinThetaPositron;
    Double_t energyPositron;
    Index_t  zElement(targetElements[ibase]);

    InteractKernel<Backend>(energyIn, zElement, energyElectron, energyPositron,
                            sinThetaElectron,sinThetaPositron);

    //need to rotate the angle with respect to the line of flight
    Double_t invp = 1./energyIn;
    Double_t xhat = px*invp;
    Double_t yhat = py*invp;
    Double_t zhat = pz*invp;

    Double_t uhat = 0.;
    Double_t vhat = 0.;
    Double_t what = 0.;


    // Kill the primary photon: use it as a placeholder for the secondary 
    // positron temporarily - change (charge, pid) accordingly 

    RotateAngle<Backend>(sinThetaPositron,xhat,yhat,zhat,uhat,vhat,what);

    energyPositron.store(&inProjectile.E[ibase]);
    Double_t pxFinal, pyFinal, pzFinal;
     
    pxFinal= energyPositron*uhat;
    pyFinal= energyPositron*vhat;
    pzFinal= energyPositron*what;
    pxFinal.store(&inProjectile.px[ibase]);
    pyFinal.store(&inProjectile.py[ibase]);
    pzFinal.store(&inProjectile.pz[ibase]);

    // create Secondary
    RotateAngle<Backend>(sinThetaElectron,xhat,yhat,zhat,uhat,vhat,what);

    Double_t secE = energyElectron; 
    Double_t pxSec= secE*uhat;
    Double_t pySec= secE*vhat;
    Double_t pzSec= secE*what;

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
void GUConversionBetheHeitler::ConvertXtoFinalState(double energyIn, 
                                                    double energyElectron, 
                                                    double energyPositron, 
                                                    double sinThetaElectron, 
                                                    double sinThetaPositron, 
                                                    GUTrack& inProjectile,
                                                    GUTrack& outSecondary) const
{
  //need to rotate the angle with respect to the line of flight
  double invp = 1./energyIn;
  double xhat = inProjectile.px*invp;
  double yhat = inProjectile.py*invp;
  double zhat = inProjectile.pz*invp;

  double uhat = 0.;
  double vhat = 0.;
  double what = 0.;


  // Kill the primary photon: use it as a placeholder for the secondary 
  // positron temporarily - change (charge, pid) accordingly 
  RotateAngle<Backend>(sinThetaPositron,xhat,yhat,zhat,uhat,vhat,what);
  inProjectile.E  = energyPositron;
  inProjectile.px = energyPositron*uhat;
  inProjectile.py = energyPositron*vhat;
  inProjectile.pz = energyPositron*what;

  //create secondary
  RotateAngle<Backend>(sinThetaElectron,xhat,yhat,zhat,uhat,vhat,what);
  outSecondary.E  = energyElectron; 
  outSecondary.px = outSecondary.E*uhat;
  outSecondary.py = outSecondary.E*vhat;
  outSecondary.pz = outSecondary.E*what;
  //fill other information
}

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUConversionBetheHeitler::InteractG4(GUTrack& inProjectile,
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
