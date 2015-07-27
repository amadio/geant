#ifndef GUConversionBetheHeitler_H
#define GUConversionBetheHeitler_H 1

#include "backend/Backend.h"

#include "GUConstants.h"
#include "GUTrack.h"
#include "PhysicalConstants.h"

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
  InteractKernel(typename Backend::double energyIn, 
                 typename Backend::Index_t   zElement,
                 typename Backend::double& energyElectron,
                 typename Backend::double& energyPositron,
                 typename Backend::double& sinThetaElectron,
                 typename Backend::double& sinThetaPositron) const;


  // Alternative Implementation method(s) - for reference/comparisons
  // ----------------------------------------------------------------
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void SampleByCompositionRejection(typename Backend::int     elementZ,
                                    typename Backend::double  energyIn,
                                    typename Backend::double& energyElectron,
                                    typename Backend::double& energyPositron,
                                    typename Backend::double& sinThetaElectron,
                                    typename Backend::double& sinThetaPositron);

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
  void
  SampleSinTheta(typename Backend::double energyElectron,
                 typename Backend::double energyPositron,
		 typename Backend::double& sinThetaElectron,
		 typename Backend::double& sinThetaPositron) const; 

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::double 
  TotalCrossSection(typename Backend::double energy,
                    typename Backend::double Z) const;

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
GUConversionBetheHeitler::InteractKernel(typename Backend::double  energyIn, 
                                         typename Backend::Index_t   zElement,
                                         typename Backend::double& energyElectron,
                                         typename Backend::double& energyPositron,
                                         typename Backend::double& sinThetaElectron,
                                         typename Backend::double& sinThetaPositron)
                                         const
{
  typedef typename Backend::Bool_t   Bool_t;
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
  
  double mininumE = electron_mass_c2;
  double deltaE = energyIn - mininumE;

  //electron energy
  double energyOut = mininumE + fAliasSampler->SampleX<Backend>(deltaE,probNA,
					        aliasInd,icol,fraction);


  double r1 = UniformRandom<Backend>(fRandomState,fThreadId);
  Bool_t condition = 0.5 > r1;

  MaskedAssign( condition, energyOut, &energyElectron);     
  MaskedAssign( condition, energyIn - energyOut, &energyPositron);     

  MaskedAssign(!condition, energyOut, &energyPositron);     
  MaskedAssign(!condition, energyIn - energyOut, &energyElectron);     
  
  SampleSinTheta<Backend>(energyElectron,energyPositron,
                          sinThetaElectron, sinThetaPositron);
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
void
GUConversionBetheHeitler::RotateAngle(typename Backend::double sinTheta,
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
void
GUConversionBetheHeitler::
SampleSinTheta(typename Backend::double energyElectron,
               typename Backend::double energyPositron,
	       typename Backend::double& sinThetaElectron,
	       typename Backend::double& sinThetaPositron) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::double double;

  //angles of the pair production (gamma -> e+e-)

  double u;
  const double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  double r1 =  UniformRandom<Backend>(fRandomState,fThreadId);
  Bool_t condition = 9./(9. + d) > r1;
  MaskedAssign( condition, -log( UniformRandom<Backend>(fRandomState,fThreadId)*
                       UniformRandom<Backend>(fRandomState,fThreadId))/a1, &u );
  MaskedAssign(!condition, -log( UniformRandom<Backend>(fRandomState,fThreadId)*
                       UniformRandom<Backend>(fRandomState,fThreadId))/a2, &u );

  double TetEl = u*electron_mass_c2/energyElectron;
  double TetPo = u*electron_mass_c2/energyPositron;

  //sinTheta - just return theta instead!
  sinThetaElectron =  sin(TetEl);
  sinThetaPositron = -sin(TetPo);
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH 
void GUConversionBetheHeitler::
SampleByCompositionRejection(typename Backend::int     elementZ,
                             typename Backend::double  energyIn,
                             typename Backend::double& energyElectron,
                             typename Backend::double& energyPositron,
                             typename Backend::double& sinThetaElectron,
                             typename Backend::double& sinThetaPositron)
{
  typedef typename Backend::double double;

  double epsil ;
  double epsil0 = electron_mass_c2/energyIn ;
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
  //  const double twopi = 2.*3.14159265358979323846;

  double TetEl = u*electron_mass_c2/energyElectron;
  double TetPo = u*electron_mass_c2/energyPositron;

  //return sinTheta
  sinThetaElectron = sin(TetEl);
  sinThetaPositron =-sin(TetPo);

  // direction
  //  double Phi  = twopi * UniformRandom<Backend>(fRandomState,fThreadId);
  //

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
  InteractKernel<Backend>(energyIn, targetElement, energyElectron, energyPositron, 
                          sinThetaElectron, sinThetaPositron);
  
  //update final states of the primary and store the secondary
  ConvertXtoFinalState<Backend>(energyIn,energyElectron, energyPositron,
                                sinThetaElectron,sinThetaPositron,
                                inProjectile,outSecondary);
}
  
#ifndef VECPHYS_NVCC
template <typename Backend>
void GUConversionBetheHeitler::Interact(GUTrack_v& inProjectile,  
                                        const int *targetElements,
                                        GUTrack_v& outSecondary ) const
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
    double sinThetaElectron;
    double energyElectron;
    double sinThetaPositron;
    double energyPositron;
    Index_t  zElement(targetElements[ibase]);

    InteractKernel<Backend>(energyIn, zElement, energyElectron, energyPositron,
                            sinThetaElectron,sinThetaPositron);

    //need to rotate the angle with respect to the line of flight
    double invp = 1./energyIn;
    double xhat = px*invp;
    double yhat = py*invp;
    double zhat = pz*invp;

    double uhat = 0.;
    double vhat = 0.;
    double what = 0.;


    // Kill the primary photon: use it as a placeholder for the secondary 
    // positron temporarily - change (charge, pid) accordingly 

    RotateAngle<Backend>(sinThetaPositron,xhat,yhat,zhat,uhat,vhat,what);

    energyPositron.store(&inProjectile.E[ibase]);
    double pxFinal, pyFinal, pzFinal;
     
    pxFinal= energyPositron*uhat;
    pyFinal= energyPositron*vhat;
    pzFinal= energyPositron*what;
    pxFinal.store(&inProjectile.px[ibase]);
    pyFinal.store(&inProjectile.py[ibase]);
    pzFinal.store(&inProjectile.pz[ibase]);

    // create Secondary
    RotateAngle<Backend>(sinThetaElectron,xhat,yhat,zhat,uhat,vhat,what);

    double secE = energyElectron; 
    double pxSec= secE*uhat;
    double pySec= secE*vhat;
    double pzSec= secE*what;

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

  Precision energyElectron;
  Precision energyPositron;
  Precision sinThetaElectron;
  Precision sinThetaPositron;

  SampleByCompositionRejection<Backend>(targetElement,energyIn,
                                        energyElectron,energyPositron,
					sinThetaElectron,sinThetaPositron);

  ConvertXtoFinalState<Backend>(energyIn,energyElectron,energyPositron,
                                sinThetaElectron,sinThetaPositron,
                                inProjectile,outSecondary);
  
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::double 
GUConversionBetheHeitler::
TotalCrossSection(typename Backend::double energy,
                  typename Backend::double Z) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::double double;

  double sigma = 0.;

  if ( Z < 0.9 || energy <= 2.0*electron_mass_c2 ) { return sigma; }

  double energySave = energy;

  //gamma energyLimit = 1.5*MeV
  double energyLimit = 1.5*MeV;
  Bool_t condition = energy < energyLimit;
  MaskedAssign( condition, energyLimit, &energy );
  
  double X = log(energy/electron_mass_c2);
  double X2 = X*X;
  double X3 =X2*X;
  double X4 =X3*X;
  double X5 =X4*X;

  //put coff's to a constant header
  double a0= 8.7842e+2*microbarn; 
  double a1=-1.9625e+3*microbarn; 
  double a2= 1.2949e+3*microbarn;
  double a3=-2.0028e+2*microbarn; 
  double a4= 1.2575e+1*microbarn; 
  double a5=-2.8333e-1*microbarn;
  
  double b0=-1.0342e+1*microbarn; 
  double b1= 1.7692e+1*microbarn; 
  double b2=-8.2381   *microbarn;
  double b3= 1.3063   *microbarn; 
  double b4=-9.0815e-2*microbarn; 
  double b5= 2.3586e-3*microbarn;
  
  double c0=-4.5263e+2*microbarn; 
  double c1= 1.1161e+3*microbarn; 
  double c2=-8.6749e+2*microbarn;
  double c3= 2.1773e+2*microbarn; 
  double c4=-2.0467e+1*microbarn; 
  double c5= 6.5372e-1*microbarn;

  double F1 = a0 + a1*X + a2*X2 + a3*X3 + a4*X4 + a5*X5;
  double F2 = b0 + b1*X + b2*X2 + b3*X3 + b4*X4 + b5*X5;
  double F3 = c0 + c1*X + c2*X2 + c3*X3 + c4*X4 + c5*X5;     

  sigma = (Z + 1.)*(F1*Z + F2*Z*Z + F3);
  Bool_t done = energySave < energyLimit;

  if(Any(done)) {
    X = (energySave - 2.*electron_mass_c2)/(energyLimit - 2.*electron_mass_c2);
    double tmpsigma = sigma*X*X;
    MaskedAssign( done, tmpsigma, &sigma );
  }

  Bool_t check = sigma < 0.;
  MaskedAssign( check, 0., &sigma );

  return sigma;
}

} // end namespace impl
} // end namespace vecphys

#endif
