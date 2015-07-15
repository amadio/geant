#ifndef EmModelBase_H
#define EmModelBase_H 1

#include "backend/Backend.h"

#include "GUConstants.h"
#include "GUTrack.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

template <class EmModel>
class EmModelBase {

public:

  VECPHYS_CUDA_HEADER_BOTH 
  EmModelBase(EmModel * model, Random_t* states, int tid) 
    : fModel(model), fRandomState(states), fThreadId(tid) {};

  //scalar
  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH 
  void Interact(GUTrack& projectile,   
                int      targetElement,
                GUTrack& secondary );

  //vector
#ifndef VECPHYS_NVCC
  template <typename Backend>
  void Interact(GUTrack_v& inProjectile,  
                const int *targetElements,
                GUTrack_v& outSecondaryV );     
#endif

  //validation 
  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void InteractG4(GUTrack& inProjectile,
                  int      targetElement,
                  GUTrack& outSecondary );

private:

  // Implementation methods
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void RotateAngle(typename Backend::Double_t sinTheta,
                   typename Backend::Double_t xhat,
                   typename Backend::Double_t yhat,
                   typename Backend::Double_t zhat,
                   typename Backend::Double_t &xr,
                   typename Backend::Double_t &yr,
                   typename Backend::Double_t &zr);

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void ConvertXtoFinalState(double energyIn, 
                            double energyOut, 
                            double sinTheta, 
                            GUTrack& primary, 
                            GUTrack& secondary);

  //data members

private:
  EmModel *const fModel; // to keep a pointer of implementation class

protected:
  Random_t* fRandomState;
  int       fThreadId;
};

//Implementation
template <class EmModel>
template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void EmModelBase<EmModel>::Interact(GUTrack& inProjectile,
			            int      targetElement,
			            GUTrack& outSecondary ) 
{
  double energyIn= inProjectile.E;
  double energyOut, sinTheta;

  fModel->EmModel::template InteractKernel<Backend>(energyIn,targetElement,energyOut,sinTheta);

  //update final states of the primary and store the secondary
  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
}
  
#ifndef VECPHYS_NVCC
template <class EmModel>
template <typename Backend>
void EmModelBase<EmModel>::Interact( GUTrack_v& inProjectile,  
			             const int *targetElements,
			             GUTrack_v& outSecondary) 
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::Double_t Double_t;

  for(int j = 0; j < inProjectile.numTracks  ; ++j) {
    assert( (targetElements[j] > 0)  && (targetElements[j] <= maximumZ ) );
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

    fModel->EmModel::template InteractKernel<Backend>(energyIn,zElement,energyOut,sinTheta);

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

template <class EmModel>
template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH 
void EmModelBase<EmModel>::InteractG4(GUTrack& inProjectile,
                                      int      targetElement,
                                      GUTrack& outSecondary )
{

  Precision energyIn = inProjectile.E;
  Precision energyOut;
  Precision sinTheta;

  fModel->EmModel::SampleByCompositionRejection(energyIn,energyOut,sinTheta);

  //update final states of the primary and store the secondary
  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
}

template <class EmModel>
template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void
EmModelBase<EmModel>::RotateAngle(typename Backend::Double_t sinTheta,
                                  typename Backend::Double_t xhat,
                                  typename Backend::Double_t yhat,
                                  typename Backend::Double_t zhat,
                                  typename Backend::Double_t &xr,
                                  typename Backend::Double_t &yr,
                                  typename Backend::Double_t &zr)
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

} // end namespace impl
} // end namespace vecphys

#endif
