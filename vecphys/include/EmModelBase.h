#ifndef EmModelBase_H
#define EmModelBase_H

#include "base/Global.h"
#include "base/SystemOfUnits.h"

#include "GUAuxFunctions.h"    // Define sincos if needed
#include "GUConstants.h"

#include "GUTrack.h"
#include "GUAliasSampler.h"
#include "SamplingMethod.h"
#include "MaterialHandler.h"

#ifndef VECCORE_NVCC
#include <bitset>
#include <vector>
#endif

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

template <class EmModel>
class EmModelBase {

public:

  VECCORE_CUDA_HOST
  EmModelBase(Random_t* states, int tid);

  VECCORE_CUDA_HOST_DEVICE
  EmModelBase(Random_t* states, int tid, GUAliasSampler* sampler);

  VECCORE_CUDA_HOST_DEVICE
  ~EmModelBase();

  VECCORE_CUDA_HOST
  void Initialization();

  VECCORE_CUDA_HOST
  void BuildCrossSectionTable();

  VECCORE_CUDA_HOST
  void BuildAliasTable(bool atomicDependentModel = false);

  //scalar
  template <typename Backend>
  VECCORE_CUDA_HOST_DEVICE
  void AtomicCrossSection(GUTrack&  projectile,
                          const int targetElement,
                          double&   sigma);

  template <typename Backend>
  VECCORE_CUDA_HOST_DEVICE
  void Interact(GUTrack&  projectile,
                const int targetElement,
                GUTrack&  secondary );

  //vector
#ifndef VECCORE_NVCC
  template <typename Backend>
  void AtomicCrossSection(GUTrack_v& inProjectile,
                          const int* targetElements,
                          double*    sigma);

  template <typename Backend>
  void Interact(GUTrack_v& inProjectile,
                const int* targetElements,
                GUTrack_v& outSecondaryV);

  //temporary method for testing
  template <typename Backend>
  void InteractUnpack(GUTrack_v& inProjectile,
                      const int* targetElements,
                      GUTrack_v& outSecondary);

#endif

  //validation
  template <typename Backend>
  VECCORE_CUDA_HOST_DEVICE
  void AtomicCrossSectionG4(GUTrack&  inProjectile,
                            const int targetElement,
                            double&   sigma);

  template <typename Backend>
  VECCORE_CUDA_HOST_DEVICE
  void InteractG4(GUTrack&  inProjectile,
                  const int targetElement,
                  GUTrack&  outSecondary);

  VECCORE_CUDA_HOST_DEVICE
  GUAliasSampler* GetSampler() {return fAliasSampler;}

  VECCORE_CUDA_HOST_DEVICE
  void SetSampler(GUAliasSampler* sampler) { fAliasSampler = sampler ;}

  VECCORE_CUDA_HOST_DEVICE
  void SetSamplingMethod(SamplingMethod type) { fSampleType = type; }

  VECCORE_CUDA_HOST_DEVICE
  SamplingMethod GetSamplingMethod() { return fSampleType; }

protected:
  // Auxiliary methods
  VECCORE_CUDA_HOST_DEVICE
  void SetLowEnergyLimit(double lowLimit) { fLowEnergyLimit = lowLimit; }

  VECCORE_CUDA_HOST_DEVICE
  void SetHighEnergyLimit(double highLimit) { fHighEnergyLimit = highLimit; }

  VECCORE_CUDA_HOST_DEVICE double
  ComputeCoulombFactor(double fZeff);

protected:
  // Implementation methods
  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  void RotateAngle(typename Backend::Double_v sinTheta,
                   typename Backend::Double_v xhat,
                   typename Backend::Double_v yhat,
                   typename Backend::Double_v zhat,
                   typename Backend::Double_v &xr,
                   typename Backend::Double_v &yr,
                   typename Backend::Double_v &zr);

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  void ConvertXtoFinalState(double energyIn,
                            double energyOut,
                            double sinTheta,
                            GUTrack& primary,
                            GUTrack& secondary);

#ifndef VECCORE_NVCC
  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  void ConvertXtoFinalState(typename Backend::Double_v energyIn,
                            typename Backend::Double_v energyOut,
                            typename Backend::Double_v sinTheta,
                            int index,
                            GUTrack_v& primary,
                            GUTrack_v& secondary);

  //this inner template cannot be specialized unless template <class EmModel>
  //is also explicitly specialized
  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  void ConvertXtoFinalState_Scalar(typename Backend::Double_v energyIn,
                                   typename Backend::Double_v energyOut,
                                   typename Backend::Double_v sinTheta,
                                   int index,
                                   GUTrack_v& primary,
                                   GUTrack_v& secondary);
#endif

  //data members
protected:
  Random_t* fRandomState;
  int       fThreadId;

  bool      fAtomicDependentModel;
  double    fLowEnergyLimit;
  double    fHighEnergyLimit;

  //Sampling methods
  SamplingMethod  fSampleType;

  //Alias sampling Tables
  GUAliasSampler* fAliasSampler;
};

//Implementation
template <class EmModel>
VECCORE_CUDA_HOST
EmModelBase<EmModel>::EmModelBase(Random_t* states, int tid)
  : fRandomState(states), fThreadId(tid),
    fAtomicDependentModel(false),
    fLowEnergyLimit(0.1*keV),
    fHighEnergyLimit(1.0*TeV),
    fSampleType(kAlias),
    fAliasSampler(0)
{
}

template <class EmModel>
VECCORE_CUDA_HOST_DEVICE
EmModelBase<EmModel>::EmModelBase(Random_t* states, int tid, GUAliasSampler* sampler)
  : fRandomState(states), fThreadId(tid),
    fAtomicDependentModel(false),
    fLowEnergyLimit(0.1*keV),
    fHighEnergyLimit(1.0*TeV),
    fSampleType(kAlias)
{
  fAliasSampler = sampler;
}

template <class EmModel>
VECCORE_CUDA_HOST_DEVICE
EmModelBase<EmModel>::~EmModelBase()
{
  //  if(fAliasSampler) delete fAliasSampler;
}

template <class EmModel>
VECCORE_CUDA_HOST
void EmModelBase<EmModel>::BuildCrossSectionTable()
{
  //dummy interface for now
  for( int z= 1 ; z < maximumZ ; ++z)
  {
    static_cast<EmModel*>(this)->BuildCrossSectionTablePerAtom(z);
  }
}

template <class EmModel>
VECCORE_CUDA_HOST
  void EmModelBase<EmModel>::BuildAliasTable(bool atomicDependentModel)
{
  //size of the array for the alias table data
  size_t sizeOfTable = (fAliasSampler->GetNumEntries()+1)*fAliasSampler->GetSamplesPerEntry();
  double *pdf = new double [sizeOfTable];

  int z = -1;
  if(fAtomicDependentModel) {
    MaterialHandler* materialHander = MaterialHandler::Instance();
    int numberOfElements = materialHander->GetNumberOfElements();
    for( int i = 0 ; i < numberOfElements ; ++i) {
      int z = (materialHander->GetElementArray())[i];
      static_cast<EmModel*>(this)->BuildPdfTable(z,pdf);
      fAliasSampler->BuildAliasTable(z,pdf);
    }
  }
  else {
    z = 0 ; //Z=0 is convention for atomic independent models
    static_cast<EmModel*>(this)->BuildPdfTable(z,pdf);
    fAliasSampler->BuildAliasTable(z,pdf);
  }
  delete [] pdf;
}

template <class EmModel>
template <typename Backend>
VECCORE_CUDA_HOST_DEVICE
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
VECCORE_CUDA_HOST_DEVICE
void EmModelBase<EmModel>::Interact(GUTrack&  inProjectile,
                                    const int targetElement,
                                    GUTrack&  outSecondary )
{
  double energyIn = inProjectile.E;

  //check for the validity of energy
  if(energyIn < fLowEnergyLimit || energyIn > fHighEnergyLimit) return;

  double energyOut =0;
  double sinTheta = 0;

  //eventually there will be no switch as we select one Kernel for each model
  switch(fSampleType) {
    case SamplingMethod::kAlias :
      static_cast<EmModel*>(this)-> template InteractKernel<Backend>(energyIn,targetElement,energyOut,sinTheta);
      break;
    case SamplingMethod::kRejection :
        static_cast<EmModel*>(this)-> template InteractKernelCR<Backend>(energyIn,targetElement,energyOut,sinTheta);
      break;
    case SamplingMethod::kUnpack :
      {
        bool status = false;
        do {
          static_cast<EmModel*>(this)-> template InteractKernelUnpack<Backend>(energyIn,
                                      targetElement,energyOut,sinTheta,status);
        } while ( status );
      }
    default :
      ;
  }

  //update final states of the primary and store the secondary
  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
}

#ifndef VECCORE_NVCC
template <class EmModel>
template <typename Backend>
void EmModelBase<EmModel>::AtomicCrossSection(GUTrack_v& inProjectile,
                                              const int* targetElements,
                                              double*    sigma)
{
  using Double_v = typename Backend::Double_v;

  for(int j = 0; j < inProjectile.numTracks  ; ++j) {
    assert( (targetElements[j] > 0)  && (targetElements[j] <= maximumZ ) );
  }

  int ibase= 0;
  int numChunks= (inProjectile.numTracks/VectorSize<Double_v>());

  for(int i=0; i < numChunks ; ++i) {
    Double_v energyIn(&inProjectile.E[ibase]);
    Index_v<Double_v>  zElement(targetElements[ibase]);

    Double_v sigmaOut = static_cast<EmModel*>(this)-> template CrossSectionKernel<Backend>(energyIn,zElement);

    sigmaOut.store(&sigma[ibase]);
    ibase+= VectorSize<Double_v>();
  }

  //leftover - do scalar
  for(int i = numChunks*VectorSize<Double_v>() ; i < inProjectile.numTracks ; ++i) {
    sigma[i] = static_cast<EmModel*>(this)-> template CrossSectionKernel<backend::Scalar>(inProjectile.E[i],targetElements[i]);
  }
}

template <class EmModel>
template <typename Backend>
void EmModelBase<EmModel>::Interact(GUTrack_v& inProjectile,
                                    const int* targetElements,
                                    GUTrack_v& outSecondary)
{
  //check for the validity of energy
  int nTracks = inProjectile.numTracks;

  // this inclusive check may be redundant as this model/process should not be
  // selected if energy of the track is outside the valid energy region
  //  if(inProjectile.E[0]         < fLowEnergyLimit ||
  //     inProjectile.E[nTracks-1] > fHighEnergyLimit) return;

  using Double_v = typename Backend::Double_v;

  for(int j = 0; j < nTracks  ; ++j) {
    assert( (targetElements[j] > 0)  && (targetElements[j] <= maximumZ ) );
  }

  int ibase= 0;
  int numChunks= (nTracks/VectorSize<Double_v>());

  for(int i= 0; i < numChunks ; ++i) {
    Index_v<Double_v>  zElement(targetElements[ibase]);
    Double_v energyIn(&inProjectile.E[ibase]);

    Double_v sinTheta(0.);
    Double_v energyOut;

    //eventually there will be no switch as we select one Kernel for each model
    switch(fSampleType) {
      case SamplingMethod::kAlias :
        static_cast<EmModel*>(this)-> template InteractKernel<Backend>(energyIn,zElement,energyOut,sinTheta);
        break;
      case SamplingMethod::kRejection :
        static_cast<EmModel*>(this)-> template InteractKernelCR<Backend>(energyIn,zElement,energyOut,sinTheta);
        break;
      case SamplingMethod::kUnpack :
        ; //dummy - see InteractUnpack for Vc
        break;
      default :
        ;
    }

    ConvertXtoFinalState<Backend>(energyIn, energyOut, sinTheta, ibase, inProjectile, outSecondary);
    ibase+= VectorSize<Double_v>();
  }

  //leftover - do scalar (temporary)
  for(int i = numChunks*VectorSize<Double_v>() ; i < nTracks ; ++i) {

    double senergyIn= inProjectile.E[i];
    double senergyOut, ssinTheta;

    static_cast<EmModel*>(this)-> template InteractKernel<backend::Scalar>(senergyIn,targetElements[i],senergyOut,ssinTheta);
    ConvertXtoFinalState_Scalar<backend::Scalar>(senergyIn, senergyOut, ssinTheta, i, inProjectile, outSecondary);
  }
}

//this is temporary implementation for the purpose of testing and will be removed
//from the code once performance and validation studies are done
template <class EmModel>
template <typename Backend>
void EmModelBase<EmModel>::InteractUnpack(GUTrack_v& inProjectile,
                                          const int* targetElements,
                                          GUTrack_v& outSecondary)
{
  //check for the validity of energy
  int sizeOfInputTracks = inProjectile.numTracks;
  if(inProjectile.E[0]                   < fLowEnergyLimit ||
     inProjectile.E[sizeOfInputTracks-1] > fHighEnergyLimit) return;

  using Double_v = typename Backend::Double_v;

  for(int j = 0; j < sizeOfInputTracks  ; ++j) {
    assert( (targetElements[j] > 0)  && (targetElements[j] <= maximumZ ) );
  }

  int ibase= 0;
  int numChunks= sizeOfInputTracks/VectorSize<Double_v>();

  //creat a bit set to hold the status of sampling

  const int maxSizeOfTracks = 160000; // the maximum number of tracks - use ~2496*64 for now
  if(sizeOfInputTracks > maxSizeOfTracks ) {
    std::cout << "Size of input tracks > maxSizeOfTracks = " << maxSizeOfTracks << std::endl;
    exit(0);
  }
  std::bitset<maxSizeOfTracks> flag;

  //working arrays for the shuffling loop: using std::vector for now
  std::vector<double> wenergyIn;
  std::vector<double> tenergyIn;

  std::vector<double> wenergyOut;
  std::vector<double> tenergyOut;

  std::vector<double> wsinTheta;
  std::vector<double> tsinTheta;

  std::vector<int> windex;
  std::vector<int> tindex;

  //initial copy
  for(int i = 0; i < sizeOfInputTracks ; ++i) {
    wenergyIn.push_back(inProjectile.E[i]);
    wenergyOut.push_back(0);
    wsinTheta.push_back(0);
    windex.push_back(i);
  }

  Mask_v<Double_v> status(false);

  //shuffling loop

  do {
    ibase= 0;
    flag.reset();

    //vectorization loop
    int currentSize = wenergyOut.size();
    numChunks = currentSize/VectorSize<Double_v>();

    for(int i= 0; i < numChunks ; ++i) {

      Double_v energyIn(&wenergyIn[ibase]);
      Index_v<Double_v>  zElement(targetElements[ibase]);

      Double_v energyOut;
      Double_v sinTheta;

      //kernel
      static_cast<EmModel*>(this)-> template InteractKernelUnpack<Backend>(energyIn,
                                             zElement,energyOut,sinTheta,status);

      //store energy and sinTheta of the secondary
      Double_v secE = energyIn - energyOut;

      secE.store(&wenergyOut[ibase]);
      sinTheta.store(&wsinTheta[ibase]);

      //set the status bit
      for(size_t j = 0; j < VectorSize<Double_v>() ; ++j) {
        flag.set(ibase+j,status[j]);
      }

      ibase+= VectorSize<Double_v>();
    }
    //end of vectorization loop

    //@@@syjun - need to add a cleanup task for the leftover array of which
    //size < VectorSize<Double_v>() - see the bottom of this method for a similar task

    //scatter
    for(int i = 0; i < currentSize ; ++i) {

      outSecondary.E[windex[i]]  = wenergyOut[i];
      outSecondary.px[windex[i]] = wsinTheta[i];

      tenergyIn.push_back(wenergyIn[i]);
      tenergyOut.push_back(wenergyOut[i]);
      tsinTheta.push_back(wsinTheta[i]);
      tindex.push_back(windex[i]);
    }

    //clear working arrays
    wenergyIn.clear();
    wenergyOut.clear();
    wsinTheta.clear();
    windex.clear();

    for(int i = 0; i < currentSize ; ++i) {
      if(flag.test(i)) {
        wenergyIn.push_back(tenergyIn[i]);
        wenergyOut.push_back(tenergyOut[i]);
        wsinTheta.push_back(tsinTheta[i]);
        windex.push_back(tindex[i]);
      }
    }

    //clear temporary arrays
    tenergyIn.clear();
    tenergyOut.clear();
    tsinTheta.clear();
    tindex.clear();

  } while ( wenergyOut.size() > 0 );
  //end of shuffling loop

  //reset ibase and number of chunks
  ibase= 0;
  numChunks= (sizeOfInputTracks/VectorSize<Double_v>());

  for(int i= 0; i < numChunks ; ++i) {
    Double_v energyIn(&inProjectile.E[ibase]);
    Double_v sinTheta(&outSecondary.px[ibase]);
    Double_v energyOut(&outSecondary.E[ibase]);

    ConvertXtoFinalState<Backend>(energyIn, energyOut, sinTheta, ibase, inProjectile, outSecondary);

    ibase+= VectorSize<Double_v>();
  }

  //leftover - do scalar (temporary)
  for(int i = numChunks*VectorSize<Double_v>() ; i < inProjectile.numTracks ; ++i) {

    double senergyIn= inProjectile.E[i];
    double senergyOut, ssinTheta;

    static_cast<EmModel*>(this)-> template InteractKernel<backend::Scalar>(senergyIn,targetElements[i],senergyOut,ssinTheta);
    ConvertXtoFinalState_Scalar<backend::Scalar>(senergyIn, senergyOut, ssinTheta, i, inProjectile, outSecondary);
  }
}

#endif

template <class EmModel>
template <typename Backend>
VECCORE_CUDA_HOST_DEVICE
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
VECCORE_CUDA_HOST_DEVICE
void EmModelBase<EmModel>::InteractG4(GUTrack&  inProjectile,
                                      const int targetElement,
                                      GUTrack&  outSecondary)
{

  Real_t energyIn = inProjectile.E;
  Real_t energyOut;
  Real_t sinTheta;

  static_cast<EmModel*>(this)->SampleByCompositionRejection(targetElement,energyIn,energyOut,sinTheta);

  //update final states of the primary and store the secondary
  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
}

template <class EmModel>
template<class Backend>
VECCORE_CUDA_HOST_DEVICE void
EmModelBase<EmModel>::RotateAngle(typename Backend::Double_v sinTheta,
                                  typename Backend::Double_v xhat,
                                  typename Backend::Double_v yhat,
                                  typename Backend::Double_v zhat,
                                  typename Backend::Double_v &xr,
                                  typename Backend::Double_v &yr,
                                  typename Backend::Double_v &zr)
{
  typedef typename Backend::Int_t    Int_t;
  using Double_v = typename Backend::Double_v;

  Double_v phi = UniformRandom<Backend>(fRandomState,Int_t(fThreadId));
  Double_v pt = xhat*xhat + yhat*yhat;

  Double_v cosphi, sinphi;
  sincos(phi, &sinphi, &cosphi);

  Double_v uhat = sinTheta*cosphi; // cos(phi);
  Double_v vhat = sinTheta*sinphi; // sin(phi);
  Double_v what = Sqrt((1.-sinTheta)*(1.+sinTheta));

  Mask_v<Double_v> positive = ( pt > 0. );
  Mask_v<Double_v> negativeZ = ( zhat < 0. );

  Double_v phat;

  MaskedAssign(positive, Sqrt(pt) , &phat);
  MaskedAssign(positive, (xhat*zhat*uhat - yhat*vhat)/phat + xhat*what , &xr);
  MaskedAssign(positive, (yhat*zhat*uhat - xhat*vhat)/phat + yhat*what , &yr);
  MaskedAssign(positive, -phat*uhat + zhat*what , &zr);

  MaskedAssign(negativeZ,-xhat, &xr);
  MaskedAssign(negativeZ, yhat, &yr);
  MaskedAssign(negativeZ,-zhat, &zr);

  MaskedAssign(!positive && !negativeZ, xhat, &xr);
  MaskedAssign(!positive && !negativeZ, yhat, &yr);
  MaskedAssign(!positive && !negativeZ, zhat, &zr);

  //mask operation???
 /* if(positive) {
    Double_v phat = Sqrt(pt);
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
  }*/
}

template <class EmModel>
template <typename Backend>
VECCORE_CUDA_HOST_DEVICE
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

#ifndef VECCORE_NVCC
template <class EmModel>
template <typename Backend>
VECCORE_CUDA_HOST_DEVICE void
EmModelBase<EmModel>::ConvertXtoFinalState(typename Backend::Double_v energyIn,
                                           typename Backend::Double_v energyOut,
                                           typename Backend::Double_v sinTheta,
                                           int ibase,
                                           GUTrack_v& primary,
                                           GUTrack_v& secondary) // const
{
    using Double_v = typename Backend::Double_v;

    //need to rotate the angle with respect to the line of flight
    Double_v px(&primary.px[ibase]);
    Double_v py(&primary.py[ibase]);
    Double_v pz(&primary.pz[ibase]);

    Double_v invp = 1./energyIn;
    Double_v xhat = px*invp;
    Double_v yhat = py*invp;
    Double_v zhat = pz*invp;

    Double_v uhat = 0.;
    Double_v vhat = 0.;
    Double_v what = 0.;

    RotateAngle<Backend>(sinTheta,xhat,yhat,zhat,uhat,vhat,what);

    // Update primary
    energyOut.store(&primary.E[ibase]);
    Double_v pxFinal, pyFinal, pzFinal;

    pxFinal= energyOut*uhat;
    pyFinal= energyOut*vhat;
    pzFinal= energyOut*what;
    pxFinal.store(&primary.px[ibase]);
    pyFinal.store(&primary.py[ibase]);
    pzFinal.store(&primary.pz[ibase]);

    // create Secondary
    Double_v secE = energyIn - energyOut;
    Double_v pxSec= secE*(xhat-uhat);
    Double_v pySec= secE*(yhat-vhat);
    Double_v pzSec= secE*(zhat-what);

    secE.store(&secondary.E[ibase]);
    pxSec.store(&secondary.px[ibase]);
    pySec.store(&secondary.py[ibase]);
    pzSec.store(&secondary.pz[ibase]);

  //fill other information

}

template<class EmModel>
template <typename Backend>
VECCORE_CUDA_HOST_DEVICE void
EmModelBase<EmModel>::ConvertXtoFinalState_Scalar(typename Backend::Double_v energyIn,
                                                  typename Backend::Double_v energyOut,
                                                  typename Backend::Double_v sinTheta,
                                                  int ibase,
                                                  GUTrack_v& primary,
                                                  GUTrack_v& secondary)
{
  typedef typename kScalar::Double_v Double_v;

  //need to rotate the angle with respect to the line of flight
  Double_v invp = 1./energyIn;
  Double_v xhat = primary.px[ibase]*invp;
  Double_v yhat = primary.py[ibase]*invp;
  Double_v zhat = primary.pz[ibase]*invp;

  Double_v uhat = 0.;
  Double_v vhat = 0.;
  Double_v what = 0.;

  RotateAngle<backend::Scalar>(sinTheta,xhat,yhat,zhat,uhat,vhat,what);

  // Update primary
  primary.E[ibase]  = energyOut;

  primary.px[ibase] = energyOut*uhat;
  primary.py[ibase] = energyOut*vhat;
  primary.pz[ibase] = energyOut*what;

  // create Secondary
  Double_v secE = energyIn - energyOut;
  secondary.E[ibase]  = secE;
  secondary.px[ibase] = secE*(xhat-uhat);
  secondary.py[ibase] = secE*(yhat-vhat);
  secondary.pz[ibase] = secE*(zhat-what);

  //fill other information
}
#endif

template <class EmModel>
VECCORE_CUDA_HOST_DEVICE double
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
