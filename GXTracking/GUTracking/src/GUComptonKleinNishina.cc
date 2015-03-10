#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "GUComptonKleinNishina.h"
#include <iostream>

#include "backend/Backend.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_BOTH 
GUComptonKleinNishina::GUComptonKleinNishina(Random_t* states, int threadId) 
  :
  fRandomState(states), fThreadId(threadId),
  fMinX(1.), fMaxX(1001.), fDeltaX(0.1), 
  fMinY(0.), fMaxY(1001.), fDeltaY(0.1),
  fNrow(100), fNcol(100) 
{
  //replace hard coded numbers by default constants

  fAliasSampler = new GUAliasSampler(fRandomState,fThreadId,
                                     10,fMinX,fMaxX,fNrow,fNcol);

  //eventually arguments of BuildTable should be replaced by members of *this
  //and dropped from the function signature. Same for BuildPdfTable
  BuildTable(10,fMinX,fMaxX,fNrow,fNcol);   
}

VECPHYS_CUDA_HEADER_BOTH 
GUComptonKleinNishina::GUComptonKleinNishina(Random_t* states, int threadId,
                                             GUAliasSampler* sampler) 
  :
  fRandomState(states), fThreadId(threadId),
  fMinX(1.), fMaxX(1001.), fDeltaX(0.1), 
  fMinY(0.), fMaxY(1001.), fDeltaY(0.1),
  fNrow(100), fNcol(100) 
{
  //replace hard coded numbers by default constants
  fAliasSampler = sampler;
}

//need another Ctor with setable parameters

VECPHYS_CUDA_HEADER_BOTH 
GUComptonKleinNishina::~GUComptonKleinNishina() 
{
  if(fAliasSampler) delete fAliasSampler;
}

#ifdef VECPHYS_VC 

void 
GUComptonKleinNishina::Interact( GUTrack_v& inProjectile,    // In/Out
          const int *targetElements,  // Number equal to num of tracks
          GUTrack_v& outSecondary    // Empty vector for secondaries
          ) const
{
  Vc::double_v energyIn;
  Vc::double_v deltaE; //temporary - this should be dy in BuildPdfTable
  Vc::Vector<Precision> index;
  Vc::Vector<Precision> icol;
  Vc::double_v fraction;

  Vc::double_v px;
  Vc::double_v py;
  Vc::double_v pz;

  for(int i=0; i < inProjectile.numTracks/Vc::double_v::Size ; ++i) {

    //gather
    // loads energies into a VC-"register" type called energyIn
    for(int j = 0; j < Vc::double_v::Size ; ++j) {
      energyIn[j] = inProjectile.E[ i*Vc::double_v::Size + j];
      deltaE[j] = energyIn[j] - energyIn[j]/(1+2.0*energyIn[j]*inv_electron_mass_c2);

      px[j] =  inProjectile.px[ i*Vc::double_v::Size + j];
      py[j] =  inProjectile.py[ i*Vc::double_v::Size + j];
      pz[j] =  inProjectile.pz[ i*Vc::double_v::Size + j];
    }

    fAliasSampler->SampleBin<kVc>(energyIn,index,icol,fraction);

    Vc::Vector<Precision> probNA;   // Non-alias probability
    Vc::Vector<Precision> aliasInd; // This is really an integer -- could be Index_t !?

    //gather for alias table lookups
    fAliasSampler->GatherAlias<kVc>(index,probNA,aliasInd);

    Vc::double_v energyOut = fAliasSampler->SampleX<kVc>(deltaE,probNA,aliasInd,
						   icol,fraction);

    //TODO: write back result energyOut somewhere

    //store only secondary energy for now
    //evaluate the scattered angle based on xV

    Vc::double_v sinTheta = SampleSinTheta<kVc>(energyIn,energyOut);

    //need to rotate the angle with respect to the line of flight
    Vc::double_v invp = 1./energyIn;
    Vc::double_v xhat = px*invp;
    Vc::double_v yhat = py*invp;
    Vc::double_v zhat = pz*invp;

    Vc::double_v uhat = 0.;
    Vc::double_v vhat = 0.;
    Vc::double_v what = 0.;

    RotateAngle<kVc>(sinTheta,xhat,yhat,zhat,uhat,vhat,what);

    //scatter 
    for(int j = 0; j < Vc::double_v::Size ; ++j) {
      int it = i*Vc::double_v::Size + j;

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

VECPHYS_CUDA_HEADER_BOTH void 
GUComptonKleinNishina::BuildTable( int Z, 
                                   const double xmin, 
                                   const double xmax,
                                   const int nrow,
                                   const int ncol)
{
  //for now, the model does not own pdf.  Otherwise, pdf should be the 
  //data member of *this and set its point to the fpdf of fAliasSampler 
  double *pdf = new double [nrow*ncol];

  BuildPdfTable(Z,xmin,xmax,nrow,ncol,pdf); 
  fAliasSampler->BuildAliasTables(nrow,ncol,pdf);

  delete [] pdf;
}

VECPHYS_CUDA_HEADER_BOTH void 
GUComptonKleinNishina::BuildPdfTable(int Z, 
                                     const double xmin, 
                                     const double xmax,
                                     const int nrow,
                                     const int ncol,
                                     double *p
                                     )
{
  // Build the probability density function (KleinNishina pdf) 
  // in the energy randge [xmin,xmax] with an equal bin size
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  //           xmin (miminum energy)
  //           xmax (maxinum energy)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  //build pdf  
  double dx = (xmax - xmin)/nrow;
  double xo =  xmin + 0.5*dx;


  for(int i = 0; i < nrow ; ++i) {
    //for each input energy bin
    double x = xo + dx*i;

    double ymin = x/(1+2.0*x*inv_electron_mass_c2);
    double dy = (x - ymin)/(ncol-1);
    double yo = ymin + 0.5*dy;
  
    double sum = 0.;

    for(int j = 0; j < ncol ; ++j) {
      //for each output energy bin
      double y = yo + dy*j;
      double xsec = CalculateDiffCrossSection(0,x,y);
      p[i*ncol+j] = xsec;
      sum += xsec;
    }

    //normalization
    sum = 1.0/sum;

    for(int j = 0; j < ncol ; ++j) {
      p[i*ncol+j] *= sum;
    }
  }
}

// function implementing the cross section for KleinNishina
// TODO: need to get electron properties from somewhere

VECPHYS_CUDA_HEADER_BOTH double 
GUComptonKleinNishina::CalculateDiffCrossSection( int Zelement, 
                                                  double energy0, 
                                                  double energy1 ) const
{
  // based on Geant4 : KleinNishinaCompton
  // input  : energy0 (incomming photon energy)
  //          energy1 (scattered photon energy)
  // output : dsigma  (differential cross section) 

  double E0_m = energy0*inv_electron_mass_c2 ;
  double epsilon = energy1/energy0;

  double onecost = (1.- epsilon)/(epsilon*E0_m);
  double sint2   = onecost*(2.-onecost);
  double greject = 1. - epsilon*sint2/(1.+ epsilon*epsilon);
  double dsigma = (epsilon + 1./epsilon)*greject;

  return dsigma;
}

} // end namespace impl
} // end namespace vecphys
