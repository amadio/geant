#include "GUAliasSampler.h"
#include "GUComptonKleinNishina.h"
#include <iostream>

#include "backend/Backend.h"
#include "backend/vc/Backend.h"

using namespace VECGEOM_NAMESPACE;

FQUALIFIER 
GUComptonKleinNishina::GUComptonKleinNishina() 
  :
  fRandomState(0), fThreadId(-1),
  fMinX(1.0), fMaxX(1001), fDeltaX(0.1),
  fMinY(0.), fMaxY(0.), fDeltaY(0.),
  fNrow(100), fNcol(100)
{
  fAliasSampler = new GUAliasSampler(10,fMinX,fMaxX,fNrow,fNcol);

  //eventually arguments of BuildTable should be replaced by members of *this
  //and dropped from the function signature. Same for BuildPdfTable
  BuildTable(10,fMinX,fMaxX,fNrow,fNcol);   

}

//need another Ctor with setable parameters

FQUALIFIER 
GUComptonKleinNishina::~GUComptonKleinNishina() 
{
  if(fAliasSampler) delete fAliasSampler;
}

FQUALIFIER 
void GUComptonKleinNishina::Interact(GUTrack& projectile,
                                     int      targetElement,
                                     GUTrack* secondary ) const
{
  //  double delta = 0.1;
  //  double x = fAliasSampler->Sample<vecgeom::kScalar>(projectile.E,delta);
}


FQUALIFIER void 
GUComptonKleinNishina::Interact( GUTrack_v& inProjectile,    // In/Out
          const int *targetElements,  // Number equal to num of tracks
          GUTrack_v* outSecondaryV    // Empty vector for secondaries
          ) const
{
  Vc::double_v EVc;
  Vc::double_v deltaVc; //temporary - this should be dy in BuildPdfTable

  for(int i=0; i < inProjectile.numTracks/Vc::double_v::Size ; ++i) {

    // loads energies into a VC-"register" type called EVc
    // Vc::double_v EVc( &inProjectile.E[i*Vc::double_v::Size]); 

    //gather
    for(int j = 0; j < Vc::double_v::Size ; ++j) {
      EVc[j] = inProjectile.E[ i*Vc::double_v::Size + j];
      deltaVc[j] =  EVc[j] - EVc[j]/(1+2.0*EVc[j]*inv_electron_mass_c2);
    }

    Vc::double_v xVc = fAliasSampler->Sample<kVc>(EVc,deltaVc);

    //TODO: write back result xVc somewhere
    // xVc.store(/* some address*/);

    //store only secondary energy for now
    //evaluate the scattered angle based on xV
    Vc::double_v angleVc = SampleSinTheta<kVc>(EVc,xVc);

    //need to rotate the angle with respect to the line of flight

    //scatter 
    for(int j = 0; j < Vc::double_v::Size ; ++j) {
      outSecondaryV->E[ i*Vc::double_v::Size + j] = xVc[j]; 

      //fill also (x,y,z) and (px,py,pz), q and etc
    }
  }
}    

FQUALIFIER void 
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

FQUALIFIER void 
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

FQUALIFIER double 
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

FQUALIFIER void
GUComptonKleinNishina::SampleByCompositionRejection(double energyIn,
						    double& energyOut,
						    double& sinTheta)
{
  double epsilon, epsilonsq, onecost, sint2, greject ;
  
  double E0_m = energyIn*inv_electron_mass_c2 ;
  double eps0       = 1./(1. + 2.*E0_m);
  double epsilon0sq = eps0*eps0;
  double alpha1     = - log(eps0);
  double alpha2     = 0.5*(1.- epsilon0sq);
  
  do {
    if( alpha1/(alpha1+alpha2) > GUUniformRand(0, -1) ) {
      epsilon   = exp(-alpha1*GUUniformRand(0, -1));
      epsilonsq = epsilon*epsilon; 
    } 
    else {
      epsilonsq = epsilon0sq+(1.- epsilon0sq)*GUUniformRand(0,-1)
;
      epsilon   = sqrt(epsilonsq);
    }
    
    onecost = (1.- epsilon)/(epsilon*E0_m);
    sint2   = onecost*(2.-onecost);
    greject = 1. - epsilon*sint2/(1.+ epsilonsq);
    
  } while (greject < GUUniformRand(0, -1));
  
  energyOut = epsilon*energyIn;
  sinTheta = (sint2 < 0.0) ? 0.0 : sqrt(sint2);
}
