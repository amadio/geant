#include "GXKleinNishina.h"
#include "GPRandom.h"
#include "GPThreeVector.h"

FQUALIFIER
GXKleinNishina::GXKleinNishina(curandState* devState,
			       int threadId,
			       G4double minEnergy,
			       G4double maxEnergy,
			       G4int dimRow,
			       G4int dimCol,
			       G4double *pdfX,
			       G4double *pdfY) 
  :  GXVDiscreteSampling(devState,threadId,dimRow,dimCol,pdfX,pdfY)
{
  fThreadId = threadId;
  fDevState = devState; 

  //Sampling Tables
  fNrow = dimRow;
  fNcol = dimCol;

  fMinX = minEnergy;
  fMaxX = maxEnergy;
  fDeltaX = (fMaxX - fMinX)/fNrow;
}

FQUALIFIER
GXKleinNishina::GXKleinNishina(curandState* devState,
			       int threadId,
			       G4double minEnergy,
			       G4double maxEnergy,
			       G4int dimRow,
			       G4int dimCol,
			       G4int    *pdfA,
			       G4double *pdfX,
			       G4double *pdfY) 
  : GXVDiscreteSampling(devState,threadId,dimRow,dimCol,pdfA,pdfX,pdfY)
{
  fThreadId = threadId;
  fDevState = devState; 

  //Sampling Tables
  fNrow = dimRow;
  fNcol = dimCol;

  fMinX = minEnergy;
  fMaxX = maxEnergy;
  fDeltaX = (fMaxX - fMinX)/fNrow;

  fMinY = 0;
  fMaxY = 0;
  fDeltaY = 0;
}

FQUALIFIER GXKleinNishina::~GXKleinNishina()
{
  ;
}

FQUALIFIER void
GXKleinNishina::GetSampleParameters(G4double kineticEnergy,
				    G4int &irow, G4int &icol,
				    G4double &t) 
{
  irow = G4int((kineticEnergy - fMinX)/fDeltaX);
  G4double r1 = (fNcol-1)*GPUniformRand(fDevState, fThreadId);
  icol = int(r1);
  t = r1 - 1.0*icol;
  
  fMinY = kineticEnergy/(1+2.0*kineticEnergy/electron_mass_c2);
  fMaxY = kineticEnergy;
  fDeltaY = fMaxY - fMinY;
}

FQUALIFIER void
GXKleinNishina::SampleByInversePDF(G4double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(kineticEnergy,irow,icol,t);

  G4double yhat = InversePDF(irow,icol,fDeltaY,t);

  fGammaEnergy = fMinY + yhat;
  fGammaSinTheta = SampleDirection(kineticEnergy,fGammaEnergy);
}


FQUALIFIER void
GXKleinNishina::SampleByInversePDFLinearInterpolation(G4double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(kineticEnergy,irow,icol,t);

  //sample based on the inverse cumulataive pdf
  G4double yhat = InversePDFLinearInterpolation(irow,icol,fDeltaY,t);

  fGammaEnergy = fMinY + yhat;
  fGammaSinTheta = SampleDirection(kineticEnergy,fGammaEnergy);
  
}

FQUALIFIER void
GXKleinNishina::SampleByAlias(G4double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(kineticEnergy,irow,icol,t);

  //sample based on the alias method
  G4double yhat = Alias(irow,icol,fDeltaY,t);

  fGammaEnergy = fMinY + yhat;
  fGammaSinTheta = SampleDirection(kineticEnergy,fGammaEnergy);
}

FQUALIFIER void
GXKleinNishina::SampleByAliasLinearInterpolation(G4double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(kineticEnergy,irow,icol,t);

  //sample based on the alias method
  G4double yhat = AliasLinearInterpolation(irow,icol,fDeltaY,t);

  fGammaEnergy = fMinY + yhat;
  fGammaSinTheta = SampleDirection(kineticEnergy,fGammaEnergy);
}

FQUALIFIER G4double 
GXKleinNishina::SampleDirection(G4double energyIn, G4double energyOut) {

  //angle of the scatterred photon

  G4double epsilon = energyOut/energyIn;
  if(epsilon>1) epsilon = 1;
    
  G4double E0_m = energyIn/electron_mass_c2 ;
  G4double onecost = (1.- epsilon)/(epsilon*E0_m);
  G4double sint2   = onecost*(2.-onecost);
  G4double sinTheta = (sint2 < 0.0) ? 0.0 : sqrt(sint2); 

  return sinTheta;
}

FQUALIFIER void
GXKleinNishina::SampleByCompositionRejection(G4double gamEnergy0,
					     G4int& ntrial)
{
  G4double epsilon, epsilonsq, onecost, sint2, greject ;
  
  G4double E0_m = gamEnergy0/electron_mass_c2 ;
  G4double eps0       = 1./(1. + 2.*E0_m);
  G4double epsilon0sq = eps0*eps0;
  G4double alpha1     = - log(eps0);
  G4double alpha2     = 0.5*(1.- epsilon0sq);
  
  ntrial = 0 ;
  
  do {
    if( alpha1/(alpha1+alpha2) > GPUniformRand(fDevState, fThreadId) ) {
      epsilon   = exp(-alpha1*GPUniformRand(fDevState, fThreadId));
      epsilonsq = epsilon*epsilon; 
    } 
    else {
      epsilonsq = epsilon0sq+(1.- epsilon0sq)*GPUniformRand(fDevState,fThreadId);
      epsilon   = sqrt(epsilonsq);
    }
    
    onecost = (1.- epsilon)/(epsilon*E0_m);
    sint2   = onecost*(2.-onecost);
    greject = 1. - epsilon*sint2/(1.+ epsilonsq);
    
    ++ntrial;
    
  } while (greject < GPUniformRand(fDevState, fThreadId));
  
  fGammaEnergy = epsilon*gamEnergy0;
  fGammaSinTheta = (sint2 < 0.0) ? 0.0 : sqrt(sint2);
}

FQUALIFIER
G4double GXKleinNishina::GetSecondaryEnergy()
{
  return fGammaEnergy;
}

FQUALIFIER
G4double GXKleinNishina::GetSecondarySinTheta()
{
  return fGammaSinTheta;
}

FQUALIFIER
void GXKleinNishina::CreateSecondary(GXTrack* secondary,
				     GXTrack* track, 
				     G4double charge)
{
  G4double phi  = twopi*GPUniformRand(fDevState, fThreadId); 
  G4double sinTheta = fGammaSinTheta;
  G4double cosTheta = sqrt((1-sinTheta)*(1+sinTheta));
  
  GPThreeVector direction = GPThreeVector_create(sinTheta*cos(phi), 
                                                 sinTheta*sin(phi),
                                                 cosTheta);
  
  GPThreeVector refDirection = 
    GPThreeVector_unit(GPThreeVector_create(track->px,track->py,track->pz));
  
  //rotate along the refrefDirection;
  GPThreeVector_rotateUz(&direction,refDirection);
  
  //set secondary 
  G4double energy = fGammaEnergy;
  G4double p = sqrt(energy*(energy+2.0*electron_mass_c2*charge*charge));
  secondary->x  = track->x;
  secondary->y  = track->y;
  secondary->z  = track->z;
  secondary->s  = 0;
  secondary->px = p*direction.x;
  secondary->py = p*direction.y;
  secondary->pz = p*direction.z;
  secondary->E  = energy;
  secondary->q  = charge;
}
