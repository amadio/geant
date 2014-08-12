#include "GXMollerBhabha.h"
#include "GPRandom.h"
#include "GPThreeVector.h"

FQUALIFIER
GXMollerBhabha::GXMollerBhabha(curandState* devState,
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

  fMinY = 1.0;
  fMaxY = 0;
  fDeltaY = 0;

  fElectronFlag = true;
  fLimit = (fElectronFlag) ? 0.5 : 1.0;
}

FQUALIFIER
GXMollerBhabha::GXMollerBhabha(curandState* devState,
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

  fMinY = 1.0;
  fMaxY = 0;
  fDeltaY = 0;

  fElectronFlag = true;
  fLimit = (fElectronFlag) ? 0.5 : 1.0;
}

FQUALIFIER GXMollerBhabha::~GXMollerBhabha()
{
  ;
}

FQUALIFIER void
GXMollerBhabha::GetSampleParameters(G4double kineticEnergy,
				    G4int &irow, G4int &icol,
				    G4double &t) 
{
  irow = G4int((kineticEnergy - fMinX)/fDeltaX);
  G4double r1 = (fNcol-1)*GPUniformRand(fDevState, fThreadId);
  icol = int(r1);
  t = r1 - 1.0*icol;

  fMaxY = fLimit*kineticEnergy;
  fDeltaY = fMaxY - fMinY;
}

FQUALIFIER void
GXMollerBhabha::SampleByRandom(G4double kineticEnergy)
{
  fEnergy = fMinY + kineticEnergy*GPUniformRand(fDevState, fThreadId);
  fSinTheta = SampleDirection(kineticEnergy,fEnergy);
}

FQUALIFIER void
GXMollerBhabha::SampleByInversePDF(G4double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(kineticEnergy,irow,icol,t);

  G4double yhat = InversePDF(irow,icol,fDeltaY,t);

  fEnergy = fMinY + yhat;
  fSinTheta = SampleDirection(kineticEnergy,fEnergy);
  /*
  fEnergy = kineticEnergy*GPUniformRand(fDevState, fThreadId);
  fSinTheta = SampleDirection(kineticEnergy,fEnergy);
  */
}

FQUALIFIER void
GXMollerBhabha::SampleByInversePDFTexture(G4double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(kineticEnergy,irow,icol,t);

  G4double yhat = InversePDFTexture(irow,icol,fDeltaY,t);

  fEnergy = fMinY + yhat;
  fSinTheta = SampleDirection(kineticEnergy,fEnergy);
  /*
  fEnergy = kineticEnergy*GPUniformRand(fDevState, fThreadId);
  fSinTheta = SampleDirection(kineticEnergy,fEnergy);
  */
}

FQUALIFIER void
GXMollerBhabha::SampleByInversePDFLinearInterpolation(G4double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(kineticEnergy,irow,icol,t);

  //sample based on the inverse cumulataive pdf
  G4double yhat = InversePDFLinearInterpolation(irow,icol,fDeltaY,t);

  fEnergy = fMinY + yhat;
  fSinTheta = SampleDirection(kineticEnergy,fEnergy);
  
}

FQUALIFIER void
GXMollerBhabha::SampleByInversePDFTextureLinearInterpolation(G4double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(kineticEnergy,irow,icol,t);

  //sample based on the inverse cumulataive pdf
  G4double yhat = InversePDFTextureLinearInterpolation(irow,icol,fDeltaY,t);

  fEnergy = fMinY + yhat;
  fSinTheta = SampleDirection(kineticEnergy,fEnergy);
  
}

FQUALIFIER void
GXMollerBhabha::SampleByAlias(G4double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(kineticEnergy,irow,icol,t);

  //sample based on the alias method
  G4double yhat = Alias(irow,icol,fDeltaY,t);

  fEnergy = fMinY + yhat;
  fSinTheta = SampleDirection(kineticEnergy,fEnergy);
}

FQUALIFIER void
GXMollerBhabha::SampleByAliasLinearInterpolation(G4double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(kineticEnergy,irow,icol,t);

  //sample based on the alias method
  G4double yhat = AliasLinearInterpolation(irow,icol,fDeltaY,t);

  fEnergy = fMinY + yhat;
  fSinTheta = SampleDirection(kineticEnergy,fEnergy);
}

FQUALIFIER G4double 
GXMollerBhabha::SampleDirection(G4double energyIn, G4double energyOut) {

  //angle of the scatterred electron
  G4double energy = energyIn + electron_mass_c2;
  G4double totalMomentum = sqrt(energyIn*(energy + electron_mass_c2));

  G4double deltaMomentum = sqrt(energyOut * (energyOut + 2.0*electron_mass_c2));
  G4double cost =  energyOut * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum);
  G4double sint = (1.0 - cost)*(1. + cost);
  G4double sinTheta =  (sint < 0.0) ? 0.0 : sqrt(sint);

  return sinTheta;
}

FQUALIFIER void
GXMollerBhabha::SampleByCompositionRejection(G4double kineticEnergy,
					     G4int& counter)
{
  //cut energy
  G4double tmin = 1.0;  
  G4double tmax = fLimit*kineticEnergy; 

  if(fMaxX < tmax) { tmax = fMaxX; }
  if(tmin >= tmax) { return; }

  G4double energy = kineticEnergy + electron_mass_c2;
  G4double totalMomentum = sqrt(kineticEnergy*(energy + electron_mass_c2));

  G4double xmin   = tmin/kineticEnergy;
  G4double xmax   = tmax/kineticEnergy;
  G4double gam    = energy/electron_mass_c2;
  G4double gamma2 = gam*gam;
  G4double beta2  = 1.0 - 1.0/gamma2;
  G4double x, z, q, grej;

  counter = 0 ;

  if (fElectronFlag) {
    //Moller (e-e-) scattering
    G4double gg = (2.0*gam - 1.0)/gamma2;
    G4double y = 1.0 - xmax;
    grej = 1.0 - gg*xmax + xmax*xmax*(1.0 - gg + (1.0 - gg*y)/(y*y));
  
    do {
      q = GPUniformRand(fDevState, fThreadId);
      x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
      y = 1.0 - x;
      z = 1.0 - gg*x + x*x*(1.0 - gg + (1.0 - gg*y)/(y*y));
      ++counter;    
    } while(grej * GPUniformRand(fDevState, fThreadId) > z);
  } else {
    //Bhabha (e+e-) scattering
    G4double y   = 1.0/(1.0 + gam);
    G4double y2  = y*y;
    G4double y12 = 1.0 - 2.0*y;
    G4double b1  = 2.0 - y2;
    G4double b2  = y12*(3.0 + y2);
    G4double y122= y12*y12;
    G4double b4  = y122*y12;
    G4double b3  = b4 + y122;

    y    = xmax*xmax;
    grej = 1.0 + (y*y*b4 - xmin*xmin*xmin*b3 + y*b2 - xmin*b1)*beta2; 

    do {
      q = GPUniformRand(fDevState, fThreadId);
      x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
      y = x*x;
      z = 1.0 + (y*y*b4 - x*y*b3 + y*b2 - x*b1)*beta2; 
      ++counter;    
    } while(grej * GPUniformRand(fDevState, fThreadId) > z);
  }

  fEnergy = x * kineticEnergy;
  G4double deltaMomentum = sqrt(fEnergy * (fEnergy + 2.0*electron_mass_c2));
  G4double cost = fEnergy * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum);
  G4double sint = (1.0 - cost)*(1. + cost);

  fSinTheta =  (sint < 0.0) ? 0.0 : sqrt(sint);
}

FQUALIFIER void
GXMollerBhabha::SampleByAverageTrials(G4int ntrials,
				      G4double kineticEnergy,
				      G4int& counter)
{
  //cut energy
  G4double tmin = 1.0;  
  G4double tmax = fLimit*kineticEnergy; 

  if(fMaxX < tmax) { tmax = fMaxX; }
  if(tmin >= tmax) { return; }

  G4double energy = kineticEnergy + electron_mass_c2;
  G4double totalMomentum = sqrt(kineticEnergy*(energy + electron_mass_c2));

  G4double xmin   = tmin/kineticEnergy;
  G4double xmax   = tmax/kineticEnergy;
  G4double gam    = energy/electron_mass_c2;
  G4double gamma2 = gam*gam;
  G4double beta2  = 1.0 - 1.0/gamma2;
  G4double x, z, q, grej;

  counter = 0 ;

  if (fElectronFlag) {
    //Moller (e-e-) scattering
    G4double gg = (2.0*gam - 1.0)/gamma2;
    G4double y = 1.0 - xmax;
    grej = 1.0 - gg*xmax + xmax*xmax*(1.0 - gg + (1.0 - gg*y)/(y*y));
  
    for(G4int i = 0 ; i < ntrials ; ++i) {
      q = GPUniformRand(fDevState, fThreadId);
      x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
      y = 1.0 - x;
      z = 1.0 - gg*x + x*x*(1.0 - gg + (1.0 - gg*y)/(y*y));
      //dummy operation to ensure all calculation step
      grej += z;
      ++counter;    
    }
  } else {
    //Bhabha (e+e-) scattering
    G4double y   = 1.0/(1.0 + gam);
    G4double y2  = y*y;
    G4double y12 = 1.0 - 2.0*y;
    G4double b1  = 2.0 - y2;
    G4double b2  = y12*(3.0 + y2);
    G4double y122= y12*y12;
    G4double b4  = y122*y12;
    G4double b3  = b4 + y122;

    y    = xmax*xmax;
    grej = 1.0 + (y*y*b4 - xmin*xmin*xmin*b3 + y*b2 - xmin*b1)*beta2; 

    for(G4int i = 0 ; i < ntrials ; ++i) {
      q = GPUniformRand(fDevState, fThreadId);
      x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
      y = x*x;
      z = 1.0 + (y*y*b4 - x*y*b3 + y*b2 - x*b1)*beta2; 
      //dummy operation to ensure all calculation steps
      grej += z;
      ++counter;    
    }
  }

  //only for testing purpose - simulation result will be wrong !!!
  //  fEnergy = x * kineticEnergy;
  fEnergy = x * kineticEnergy + grej;

  G4double deltaMomentum = sqrt(fEnergy * (fEnergy + 2.0*electron_mass_c2));
  G4double cost = fEnergy * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum);
  G4double sint = (1.0 - cost)*(1. + cost);

  fSinTheta =  (sint < 0.0) ? 0.0 : sqrt(sint);
}

FQUALIFIER void 
GXMollerBhabha::SetElectronFlag(bool flag)
{
  fElectronFlag = flag;
  fLimit = (fElectronFlag) ? 0.5 : 1.0;
}

FQUALIFIER
G4double GXMollerBhabha::GetSecondaryEnergy()
{
  return fEnergy;
}

FQUALIFIER
G4double GXMollerBhabha::GetSecondarySinTheta()
{
  return fSinTheta;
}

FQUALIFIER
void GXMollerBhabha::CreateSecondary(GXTrack* secondary,
				     GXTrack* track, 
				     G4double charge)
{
  G4double phi  = twopi*GPUniformRand(fDevState, fThreadId); 
  G4double sinTheta = fSinTheta;
  G4double cosTheta = sqrt((1-sinTheta)*(1+sinTheta));
  
  GPThreeVector direction = GPThreeVector_create(sinTheta*cos(phi), 
                                                 sinTheta*sin(phi),
                                                 cosTheta);
  
  GPThreeVector refDirection = 
    GPThreeVector_unit(GPThreeVector_create(track->px,track->py,track->pz));
  
  //rotate along the refrefDirection;
  GPThreeVector_rotateUz(&direction,refDirection);
  
  //set secondary 
  G4double energy = fEnergy;

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
