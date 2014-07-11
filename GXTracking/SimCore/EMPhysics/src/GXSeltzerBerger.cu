#include "GXSeltzerBerger.h"
#include "GPThreeVector.h"

FQUALIFIER
GXSeltzerBerger::GXSeltzerBerger(curandState* devState,
				 G4int threadId,
				 GPPhysics2DVector* data,
				 G4double minEnergy,
				 G4double maxEnergy,
				 G4int dimRow,
				 G4int dimCol,
				 G4double *pdfX,
				 G4double *pdfY)
  :  GXVDiscreteSampling(devState,threadId,
			 dimRow,dimCol,pdfX,pdfY)
{
  fThreadId = threadId;
  fDevState = devState;

  //cross section data for G4SeltzerBergerModel
  fMass = electron_mass_c2; 
  fDataSB = data;
  fDensityCorr = 0.;

  //Sampling tables
  fNrow = dimRow;
  fNcol = dimCol;

  fMinX = minEnergy;
  fMaxX = maxEnergy;
  fDeltaX = (log(fMaxX) - log(fMinX))/fNrow;

  fMinY = 0;
  fMaxY = 0;
  fDeltaY = 0;
}

FQUALIFIER
GXSeltzerBerger::GXSeltzerBerger(curandState* devState,
				 G4int threadId,
				 GPPhysics2DVector* data,
				 G4double minEnergy,
				 G4double maxEnergy,
				 G4int dimRow,
				 G4int dimCol,
				 G4int    *pdfA,
				 G4double *pdfX,
				 G4double *pdfY)
  : GXVDiscreteSampling(devState,threadId,
			dimRow,dimCol,pdfA,pdfX,pdfY)
{
  fThreadId = threadId;
  fDevState = devState;

  //cross section data for G4SeltzerBergerModel
  fMass = electron_mass_c2; 
  fDataSB = data;
  fDensityCorr = 0.;

  //Sampling tables
  fNrow = dimRow;
  fNcol = dimCol;

  fMinX = minEnergy;
  fMaxX = maxEnergy;
  fDeltaX = (log(fMaxX) - log(fMinX))/fNrow;

  fMinY = 0;
  fMaxY = 0;
  fDeltaY = 0;
}

FQUALIFIER GXSeltzerBerger::~GXSeltzerBerger()
{
  ;
}

FQUALIFIER void
GXSeltzerBerger::GetSampleParameters(G4int Z,
				     G4double kineticEnergy,
				     G4int &irow, G4int &icol,
				     G4double &t) 
{
  irow = G4int((log(kineticEnergy) - log(fMinX))/fDeltaX);
  G4double r1 = (fNcol-1)*GPUniformRand(fDevState, fThreadId);
  icol = int(r1);
  t = r1 - 1.0*icol;

  //Z dependency goes here - dummy for now
  G4double densityFactor = (1.0*Z)/Z; 
  
  G4double emin = fmin(fMinX, kineticEnergy);
  G4double emax = fmin(fMaxX, kineticEnergy);

  G4double totalEnergy = kineticEnergy + fMass;
  fDensityCorr = densityFactor*totalEnergy*totalEnergy;
  fMinY = log(emin*emin + fDensityCorr);
  fMaxY = log(emax*emax + fDensityCorr);
  fDeltaY = fMaxY - fMinY;
}

FQUALIFIER void
GXSeltzerBerger::SampleByInversePDF(int Z, double kineticEnergy)
{
  //sample based on the inverse cumulative pdf
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(Z,kineticEnergy,irow,icol,t);

  //sample based on the inverse cumulataive pdf
  G4double yhat = InversePDF(irow,icol,fDeltaY,t);

  fGammaEnergy =  sqrt(fmax(exp(fMinY +yhat) -fDensityCorr,0.0));
  fGammaSinTheta = SampleDirection(kineticEnergy + fMass - fGammaEnergy);
}


FQUALIFIER void
GXSeltzerBerger::SampleByInversePDFLinearInterpolation(int Z,
						       double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(Z,kineticEnergy,irow,icol,t);

  //sample based on the inverse cumulataive pdf
  G4double yhat = InversePDFLinearInterpolation(irow,icol,fDeltaY,t);

  fGammaEnergy = sqrt(fmax(exp(fMinY+ yhat)-fDensityCorr,0.0));
  fGammaSinTheta = SampleDirection(kineticEnergy + fMass - fGammaEnergy);
}

FQUALIFIER void
GXSeltzerBerger::SampleByAlias(G4int Z, G4double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(Z,kineticEnergy,irow,icol,t);

  //sample based on the alias method
  G4double yhat = Alias(irow,icol,fDeltaY,t);

  fGammaEnergy = sqrt(fmax(exp(fMinY + yhat)-fDensityCorr,0.0));
  fGammaSinTheta = SampleDirection(kineticEnergy + fMass - fGammaEnergy);
}

FQUALIFIER void
GXSeltzerBerger::SampleByAliasLinearInterpolation(G4int Z, 
						  G4double kineticEnergy)
{
  G4int irow = -1;
  G4int icol = -1;
  G4double t = 0.;

  GetSampleParameters(Z,kineticEnergy,irow,icol,t);

  //sample based on the alias method
  G4double yhat = AliasLinearInterpolation(irow,icol,fDeltaY,t);

  fGammaEnergy = sqrt(fmax(exp(fMinY + yhat)-fDensityCorr,0.0));
  fGammaSinTheta = SampleDirection(kineticEnergy + fMass - fGammaEnergy);
}

FQUALIFIER void
GXSeltzerBerger::SampleByCompositionRejection(int Z, 
					      double energy,
					      G4int& ntrial)
{
  // G4SeltzerBergerModel::SampleSecondaries

  G4double kineticEnergy = energy ;
  G4double emin = fmin(fMinX, kineticEnergy);
  G4double emax = fmin(fMaxX, kineticEnergy);
  if(emin >= emax) { return; }

  int currentZ = Z;
  const G4double densityFactor =1.0;
  
  G4double totalEnergy = kineticEnergy + fMass;
  G4double densityCorr = densityFactor*totalEnergy*totalEnergy;
  G4double totMomentum = sqrt(kineticEnergy*(totalEnergy + electron_mass_c2));

  G4double xmin = log(emin*emin + densityCorr);
  G4double xmax = log(emax*emax + densityCorr);
  G4double y = log(kineticEnergy/MeV);

  G4double gammaEnergy, v; 

  // majoranta
  G4double x0 = emin/kineticEnergy;
  G4double vmax = fDataSB[Z].Value(x0, y)*1.02;

  const G4double epeaklimit= 300*MeV; 
  const G4double elowlimit = 10*keV; 

  // majoranta corrected for e-
  bool isElectron = true;
  if(isElectron && x0 < 0.97 && 
     ((kineticEnergy > epeaklimit) || (kineticEnergy < elowlimit))) {
    G4double ylim = fmin(fDataSB[Z].Value(0.97, 4*log(10.)),
			 1.1*fDataSB[Z].Value(0.97, y)); 
    if(ylim > vmax) { vmax = ylim; }
  }
  if(x0 < 0.05) { vmax *= 1.2; }
  
  do {
    G4double auxrand = GPUniformRand(fDevState, fThreadId);
    G4double x = exp(xmin + auxrand*(xmax - xmin))-densityCorr;
    if(x < 0.0) { x = 0.0; }
    gammaEnergy = sqrt(x);
    G4double x1 = gammaEnergy/kineticEnergy;
    v = fDataSB[Z].Value(x1, y);

    // correction for positrons        
    
    if(!isElectron) {
      G4double e1 = kineticEnergy - emin;
      G4double invbeta1 = (e1 + fMass)/sqrt(e1*(e1 + 2*fMass));
      G4double e2 = kineticEnergy - gammaEnergy;
      G4double invbeta2 = (e2 + fMass)/sqrt(e2*(e2 + 2*fMass));
      G4double xxx = twopi*fine_structure_const*currentZ*(invbeta1 - invbeta2); 
      
      if(xxx < -12. ) { v = 0.0; }
      else { v *= exp(xxx); }
    }
    ++ntrial;
  } while (v < vmax*GPUniformRand(fDevState, fThreadId));

  fGammaEnergy = gammaEnergy; 
  fGammaSinTheta = SampleDirection(totalEnergy-gammaEnergy);
}

FQUALIFIER G4double 
GXSeltzerBerger::SampleDirection(G4double eTkin)
{
  //GetAngularDistribution()->SampleDirection 
  //see the constructor of G4eBremsstrahlungRelModel where the angle model
  //is set: SetAngularDistribution(new G4DipBustGenerator());
  //
  //G4ThreeVector& 
  //G4DipBustGenerator::SampleDirection(const G4DynamicParticle* dp,
  //                                    G4double, G4int, const G4Material*)
  // change arguments with
  // G4double eTkin = dp->GetKineticEnergy();
  // refDirection = dp->GetMomentumDirection()
  
  G4double c, cosTheta, delta, cofA, signc = 1., a, power = 1./3.;

  c = 4. - 8.*GPUniformRand(fDevState, fThreadId);
  a = c;
 
  if( c < 0. )
  {
    signc = -1.;
    a     = -c;
  }
  delta  = sqrt(a*a+4.);
  delta += a;
  delta *= 0.5; 

  cofA = -signc*pow(delta, power);

  cosTheta = cofA - 1./cofA;

  G4double tau = eTkin/electron_mass_c2;
  G4double beta = sqrt(tau*(tau + 2.))/(tau + 1.);

  cosTheta = (cosTheta + beta)/(1 + cosTheta*beta);

  G4double sinTheta = sqrt((1 - cosTheta)*(1 + cosTheta));

  return sinTheta;
}

FQUALIFIER
G4double GXSeltzerBerger::GetSecondaryEnergy()
{
  return fGammaEnergy;
}

FQUALIFIER
G4double GXSeltzerBerger::GetSecondarySinTheta()
{
  return fGammaSinTheta;
}

FQUALIFIER
void GXSeltzerBerger::CreateSecondary(GXTrack* secondary,
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
