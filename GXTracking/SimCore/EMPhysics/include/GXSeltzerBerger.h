#ifndef GXSeltzerBerger_H
#define GXSeltzerBerger_H 1

#include "GPTypeDef.h"
#include "GPRandom.h"
#include "GPPhysics2DVector.h"
#include "GXTrack.h"

#include "GXVDiscreteSampling.h"

class GXSeltzerBerger : public GXVDiscreteSampling
{
public:

  FQUALIFIER GXSeltzerBerger(curandState* devState,
			     int threadId,
			     GPPhysics2DVector* data, 
			     G4double minE, 
			     G4double maxE,
			     G4int dimRow,
			     G4int dimCol,
			     G4double *pdfX,
			     G4double *pdfY);
  
  FQUALIFIER GXSeltzerBerger(curandState* devState,
			     int threadId,
			     GPPhysics2DVector* data, 
			     G4double minE, 
			     G4double maxE,
			     G4int dimRow,
			     G4int dimCol,
			     G4int    *pdfA,
			     G4double *pdfX,
			     G4double *pdfY);
  
  FQUALIFIER ~GXSeltzerBerger();



  FQUALIFIER void GetSampleParameters(G4int Z, G4double x, 
				      G4int &irow, G4int &icol,
                                      G4double &t);

  //material dependency with Z for the future implementation
  FQUALIFIER void SampleByRandom(G4double energy);
  FQUALIFIER void SampleByInversePDF(G4int Z, G4double energy);
  FQUALIFIER void SampleByInversePDFLinearInterpolation(G4int Z,
							G4double energy);
  FQUALIFIER void SampleByInversePDFTexture(G4int Z, G4double energy);
  FQUALIFIER void SampleByInversePDFTextureLinearInterpolation(G4int Z,
							G4double energy);
  FQUALIFIER void SampleByAlias(G4int Z, G4double energy);
  FQUALIFIER void SampleByAliasLinearInterpolation(G4int Z, G4double energy);
  
  FQUALIFIER void SampleByCompositionRejection(G4int Z, G4double energy,
					       G4int& ntrial);
  FQUALIFIER void SampleByAverageTrials(G4int Z,
					int ntrials, 
                                        G4double energy,
                                        int& counter);

  FQUALIFIER G4double SampleDirection(G4double eTkin);
  
  FQUALIFIER G4double GetSecondaryEnergy();
  FQUALIFIER G4double GetSecondarySinTheta();
  FQUALIFIER void CreateSecondary(GXTrack* secondary, GXTrack* track, 
				  G4double charge);
  
 private:
  int fThreadId;
  curandState* fDevState;
  
  //SeltzerBerger data
  G4double fMass;
  GPPhysics2DVector* fDataSB;
  G4double fDensityCorr;

  G4double fMinX;
  G4double fMaxX;
  G4double fDeltaX;

  G4double fMinY;
  G4double fMaxY;
  G4double fDeltaY;

  //Sampling Tables
  G4int fNrow;
  G4int fNcol;

  //secondary energy and angle
  G4double fGammaEnergy;
  G4double fGammaSinTheta;

};

#endif

