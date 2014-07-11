#ifndef GXKleinNishina_H
#define GXKleinNishina_H 1

#include "GPTypeDef.h"
#include "GPRandom.h"
#include "GPThreeVector.h"
#include "GXTrack.h"

#include "GXVDiscreteSampling.h"

class GXKleinNishina : public GXVDiscreteSampling
{
public:

  FQUALIFIER GXKleinNishina(curandState* devState,
			    int threadId,
			    G4double minE, 
			    G4double maxE,
			    G4int dimRow,
			    G4int dimCol,
			    G4double *pdfX,
			    G4double *pdfY);

  FQUALIFIER GXKleinNishina(curandState* devState,
			    int threadId,
			    G4double minE, 
			    G4double maxE,
			    G4int dimRow,
			    G4int dimCol,
			    G4int    *pdfA,
			    G4double *pdfX,
			    G4double *pdfY);

  FQUALIFIER ~GXKleinNishina();
  
  FQUALIFIER void GetSampleParameters(G4double x, G4int &irow, G4int &icol,
				      G4double &t);

  FQUALIFIER void SampleByInversePDF(G4double energy);
  FQUALIFIER void SampleByInversePDFLinearInterpolation(G4double energy);
  FQUALIFIER void SampleByAlias(G4double energy);
  FQUALIFIER void SampleByAliasLinearInterpolation(G4double energy);
  
  FQUALIFIER void SampleByCompositionRejection(G4double energy,
					       int& ntrial);
  
  FQUALIFIER G4double SampleDirection(G4double energyIn, G4double energyOut);

  FQUALIFIER G4double GetSecondaryEnergy();
  FQUALIFIER G4double GetSecondarySinTheta();
  FQUALIFIER void CreateSecondary(GXTrack* secondary, GXTrack* track, 
				  G4double charge);

private:
  G4int fThreadId;
  curandState* fDevState;

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

