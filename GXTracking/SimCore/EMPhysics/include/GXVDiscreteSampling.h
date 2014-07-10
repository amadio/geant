#ifndef GXVDiscreteSampling_H
#define GXVDiscreteSampling_H 1

#include "GPTypeDef.h"
#include "GPRandom.h"

class GXVDiscreteSampling
{
public:

  FQUALIFIER GXVDiscreteSampling(curandState* devState,
				 G4int threadId,
				 G4int dimRow,
				 G4int dimCol,
				 G4double *pdfX,
				 G4double *pdfY);

  FQUALIFIER GXVDiscreteSampling(curandState* devState,
                                 G4int threadId,
				 G4int dimRow,
				 G4int dimCol,
                                 G4int    *pdfA,
                                 G4double *pdfX,
                                 G4double *pdfY);
  
  FQUALIFIER ~GXVDiscreteSampling();

  FQUALIFIER G4double InversePDF(G4int irow, G4int icol, 
				 G4double dy, G4double tt);
  FQUALIFIER G4double InversePDFLinearInterpolation(G4int irow, G4int icol, 
						    G4double dy, G4double tt);
  FQUALIFIER G4double Alias(G4int irow, G4int icol, 
			    G4double dy, G4double tt);
  FQUALIFIER G4double AliasLinearInterpolation(G4int irow, G4int icol, 
					       G4double dy, G4double tt);

private:
  //Random state
  G4int fThreadId;
  curandState* fDevState;

  //Sampling Tables
  G4int     fNrow;
  G4int     fNcol;

  G4int    *fPDFA; 
  G4double *fPDFX; 
  G4double *fPDFY; 
};

#endif

