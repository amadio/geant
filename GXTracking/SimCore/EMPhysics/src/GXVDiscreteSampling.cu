#include "GXVDiscreteSampling.h"
#include "GXSamplingTexture.h"

FQUALIFIER
GXVDiscreteSampling::GXVDiscreteSampling(curandState* devState,
					 G4int threadId,
					 G4int dimRow,
					 G4int dimCol,
					 G4double *pdfX,
					 G4double *pdfY)
{
  //Random state
  fThreadId = threadId;
  fDevState = devState;
  
  //Sampling Tables
  fNrow = dimRow;
  fNcol = dimCol;

  fPDFA = NULL;
  fPDFX = pdfX;
  fPDFY = pdfY;
}

FQUALIFIER
GXVDiscreteSampling::GXVDiscreteSampling(curandState* devState,
                                         G4int threadId,
					 G4int dimRow,
					 G4int dimCol,
                                         G4int    *pdfA,
                                         G4double *pdfX,
                                         G4double *pdfY)  
{
  //Random state
  fThreadId = threadId;
  fDevState = devState;
  
  //Sampling Tables
  fNrow = dimRow;
  fNcol = dimCol;

  fPDFA = pdfA;
  fPDFX = pdfX;
  fPDFY = pdfY;
}

FQUALIFIER GXVDiscreteSampling::~GXVDiscreteSampling()
{
  ;
}


FQUALIFIER G4double
GXVDiscreteSampling::InversePDF(G4int irow, G4int icol, 
				G4double dy, G4double tt)
{
  G4double xd = dy*fPDFY[irow*fNcol+icol];
  G4double xu = dy*fPDFY[irow*fNcol+icol+1];
  G4double x = (1.0-tt)*xd + tt*xu;

  return x;
}

FQUALIFIER G4double
GXVDiscreteSampling::InversePDFTexture(G4int irow, G4int icol, 
				       G4double dy, G4double tt)
{
#ifdef __CUDA_ARCH__
  G4double xd = dy*FetchToDouble(texPDFY,irow*fNcol+icol);
  G4double xu = dy*FetchToDouble(texPDFY,irow*fNcol+icol+1);
#else
  G4double xd = dy*fPDFY[irow*fNcol+icol];
  G4double xu = dy*fPDFY[irow*fNcol+icol+1];
#endif
  G4double x = (1.0-tt)*xd + tt*xu;

  return x;
}

FQUALIFIER G4double
GXVDiscreteSampling::InversePDFLinearInterpolation(G4int irow, G4int icol, 
						   G4double dy, G4double tt)
{
  G4double pd = fPDFX[irow*fNcol+icol];
  G4double pu = fPDFX[irow*fNcol+icol+1];

  G4double xd = dy*fPDFY[irow*fNcol+icol];
  G4double xu = dy*fPDFY[irow*fNcol+icol+1];

  G4double p  = (1-tt)*pd + tt*pu;
  G4double r1 = GPUniformRand(fDevState, fThreadId)*(pd + pu);
  G4double x  = (p > r1) ? tt*xd + (1.0-tt)*xu : (1.0-tt)*xd + tt*xu;

  return x;
}

FQUALIFIER G4double
GXVDiscreteSampling::InversePDFTextureLinearInterpolation(G4int irow, 
							  G4int icol, 
							  G4double dy, 
							  G4double tt)
{
#ifdef __CUDA_ARCH__
  G4double pd = FetchToDouble(texPDFX,irow*fNcol+icol);
  G4double pu = FetchToDouble(texPDFX,irow*fNcol+icol+1);

  G4double xd = dy*FetchToDouble(texPDFY,irow*fNcol+icol);
  G4double xu = dy*FetchToDouble(texPDFY,irow*fNcol+icol+1);
#else
  G4double pd = fPDFX[irow*fNcol+icol];
  G4double pu = fPDFX[irow*fNcol+icol+1];

  G4double xd = dy*fPDFY[irow*fNcol+icol];
  G4double xu = dy*fPDFY[irow*fNcol+icol+1];
#endif
  G4double p  = (1-tt)*pd + tt*pu;
  G4double r1 = GPUniformRand(fDevState, fThreadId)*(pd + pu);
  G4double x  = (p > r1) ? tt*xd + (1.0-tt)*xu : (1.0-tt)*xd + tt*xu;

  return x;
}

FQUALIFIER G4double
GXVDiscreteSampling::Alias(G4int irow, G4int icol, 
			   G4double dy, G4double tt)
{
  double r1 = GPUniformRand(fDevState, fThreadId);
  
  G4double xd;
  G4double xu;
  G4double dx = dy/(fNcol-1);
  
  if(r1 <=  fPDFY[irow*fNcol+icol] ) {
    xd = dx*icol;
    xu = dx*(icol+1);
  }    
  else {
    xd = dx*fPDFA[irow*fNcol+icol];
    xu = dx*(fPDFA[irow*fNcol+icol]+1);
  }
  
  G4double x = (1-tt)*xd + tt*xu;
  
  return x;
}

FQUALIFIER G4double
GXVDiscreteSampling::AliasLinearInterpolation(G4int irow, G4int icol, 
					      G4double dy, G4double tt)
{
  double r1 = GPUniformRand(fDevState, fThreadId);
  
  G4int xa;
  G4double xd;
  G4double xu;
  G4double pd;
  G4double pu;

  G4double dx = dy/(fNcol-1);

  if(r1 <=  fPDFY[irow*fNcol+icol] ) {
    xd = dx*icol;
    xu = dx*(icol+1);
    pd = fPDFX[irow*fNcol+icol];
    pu = fPDFX[irow*fNcol+icol+1];
  }    
  else {
    xa = fPDFA[irow*fNcol+icol];    
    xd = dx*xa;
    xu = dx*(xa+1);
    pd = fPDFX[xa];
    pu = fPDFX[xa+1];
  }

  G4double r2 = GPUniformRand(fDevState, fThreadId);
  G4double x = (r2*(pd+pu) < (1-tt)*pd + tt*pu) ?
    (1-tt)*xd + tt*xu : tt*xd + (1-tt)*xu;
  
  return x;
}
