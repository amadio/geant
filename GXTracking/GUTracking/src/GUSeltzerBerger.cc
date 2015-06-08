#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "GUSeltzerBerger.h"
#include <iostream>


#include "backend/Backend.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_HOST
GUSeltzerBerger::GUSeltzerBerger(Random_t* states, int threadId) 
  :
  fRandomState(states), fThreadId(threadId),
  fMinX(1.e-8), fMaxX(1.e+3), // fDeltaX(0.1), 
  fMaxZelement(maximumZ),  
  fNrow(100), fNcol(100) 
{
  //replace hard coded numbers by default constants

  fAliasSampler = new GUAliasSampler(fRandomState,fThreadId,
                                     fMaxZelement,
                                     fMinX,fMaxX,
                                     fNrow,fNcol);

  fDataSB =
    (Physics2DVector*) malloc(maximumZ*sizeof(Physics2DVector));

  char sbDataFile[256];

  for(int iZ = 0 ; iZ < maximumZ ; iZ++) {  
    sprintf(sbDataFile,"data/brem_SB/br%d",iZ+1);
    std::ifstream fin(sbDataFile);
    bool check = RetrieveSeltzerBergerData(fin, &fDataSB[iZ]);
    if(!check) {
      printf("Failed To open SeltzerBerger Data file for Z= %d\n",iZ+1);
    }
  }

  for( int z= 1 ; z < fMaxZelement; ++z)
  {
    //eventually arguments of BuildTable should be replaced by members of *this
    //  and dropped from the function signature. Same for BuildPdfTable
    BuildOneTable(z, fMinX, fMaxX, fNrow, fNcol);
  }

}

VECPHYS_CUDA_HEADER_BOTH 
GUSeltzerBerger::GUSeltzerBerger(Random_t* states, int threadId,
                                 GUAliasSampler* sampler,
				 Physics2DVector* sbData) 
  :
  fRandomState(states), fThreadId(threadId),
  fMinX(1.e-8), fMaxX(1.e+3), // fDeltaX(0.1), 
  fNrow(100), fNcol(100) 
{
  //replace hard coded numbers by default constants
  fAliasSampler = sampler;
  fDataSB = sbData;
}

//need another Ctor with setable parameters

VECPHYS_CUDA_HEADER_BOTH 
GUSeltzerBerger::~GUSeltzerBerger() 
{
  if(fAliasSampler) delete fAliasSampler;
}


VECPHYS_CUDA_HEADER_HOST bool
GUSeltzerBerger::RetrieveSeltzerBergerData(std::ifstream& in, 
                                           Physics2DVector *vec2D)
{
  // binning
  int k;
  int dummyX; // 32 fixed up to Z = 92
  int dummyY; // 57 fixed up to Z = 92
  in >> k >> dummyX >> dummyY;
  if (in.fail())  { return false; }

  // contents
  double valx, valy, val;
  for(size_t i = 0; i< numberOfXNodes; ++i) {
    in >> valx;
    if (in.fail())  { return false; }
    vec2D->PutX(i,valx);
   }
  for(size_t j = 0; j< numberOfYNodes; ++j) {
    in >> valy;
    if (in.fail())  { return false; }
    vec2D->PutY(j,valy);
   }
  for(size_t j = 0; j< numberOfYNodes; ++j) {
    for(size_t i = 0; i< numberOfXNodes; ++i) {
      in >> val;
      if (in.fail())  { return false; }
      vec2D->PutValue(i, j, val);
     }
  }
  in.close();
  return true;

}

VECPHYS_CUDA_HEADER_HOST void 
GUSeltzerBerger::BuildOneTable( int Z, 
                               const double xmin, 
                               const double xmax,
                               const int nrow,
                               const int ncol)
{
  //for now, the model does not own pdf.  Otherwise, pdf should be the 
  //data member of *this and set its point to the fpdf of fAliasSampler 
  double *pdf = new double [(nrow+1)*ncol];

  BuildLogPdfTable(Z,xmin,xmax,nrow,ncol,pdf); 
  fAliasSampler->BuildAliasTable(Z,nrow,ncol,pdf);

  delete [] pdf;
}

VECPHYS_CUDA_HEADER_HOST void 
GUSeltzerBerger::BuildLogPdfTable(int Z, 
                                  const double xmin, 
                                  const double xmax,
                                  const int nrow,
                                  const int ncol,
                                  double *p
                                  )
{
  // Build the probability density function (SeltzerBerger pdf) in the 
  // energy range [xmin,xmax] with an equal bin size (in the log scale)
  //
  // input  :  Z    (atomic number) 
  //           xmin (miminum energy of electron)
  //           xmax (maxinum energy of electron)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  //build pdf  

  double logxmin = log(xmin);
  double dx = (log(xmax) - logxmin)/nrow;

  for(int i = 0; i <= nrow ; ++i) {
    //for each input energy bin
    double x = exp(logxmin + dx*i);

    double emin = (xmin < x) ? xmin : x;
    double emax = (xmax < x) ? xmax : x;

    //total energy
    double t = x + electron_mass_c2;

      //density correction df: should be input 
    double df = 1.0; //test 
    double dc = df*t*t;

    double ymin = log(emin*emin + dc);
    double ymax = log(emax*emax + dc);

    double dy = (ymax - ymin)/(ncol-1); 
    double yo = ymin + 0.5*dy;

    double logx = log(x);

    double sum = 0.0;

    for(int j = 0; j < ncol ; ++j) {
      //for each output energy bin
      double y = exp(yo + dy*j) - dc;
      double w = (y < 0 ) ? 0 : sqrt(y)/x;

      double xsec = CalculateDiffCrossSection(Z,w,logx);
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

// function implementing the cross section for MollerBhabha

VECPHYS_CUDA_HEADER_BOTH double 
GUSeltzerBerger::CalculateDiffCrossSection( int Zelement, 
                                            double w, 
					    double y) const
{
  // based on Geant4
  // data   : SeltzerBerger parameterization (G4LEDATA data set)
  // input  : Zelement (atomic number)
  //          w        (ratio of photon energy to electron energy)
  //          y        (log of the incident electron energy)
  // output : dsigma  (differential cross section) 

  //cross section based on the Seltzer-Berger Parameterization
  double dcross = fDataSB[Zelement].Value(w,y);
  return dcross;
}

} // end namespace impl
} // end namespace vecphys
