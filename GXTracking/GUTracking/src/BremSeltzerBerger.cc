#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "BremSeltzerBerger.h"
#include <iostream>

#include "backend/Backend.h"
#include "GUG4TypeDef.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_HOST
BremSeltzerBerger::BremSeltzerBerger(Random_t* states, int tid) 
  : EmModelBase<BremSeltzerBerger>(states,tid)
{
  SetLowEnergyLimit(10.*keV);

  fDataSB =
    (Physics2DVector*) malloc(maximumZ*sizeof(Physics2DVector));

  char sbDataFile[256];

  for(int iZ = 0 ; iZ < maximumZ ; iZ++) {  
    sprintf(sbDataFile,"data/brem_SB/br%d",iZ+1);
    std::ifstream fin(sbDataFile);
    bool check = RetrieveSeltzerBergerData(fin, &fDataSB[iZ]);
    if(!check) {
      printf("Failed To open eltzerBerger Data file for Z= %d\n",iZ+1);
    }
  }

  BuildAliasTable();
}

VECPHYS_CUDA_HEADER_BOTH 
BremSeltzerBerger::BremSeltzerBerger(Random_t* states, int tid,
                                     GUAliasSampler* sampler, 
                                     Physics2DVector* sbData) 
  : EmModelBase<BremSeltzerBerger>(states,tid,sampler)
{
  SetLowEnergyLimit(10.*keV);
  fDataSB = sbData;
}

//need another Ctor with setable parameters

VECPHYS_CUDA_HEADER_BOTH 
BremSeltzerBerger::~BremSeltzerBerger() 
{
  free(fDataSB);
}

VECPHYS_CUDA_HEADER_HOST bool
BremSeltzerBerger::RetrieveSeltzerBergerData(std::ifstream& in, 
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
BremSeltzerBerger::BuildLogPdfTable(int Z, 
                                    const double xmin, 
                                    const double xmax,
                                    const int nrow,
                                    const int ncol,
                                    double *p)
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

// function implementing the differential cross section for SeltzerBerger

VECPHYS_CUDA_HEADER_BOTH double 
BremSeltzerBerger::CalculateDiffCrossSection(int Zelement, 
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

VECPHYS_CUDA_HEADER_BOTH void 
BremSeltzerBerger::SampleByCompositionRejection(int     Z,
                                                double  kineticEnergy,
                                                double& gammaEnergy,
                                                double& sinTheta)
{
  // G4SeltzerBergerModel::SampleSecondaries
  //  G4double cut  = Min(cutEnergy, kineticEnergy);
  //  G4double emax = Min(maxEnergy, kineticEnergy);
  G4double cut  = Min(fMinX, kineticEnergy);
  G4double emax = Min(fMaxX, kineticEnergy);
  if(cut >= emax) { return; }
 

  //  SetupForMaterial(particle, couple->GetMaterial(), kineticEnergy);
  //temporary
  double densityFactor =1.0; 
 
  double totalEnergy = kineticEnergy + electron_mass_c2;
  double densityCorr = densityFactor*totalEnergy*totalEnergy;
  //G4double totMomentum = Sqrt(kineticEnergy*(totalEnergy + electron_mass_c2));
  G4double xmin = G4Log(cut*cut + densityCorr);
  G4double xmax = G4Log(emax*emax  + densityCorr);
  G4double y = G4Log(kineticEnergy/MeV);
 
  G4double v; 
 
  // majoranta
  G4double x0 = cut/kineticEnergy;
  G4double vmax = fDataSB[Z].Value(x0, y)*1.02;

  const double_t epeaklimit= 300*MeV; 
  const double_t elowlimit = 10*keV; 

  // majoranta corrected for e-
  bool isElectron = true;
  if(isElectron && x0 < 0.97 && 
     ((kineticEnergy > epeaklimit) || (kineticEnergy < elowlimit))) {
    G4double ylim = Min(fDataSB[Z].Value(0.97, 4*log(10.)),
                    1.1*fDataSB[Z].Value(0.97,y));
    if(ylim > vmax) { vmax = ylim; }
  }
  if(x0 < 0.05) { vmax *= 1.2; }
 
  do {
    G4double auxrand = UniformRandom<kScalar>(fRandomState,fThreadId);
    G4double x = G4Exp(xmin + auxrand*(xmax - xmin)) - densityCorr;
    if(x < 0.0) { x = 0.0; }
    gammaEnergy = Sqrt(x);
    G4double x1 = gammaEnergy/kineticEnergy;
    v = fDataSB[Z].Value(x1, y);
 
    // correction for positrons        
   
  } while (v < vmax*UniformRandom<kScalar>(fRandomState,fThreadId));
 
  //
  // angles of the emitted gamma. ( Z - axis along the parent particle)
  // use general interface
  //
  sinTheta = SampleSinTheta<kScalar>(gammaEnergy); 
}

VECPHYS_CUDA_HEADER_BOTH double
BremSeltzerBerger::GetG4CrossSection(double  kineticEnergy, 
                                     int Z)
{
  //temporary
  G4double cutEnergy = 1.0*keV;
  G4double maxEnergy = 1.0*TeV;

  //G4eBremsstrahlungRelModel::ComputeCrossSectionPerAtom
  if(kineticEnergy < fLowEnergyLimit) { return 0.0; }

  G4double cut  = Min(cutEnergy, kineticEnergy);
  G4double tmax = Min(maxEnergy, kineticEnergy);

  if(cut >= tmax) { return 0.0; }

  //  SetCurrentElement(Z);

  G4double cross = ComputeXSectionPerAtom(cut);

  // allow partial integration
  if(tmax < kineticEnergy) { cross -= ComputeXSectionPerAtom(tmax); }

  //constant
  G4double bremFactor = 
    fine_structure_const*classic_electr_radius*classic_electr_radius*(16./3.);
  
  cross *= Z*Z*bremFactor;

  return cross;
}

VECPHYS_CUDA_HEADER_BOTH
G4double BremSeltzerBerger::ComputeXSectionPerAtom(G4double cut)
{
  G4double cross = 0.0;

  // number of intervals and integration step 
  /*
  G4double vcut = G4Log(cut/totalEnergy);
  G4double vmax = G4Log(kinEnergy/totalEnergy);
  G4int n = (G4int)(0.45*(vmax - vcut)) + 4;
  //  n=1; //  integration test 
  G4double delta = (vmax - vcut)/G4double(n);

  G4double e0 = vcut;
  G4double xs; 

  // integration
  for(G4int l=0; l<n; l++) {

    for(G4int i=0; i<8; i++) {

      G4double eg = G4Exp(e0 + xgi[i]*delta)*totalEnergy;

      if(totalEnergy > energyThresholdLPM) {
  xs = ComputeRelDXSectionPerAtom(eg);
      } else {
  xs = ComputeDXSectionPerAtom(eg);
      }
      cross += wgi[i]*xs/(1.0 + densityCorr/(eg*eg));
    }
    e0 += delta;
  }

  cross *= delta;
  */
  return cross;
}


} // end namespace impl
} // end namespace vecphys
