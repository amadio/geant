#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "BremSeltzerBerger.h"
#include <iostream>

#include "base/Global.h"
#include "GUG4TypeDef.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

//const double


VECCORE_CUDA_HOST
BremSeltzerBerger::BremSeltzerBerger(Random_t* states, int tid)
  : EmModelBase<BremSeltzerBerger>(states,tid)
{
  fAtomicDependentModel = true;
  SetLowEnergyLimit(0.1*keV);

  Initialization();

  //Geant4
  totalEnergy = currentZ = densityFactor = densityCorr = lpmEnergy = xiLPM =
  phiLPM = gLPM = z13 = z23 = lnZ = Fel = Finel = fCoulomb = fMax = 0;

}

VECCORE_CUDA_HOST_DEVICE
BremSeltzerBerger::BremSeltzerBerger(Random_t* states, int tid,
                                     GUAliasSampler* sampler,
                                     Physics2DVector* sbData)
  : EmModelBase<BremSeltzerBerger>(states,tid,sampler)
{
  fAtomicDependentModel = true;
  SetLowEnergyLimit(0.1*keV);

  fDataSB = sbData;
}

//need another Ctor with setable parameters

VECCORE_CUDA_HOST_DEVICE
BremSeltzerBerger::~BremSeltzerBerger()
{
  free(fDataSB);
}

VECCORE_CUDA_HOST void
BremSeltzerBerger::Initialization()
{
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

  if(fSampleType == kAlias && fThreadId != -2) {
    //fThreadId = -2 is used not to build the alias table for testing the Geant4
    //composition and rejection method. This convention should be temporary
    //and is only used for the purpose of validation or other testings.
    fAliasSampler = new GUAliasSampler(fRandomState, fThreadId,
				       1.e-4, 1.e+6, 100, 100);
    BuildAliasTable(fAtomicDependentModel);
  }
}

VECCORE_CUDA_HOST bool
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

VECCORE_CUDA_HOST void
BremSeltzerBerger::BuildCrossSectionTablePerAtom(int Z)
{
  ; //dummy for now
}

VECCORE_CUDA_HOST void
BremSeltzerBerger::BuildPdfTable(int Z, double *p)
{
  // Build the probability density function (SeltzerBerger pdf) in the
  // input energy randge [minX,maxX] with an equal logarithmic bin size
  //
  // input  :  Z    (atomic number)
  // output :  p[nrow][ncol] (probability distribution)
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  const int nrow = fAliasSampler->GetNumEntries();
  const int ncol = fAliasSampler->GetSamplesPerEntry();
  const double minX = fAliasSampler->GetIncomingMin();
  const double maxX = fAliasSampler->GetIncomingMax();

  double logxmin = math::Log(minX);
  double dx = (math::Log(maxX) - logxmin)/nrow;

  for(int i = 0; i <= nrow ; ++i) {
    //for each input energy bin
    double x = math::Exp(logxmin + dx*i);

    double emin = (minX < x) ? minX : x;
    double emax = (maxX < x) ? maxX : x;

    //total energy
    double t = x + electron_mass_c2;

      //density correction df: should be input
    double df = 1.0; //test
    double dc = df*t*t;

    double ymin = math::Log(emin*emin + dc);
    double ymax = math::Log(emax*emax + dc);

    double dy = (ymax - ymin)/(ncol-1);
    double yo = ymin + 0.5*dy;

    double logx = math::Log(x);

    double sum = 0.0;

    for(int j = 0; j < ncol ; ++j) {
      //for each output energy bin
      double y = math::Exp(yo + dy*j) - dc;
      double w = (y < 0 ) ? 0 : math::Sqrt(y)/x;
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

VECCORE_CUDA_HOST_DEVICE double
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

VECCORE_CUDA_HOST_DEVICE void
BremSeltzerBerger::SampleByCompositionRejection(int     Z,
                                                double  kineticEnergy,
                                                double& gammaEnergy,
                                                double& sinTheta)
{
  // G4SeltzerBergerModel::SampleSecondaries
  //  G4double cut  = Min(cutEnergy, kineticEnergy);
  //  G4double emax = Min(maxEnergy, kineticEnergy);
  //@@@syj cutEnergy should be get from the material table (cut table).
  //other hard coded numbers are also temporary and will be replaced properly
  G4double cut  = math::Min(1.0*keV, kineticEnergy);
  G4double emax = math::Min(1.0*TeV, kineticEnergy);
  if(cut >= emax) { return; }


  //  SetupForMaterial(particle, couple->GetMaterial(), kineticEnergy);
  //temporary
  double densityFactor =1.0;

  double totalEnergy = kineticEnergy + electron_mass_c2;
  double densityCorr = densityFactor*totalEnergy*totalEnergy;
  //G4double totMomentum = math::Sqrt(kineticEnergy*(totalEnergy + electron_mass_c2));
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
    G4double ylim = math::Min(fDataSB[Z].Value(0.97, 4*math::Log(10.)),
                    1.1*fDataSB[Z].Value(0.97,y));
    if(ylim > vmax) { vmax = ylim; }
  }
  if(x0 < 0.05) { vmax *= 1.2; }

  do {
    G4double auxrand = UniformRandom<backend::Scalar>(fRandomState,fThreadId);
    G4double x = G4Exp(xmin + auxrand*(xmax - xmin)) - densityCorr;
    if(x < 0.0) { x = 0.0; }
    gammaEnergy = math::Sqrt(x);
    G4double x1 = gammaEnergy/kineticEnergy;
    v = fDataSB[Z].Value(x1, y);

    // correction for positrons

  } while (v < vmax*UniformRandom<backend::Scalar>(fRandomState,fThreadId));

  //
  // angles of the emitted gamma. ( Z - axis along the parent particle)
  // use general interface
  //
  sinTheta = SampleSinTheta<backend::Scalar>(gammaEnergy);
}

VECCORE_CUDA_HOST_DEVICE double
BremSeltzerBerger::GetG4CrossSection(double  kineticEnergy,
                                     int Z)
{
  //temporary
  G4double cutEnergy = 1.0*keV;
  G4double maxEnergy = 1.0*TeV;

  //G4eBremsstrahlungRelModel::ComputeCrossSectionPerAtom
  if(kineticEnergy < fLowEnergyLimit) { return 0.0; }

  G4double cut  = math::Min(cutEnergy, kineticEnergy);
  G4double tmax = math::Min(maxEnergy, kineticEnergy);

  if(cut >= tmax) { return 0.0; }

  SetCurrentElement(Z);

  G4double cross = ComputeXSectionPerAtom(cut,kineticEnergy);

  // allow partial integration
  if(tmax < kineticEnergy) { cross -= ComputeXSectionPerAtom(tmax, kineticEnergy); }

  //constant
  G4double bremFactor =
    fine_structure_const*classic_electr_radius*classic_electr_radius*(16./3.);

  cross *= Z*Z*bremFactor;

  return cross;
}

VECCORE_CUDA_HOST_DEVICE
void BremSeltzerBerger::SetCurrentElement(G4double Z)
{
  if(Z != currentZ) {
    currentZ = Z;

    G4int iz = G4int(Z);

    G4double x  = G4double(iz);
    //    z13 = nist->GetZ13(iz);
    //        = g4pow->Z13(iz);
    //        = pz13[Z];
    //        = pow(x,1.0/3.0);

    z13 = pow(x,1.0/3.0);

    z23 = z13*z13;
    //    lnZ = nist->GetLOGZ(iz);
    //        = g4pow->logZ(iz);
    //        = lz[Z];
    //        = math::Log(x);
    lnZ = math::Log(x);

    constexpr G4double Fel_light[5] = {0., 5.31  , 4.79  , 4.74 ,  4.71};
    constexpr G4double Finel_light[5] = {0., 6.144 , 5.621 , 5.805 , 5.924};

    G4double Fel = 0.0;
    G4double Finel = 0.0;

    if (iz <= 4) {
      Fel = Fel_light[iz];
      Finel = Finel_light[iz] ;
    }
    else {
      //    const G4double facFel = math::Log(184.15);
      //    const G4double facFinel = math::Log(1194.);
      //    Fel = facFel - lnZ/3. ;
      //    Finel = facFinel - 2.*lnZ/3. ;
      Fel = math::Log(184.15) - lnZ/3. ;
      Finel = math::Log(1194.) - 2.*lnZ/3. ;
    }

    //    fCoulomb = GetCurrentElement()->GetfCoulomb();
    fCoulomb = ComputeCoulombFactor(currentZ);
    fMax = Fel-fCoulomb + Finel/currentZ  +  (1.+1./currentZ)/12.;
  }
}

VECCORE_CUDA_HOST_DEVICE
G4double BremSeltzerBerger::ComputeXSectionPerAtom(G4double cut,
                                                   G4double kineticEnergy)
{
  G4double cross = 0.0;

  // number of intervals and integration step
  G4double totalEnergy = kineticEnergy + electron_mass_c2;
  G4double vcut = G4Log(cut/totalEnergy);
  G4double vmax = G4Log(kineticEnergy/totalEnergy);
  G4int n = (G4int)(0.45*(vmax - vcut)) + 4;
  //  n=1; //  integration test
  G4double delta = (vmax - vcut)/G4double(n);

  G4double e0 = vcut;
  G4double xs;

  // densityFactor = mat->GetElectronDensity()*fMigdalConstant;
  // energyThresholdLPM=math::Sqrt(densityFactor)*lpmEnergy;
  G4double densityFactor = 1.0;
  G4double energyThresholdLPM=1.0;
  G4double densityCorr = densityFactor*totalEnergy*totalEnergy;

  constexpr double xgi[8]={ 0.0199, 0.1017, 0.2372, 0.4083,
                            0.5917, 0.7628, 0.8983, 0.9801 };
  constexpr double wgi[8]={ 0.0506, 0.1112, 0.1569, 0.1813,
                            0.1813, 0.1569, 0.1112, 0.0506 };

  // integration
  for(G4int l=0; l<n; l++) {

    for(G4int i=0; i<8; i++) {

      G4double eg = G4Exp(e0 + xgi[i]*delta)*totalEnergy;

      if(totalEnergy > energyThresholdLPM) {
	xs = ComputeRelDXSectionPerAtom(eg);
      } else {
	//ToDo:syjun - select G4eBremsstrahlungRelModel or G4SeltzerBergerModel
	xs = ComputeDXSectionPerAtom(eg);
      }
      cross += wgi[i]*xs/(1.0 + densityCorr/(eg*eg));
    }
    e0 += delta;
  }

  cross *= delta;

  return cross;
}

VECCORE_CUDA_HOST_DEVICE
G4double BremSeltzerBerger::ComputeRelDXSectionPerAtom(G4double gammaEnergy)
// Ultra relativistic model
//   only valid for very high energies, but includes LPM suppression
//    * complete screening
{
  if(gammaEnergy < 0.0) { return 0.0; }

  G4double y = gammaEnergy/totalEnergy;
  G4double y2 = y*y*.25;
  G4double yone2 = (1.-y+2.*y2);

  // ** form factors complete screening case **

  // ** calc LPM functions: include ter-mikaelian merging with density effect **
  //  G4double xiLPM, gLPM, phiLPM;  // to be made member variables !!!

  CalcLPMFunctions(gammaEnergy);

  G4double mainLPM =
    xiLPM*(y2 * gLPM + yone2*phiLPM) * ( (Fel-fCoulomb) + Finel/currentZ );
  G4double secondTerm = (1.-y)/12.*(1.+1./currentZ);

  G4double cross = mainLPM+secondTerm;
  return cross;
}

VECCORE_CUDA_HOST_DEVICE
G4double BremSeltzerBerger::ComputeDXSectionPerAtom(G4double gammaEnergy)
// Relativistic model
//  only valid for high energies (and if LPM suppression does not play a role)
//  * screening according to thomas-fermi-Model (only valid for Z>5)
//  * no LPM effect
{
  if(gammaEnergy < 0.0) { return 0.0; }

  G4double y = gammaEnergy/totalEnergy;

  G4double main=0.,secondTerm=0.;

  //@@@ use_completescreening(false) in G4eBremsstrahlungRelModel constructor
  //if (use_completescreening|| currentZ<5) {
  if (currentZ<5) {
    // ** form factors complete screening case **
    main   = (3./4.*y*y - y + 1.) * ( (Fel-fCoulomb) + Finel/currentZ );
    secondTerm = (1.-y)/12.*(1.+1./currentZ);
  }
  else {
    // ** intermediate screening using Thomas-Fermi FF from Tsai
    // only valid for Z>=5**
    G4double dd=100.*electron_mass_c2*y/(totalEnergy-gammaEnergy);
    G4double gg=dd/z13;
    G4double eps=dd/z23;
    //    G4double phi1=Phi1(gg,currentZ);
    //    G4double phi1m2=Phi1M2(gg,currentZ);
    //    G4double psi1=Psi1(eps,currentZ);
    //    G4double psi1m2=Psi1M2(eps,currentZ);
    G4double phi1=Phi1(gg);
    G4double phi1m2=Phi1M2(gg);
    G4double psi1=Psi1(eps);
    G4double psi1m2=Psi1M2(eps);
    main   = (3./4.*y*y - y + 1.) * ( (0.25*phi1-1./3.*lnZ-fCoulomb) +
                                      (0.25*psi1-2./3.*lnZ)/currentZ );
    secondTerm = (1.-y)/8.*(phi1m2+psi1m2/currentZ);
  }
  G4double cross = main+secondTerm;
  return cross;
}

VECCORE_CUDA_HOST_DEVICE
void  BremSeltzerBerger::CalcLPMFunctions(G4double k)
{
  // *** calculate lpm variable s & sprime ***
  // Klein eqs. (78) & (79)

  G4double sprime = math::Sqrt(0.125*k*lpmEnergy/(totalEnergy*(totalEnergy-k)));

  //  G4double s1 = preS1*z23;
  G4double s1 = (1./(184.15*184.15))*z23;
  //  G4double logS1 = 2./3.*lnZ-2.*facFel;
  G4double logS1 = 2./3.*lnZ-2.*math::Log(184.15);
  //  G4double logTS1 = logTwo+logS1;
  G4double logTS1 = math::Log(2.)+logS1;

  xiLPM = 2.;

  if (sprime>1)
    xiLPM = 1.;
  else if (sprime>math::Sqrt(2.)*s1) {
    G4double h  = math::Log(sprime)/logTS1;
    //    xiLPM = 1+h-0.08*(1-h)*(1-sqr(1-h))/logTS1;
    xiLPM = 1+h-0.08*(1-h)*(1-(1-h)*(1-h))/logTS1;
  }

  G4double s0 = sprime/math::Sqrt(xiLPM);

  // *** merging with density effect***  should be only necessary in region
  // "close to" kp, e.g. k<100*kp using Ter-Mikaelian eq. (20.9)
  G4double k2 = k*k;
  s0 *= (1 + (densityCorr/k2) );

  // recalculate Xi using modified s above
  // Klein eq. (75)
  xiLPM = 1.;
  if (s0<=s1) xiLPM = 2.;
  else if ( (s1<s0) && (s0<=1) ) xiLPM = 1. + math::Log(s0)/logS1;

  // *** calculate supression functions phi and G ***
  // Klein eqs. (77)
  G4double s2=s0*s0;
  G4double s3=s0*s2;
  G4double s4=s2*s2;

  if (s0<0.1) {
    // high suppression limit
    phiLPM = 6.*s0 - 18.84955592153876*s2 + 39.47841760435743*s3
      - 57.69873135166053*s4;
    gLPM = 37.69911184307752*s2 - 236.8705056261446*s3 + 807.7822389*s4;
  }
  else if (s0<1.9516) {
    // intermediate suppression
    // using eq.77 approxim. valid s<2.
    phiLPM = 1.-math::Exp(-6.*s0*(1.+(3.-pi)*s0)
                +s3/(0.623+0.795*s0+0.658*s2));
    if (s0<0.415827397755) {
      // using eq.77 approxim. valid 0.07<s<2
      G4double psiLPM = 1-math::Exp(-4*s0-8*s2/(1+3.936*s0+4.97*s2-0.05*s3+7.50*s4));
      gLPM = 3*psiLPM-2*phiLPM;
    }
    else {
      // using alternative parametrisiation
      G4double pre =
        -0.16072300849123999 + s0*3.7550300067531581 + s2*-1.7981383069010097
        + s3*0.67282686077812381 + s4*-0.1207722909879257;
      gLPM = math::Tanh(pre);
    }
  }
  else {
    // low suppression limit valid s>2.
    phiLPM = 1. - 0.0119048/s4;
    gLPM = 1. - 0.0230655/s4;
  }

  // *** make sure suppression is smaller than 1 ***
  // *** caused by Migdal approximation in xi    ***
  if (xiLPM*phiLPM>1. || s0>0.57)  { xiLPM=1./phiLPM; }
}

VECCORE_CUDA_HOST_DEVICE
G4double BremSeltzerBerger::Phi1(G4double gg)
{
  // Thomas-Fermi FF from Tsai, eq.(3.38) for Z>=5
  //  return 20.863 - 2.*math::Log(1. + sqr(0.55846*gg) )
  return 20.863 - 2.*math::Log(1. + (0.55846*gg)*(0.55846*gg) )
    - 4.*( 1. - 0.6*math::Exp(-0.9*gg) - 0.4*math::Exp(-1.5*gg) );
}

VECCORE_CUDA_HOST_DEVICE
G4double BremSeltzerBerger::Phi1M2(G4double gg)
{
  // Thomas-Fermi FF from Tsai, eq. (3.39) for Z>=5
  // return Phi1(gg,Z) -
  return 2./(3.*(1. + 6.5*gg +6.*gg*gg) );
}

VECCORE_CUDA_HOST_DEVICE
G4double BremSeltzerBerger::Psi1(G4double eps)
{
  // Thomas-Fermi FF from Tsai, eq.(3.40) for Z>=5
  //  return 28.340 - 2.*math::Log(1. + sqr(3.621*eps) )
  return 28.340 - 2.*math::Log(1. + (3.621*eps)*(3.621*eps) )
    - 4.*( 1. - 0.7*math::Exp(-8*eps) - 0.3*math::Exp(-29.*eps) );
}

VECCORE_CUDA_HOST_DEVICE
G4double BremSeltzerBerger::Psi1M2(G4double eps)
{
  // Thomas-Fermi FF from Tsai, eq. (3.41) for Z>=5
  return  2./(3.*(1. + 40.*eps +400.*eps*eps) );
}



} // end namespace impl
} // end namespace vecphys
