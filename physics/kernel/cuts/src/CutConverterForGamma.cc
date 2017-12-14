

#include "CutConverterForGamma.h"
// from material
#include "Types.h"

#include "Material.h"
#include "MaterialProperties.h"
#include "Element.h"

#include <cmath>
#include <iostream>

namespace geantphysics {

CutConverterForGamma::CutConverterForGamma(int numebins, double mincutenergy, double maxcutenergy)
: CutConverter(0, numebins, mincutenergy, maxcutenergy)
  , fZ(-1.), fS200keV(0.), fTmin(0.), fSmin(0.), fCmin(0.), fTlow(0.), fSlow(0.), fS1keV(0.), fClow(0.), fChigh(0.) {
  if (fMinCutEnergy>=fMaxCutEnergy) {
    std::cerr << "  *** ERROR in CutConverterForGamma::CutConverterForGamma() " << std::endl
              << "       minenergy = "<< mincutenergy/geant::GeV
              << " [GeV] >= maxenergy = "
              << maxcutenergy/geant::GeV << " [GeV]"
              << std::endl;
    exit(-1);
  }
  Initialise();
}


// must be called before using the Convert method if new element has been inserted into the Element table!
void CutConverterForGamma::Initialise() {
  CutConverter::Initialise();
  BuildElossOrAbsXsecTable();
}


CutConverterForGamma::~CutConverterForGamma() {}


void CutConverterForGamma::BuildLengthVector(const Material *mat) {
  const Vector_t<Element*> elemVect       = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int   numElements = elemVect.size();

  double maxAbsLenght = -1.0;
  for (int iener=0; iener<fNumEBins; ++iener) {
    double macXsec = 0.0;
    for (int ielem=0; ielem<numElements; ++ielem) {
      int izet     = std::lrint(elemVect[ielem]->GetZ());
      double *vect = fElossOrAbsXsecTable[izet-1];
      macXsec += theAtomicNumDensityVector[ielem]*vect[iener];
    }
    double absorptionLength = 5.0/macXsec;
    fLengthVector[iener]    = absorptionLength;
    if (maxAbsLenght<absorptionLength) {
      maxAbsLenght   = absorptionLength;
      fMaxLengthIndx = iener;
    }
  }
}


// Compute the photon "absorption" cross section: sum of destructive (approximated) cross sections like
// pair production, Compton scattering and photoelectric effect (taken from Geant4)
double CutConverterForGamma::ComputeELossOrAbsXsecPerAtom(double zet, double ekin) {
  const double t1keV   =   1.0*geant::keV;
  const double t200keV = 200.0*geant::keV;
  const double t100MeV = 100.0*geant::MeV;
  //  compute Z dependent quantities if the cached Z is different than zet
  if (std::abs(zet-fZ)>0.1) {
    fZ = zet;
    double Zsquare    = fZ*fZ;
    double Zlog       = std::log(fZ);
    double Zlogsquare = Zlog*Zlog;
    // set some Z dependent variables
    fS200keV = (0.2651-0.1501*Zlog+0.02283*Zlogsquare)*Zsquare;
    fTmin    = (0.552+218.5/fZ+557.17/Zsquare)*geant::MeV;
    fSmin    = (0.01239+0.005585*Zlog-0.000923*Zlogsquare)*std::exp(1.5*Zlog);
    fCmin    = std::log(fS200keV/fSmin)/(std::log(fTmin/t200keV)*std::log(fTmin/t200keV));
    fTlow    = 0.2*std::exp(-7.355/std::sqrt(fZ))*geant::MeV;
    fSlow    = fS200keV*std::exp(0.042*fZ*std::log(t200keV/fTlow)*std::log(t200keV/fTlow));
    fS1keV   = 300.0*Zsquare;
    fClow    = std::log(fS1keV/fSlow)/std::log(fTlow/t1keV);
    fChigh   = (7.55e-5-0.0542e-5*fZ)*Zsquare*fZ/std::log(t100MeV/fTmin);
  }
  // calculate the absorption cross section (using an approximate empirical formula)
  double xs = 0.0;
  if (ekin<fTlow) {
    if (ekin<t1keV)
      xs = fSlow*std::exp(fClow*std::log(fTlow/t1keV));
    else
      xs = fSlow*std::exp(fClow*std::log(fTlow/ekin));
  } else if (ekin<t200keV) {
    xs = fS200keV*std::exp(0.042*fZ*std::log(t200keV/ekin)*std::log(t200keV/ekin));
  } else if (ekin<fTmin) {
    double dum = std::log(fTmin/ekin);
    xs = fSmin*std::exp(fCmin*dum*dum);
  } else {
    xs = fSmin+fChigh*std::log(ekin/fTmin);
  }
  return xs*geant::barn;
}

} // namespace geantphysics
