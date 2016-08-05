#include "GUAliasSampler.h"
#include "GUAliasTable.h"

#include "ElectronProcess.h"

// from VecGeom
#include "materials/Material.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_ATT_HOST
ElectronProcess::ElectronProcess(Random_t *states, int tid) : EmProcess<ElectronProcess>(states, tid)
{
  fNumberOfProcess = 3;
  fIonisation = 0;
  fBremsstrahlung = 0;
  fMSC = 0;

  fLogEnergyLowerBound = log(fEnergyLowerBound);
  fInverseLogEnergyBin = fNumberOfEnergyBin / (log(fEnergyUpperBound) - fLogEnergyLowerBound);

  Initialization();
}

VECCORE_ATT_HOST_DEVICE
ElectronProcess::ElectronProcess(Random_t *states, int tid, CrossSectionData *data)
    : EmProcess<ElectronProcess>(states, tid, data)
{
  fNumberOfProcess = 3;
  fIonisation = 0;
  fBremsstrahlung = 0;
  fMSC = 0;

  fLogEnergyLowerBound = log(fEnergyLowerBound);
  fInverseLogEnergyBin = fNumberOfEnergyBin / (log(fEnergyUpperBound) - fLogEnergyLowerBound);

  fNumberOfMaterialBin = 3; //(vecgeom::Material::GetMaterials()).size();
}

VECCORE_ATT_HOST
ElectronProcess::~ElectronProcess()
{
  delete fIonisation;
  delete fBremsstrahlung;
  delete fMSC;
}

VECCORE_ATT_HOST void ElectronProcess::Initialization()
{
  fIonisation = new IonisationMoller(0, -1);
  fBremsstrahlung = new BremSeltzerBerger(0, -1);
  fMSC = new UrbanWentzelVI(0, -1);

  fNumberOfMaterialBin = (vecgeom::Material::GetMaterials()).size();

  fCrossSectionData = (CrossSectionData *)malloc(sizeof(CrossSectionData) * fNumberOfEnergyBin * fNumberOfMaterialBin);
  // initialize table
  for (int i = 0; i < fNumberOfMaterialBin; ++i) {
    for (int j = 0; j < fNumberOfEnergyBin; ++j) {
      int ibin = i * fNumberOfEnergyBin + j;
      fCrossSectionData[ibin].fSigma = 0.0;

      for (int k = 0; k < fNumberOfProcess; ++k)
        fCrossSectionData[ibin].fAlias[k] = 0;
      for (int k = 0; k < fNumberOfProcess - 1; ++k)
        fCrossSectionData[ibin].fWeight[k] = 0.0;
    }
  }

  BuildCrossSectionTable();
}

VECCORE_ATT_HOST void ElectronProcess::BuildCrossSectionTable()
{
  // Get the material DB
  std::vector<vecgeom::Material *> &mtable = vecgeom::Material::GetMaterials();

  double cross[3] = {
      0., 0., 0.,
  };
  double energy = 0;

  double logBinInterval = (log(fEnergyUpperBound) - fLogEnergyLowerBound) / fNumberOfEnergyBin;

  for (int i = 0; i < fNumberOfMaterialBin; ++i) {
    for (int j = 0; j < fNumberOfEnergyBin; ++j) {
      energy = exp(fLogEnergyLowerBound + logBinInterval * (j + 0.5));
      int ibin = i * fNumberOfEnergyBin + j;

      cross[0] = fIonisation->G4CrossSectionPerVolume((mtable)[i], energy);
      cross[1] = fBremsstrahlung->G4CrossSectionPerVolume((mtable)[i], energy);
      cross[2] = fMSC->G4CrossSectionPerVolume((mtable)[i], energy);

      // fill cross section information (total and weights)
      double totalCrossSection = cross[0] + cross[1] + cross[2];
      fCrossSectionData[ibin].fSigma = totalCrossSection;

      // fill cross section information (total and weights)
      if (totalCrossSection > 0.0) {
        fCrossSectionData[ibin].fWeight[0] = cross[0] / totalCrossSection;
        fCrossSectionData[ibin].fWeight[1] = cross[1] / totalCrossSection;
      }
    }
  }
}

} // end namespace impl
} // end namespace vecphys
