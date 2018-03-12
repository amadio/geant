#include <Geant/PhysicalConstants.h>
#include <TString.h>
#include "Geant/VecKleinNishinaComptonModel.h"
#include "Geant/AliasTable.h"

namespace geantphysics {

void VecKleinNishinaComptonModel::Initialize()
{
  KleinNishinaComptonModel::Initialize();
  for (int i = 0; i < fSamplingTables.size(); ++i) {
    LinAlias *alias = fSamplingTables[i];
    int size        = (int)alias->fAliasIndx.size();

    LinAliasCached aliasCached(size);
    for (int j = 0; j < size - 1; ++j) { // TODO: What if sampled index is really numdata - 1?
      aliasCached.fAliasIndx[j]               = alias->fAliasIndx[j];
      aliasCached.fAliasW[j]                  = alias->fAliasW[j];
      aliasCached.fLinAliasData[j].X          = alias->fXdata[j];
      aliasCached.fLinAliasData[j].Ydata      = alias->fYdata[j];
      aliasCached.fLinAliasData[j].Xdelta     = alias->fXdata[j + 1] - alias->fXdata[j];
      aliasCached.fLinAliasData[j].YdataDelta = (alias->fYdata[j + 1] - alias->fYdata[j]) / alias->fYdata[j];
      aliasCached.fLinAliasData[j].XdivYdelta =
          aliasCached.fLinAliasData[j].X / aliasCached.fLinAliasData[j].YdataDelta;
    }
    fCachedAliasTable.push_back(aliasCached);
  }
}

void VecKleinNishinaComptonModel::SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td)
{

  EMModel::SampleSecondariesVector(tracks, td);
}

void VecKleinNishinaComptonModel::SampleSecondariesVectorAlias(LightTrack_v &tracks, geant::TaskData *td)
{
}

void VecKleinNishinaComptonModel::SampleReducedPhotonEnergyVec(const double *egamma, const double *r1, const double *r2,
                                                               const double *r3, double *out, int N)
{

  for (int i = 0; i < N; i += kPhysDVWidth) {
    // determine electron energy lower grid point
    PhysDV legamma;
    vecCore::Load(legamma, &egamma[i]);
    legamma = vecCore::math::Log(legamma);
    //
    PhysDV val        = (legamma - fSTLogMinPhotonEnergy) * fSTILDeltaPhotonEnergy;
    PhysDI indxEgamma = (PhysDI)val; // lower electron energy bin index
    PhysDV pIndxHigh  = val - indxEgamma;
    PhysDV R1;
    vecCore::Load(R1, &r1[i]);
    PhysDM mask = R1 < pIndxHigh;
    if (!mask.isEmpty()) {
      vecCore::MaskedAssign(indxEgamma, mask, indxEgamma + 1);
    }

    PhysDV xiV;
    for (int l = 0; l < kPhysDVWidth; ++l) {
      int idx             = (int)vecCore::Get(indxEgamma, l);
      LinAliasCached &als = fCachedAliasTable[idx];
      // sample the transformed variable xi=[\alpha-ln(ep)]/\alpha (where \alpha=ln(1/(1+2\kappa)))
      // that is in [0,1] when ep is in [ep_min=1/(1+2\kappa),ep_max=1] (that limits comes from energy and momentum
      // conservation in case of scattering on free electron at rest).
      // where ep = E_1/E_0 and kappa = E_0/(mc^2)
      int aliasBin = fAliasSampler->SampleLinearGetIndex(r2[i + l], fSTNumDiscreteEnergyTransferVals,
                                                         als.fAliasW.data(), als.fAliasIndx.data());
      LinAliasCached::LinAliasData &alsData = als.fLinAliasData[aliasBin];
      double xi = fAliasSampler->SampleLinearCached(alsData.X, alsData.Xdelta, alsData.Ydata, alsData.YdataDelta,
                                                    r3[i + l], alsData.XdivYdelta);

      //      const LinAlias *als = fSamplingTables[idx];
      //      // sample the transformed variable xi=[\alpha-ln(ep)]/\alpha (where \alpha=ln(1/(1+2\kappa)))
      //      // that is in [0,1] when ep is in [ep_min=1/(1+2\kappa),ep_max=1] (that limits comes from energy and
      //      momentum
      //      // conservation in case of scattering on free electron at rest).
      //      // where ep = E_1/E_0 and kappa = E_0/(mc^2)
      //      const double xi = fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]),
      //                                                    &(als->fAliasIndx[0]), fSTNumDiscreteEnergyTransferVals,
      //                                                    r2[i+l], r3[i+l]);
      // transform it back to eps = E_1/E_0
      vecCore::Set(xiV, l, xi);
    }
    // transform it back to eps = E_1/E_0
    // \epsion(\xi) = \exp[ \alpha(1-\xi) ] = \exp [\ln(1+2\kappa)(\xi-1)]
    PhysDV kappa;
    vecCore::Load(kappa, &egamma[i]);
    kappa = kappa / geant::units::kElectronMassC2;

    PhysDV outV = vecCore::math::Exp(vecCore::math::Log(1. + 2. * kappa) * (xiV - 1.));
    vecCore::Store(outV, &out[i]);
  }
}
}
