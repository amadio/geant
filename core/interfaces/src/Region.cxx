
#include "Geant/Region.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

std::vector<Region *> Region::gTheRegionTable;

Region::Region(const std::string &name, bool iscutinlength, double gcut, double emcut, double epcut, double pcut)
    : fName(name), fIsProductionCutInLength(iscutinlength), fGammaCutValue(gcut), fElectronCutValue(emcut),
      fPositronCutValue(epcut), fProtonCutValue(pcut)
{
  fIndex = gTheRegionTable.size();
  gTheRegionTable.push_back(this);
}

} // namespace vecgeom
} // namespace VECGEOM_IMPL_NAMESPACE
