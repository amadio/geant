#include "MaterialHandler.h"

#include "base/PhysicalConstants.h"
#include "materials/Material.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

MaterialHandler *MaterialHandler::fInstance = 0;

VECCORE_CUDA_HOST
MaterialHandler *MaterialHandler::Instance()
{
  if (fInstance == 0)
    fInstance = new MaterialHandler();
  return fInstance;
}

VECCORE_CUDA_HOST
MaterialHandler::MaterialHandler()
{
  // test mode: 0 all elements in the element table, 1 single element
  fElementMode = 0;

  // initialize the element array
  fNumberOfElements = 0;
  for (int i = 0; i < maximumZ; ++i)
    fElementArray[i] = 0;

  // build the element array
  BuildElementTable();
}

VECCORE_CUDA_HOST
MaterialHandler::~MaterialHandler() { ; }

VECCORE_CUDA_HOST
void MaterialHandler::BuildElementTable()
{
  // This should interface with the global material manager of GeantV so that
  // the element arrary is properly filled with all elements of detector
  // materials. Temporarily, build a table based on John's arrary

  constexpr int NumFx = 16;
  int Element[NumFx] = {82, 74, 8, 7, 6, 13, 18, 22, 26, 27, 30, 48, 54, 64, 79, 91};
  // Pb   W  O  N  C  Al  Ar, Ti  Fe  Cu  Zn  Cd  Xe  Gd  Au  Pa

  for (int ie = 0; ie < NumFx; ++ie)
    AddElement(Element[ie]);
}

VECCORE_CUDA_HOST
void MaterialHandler::AddElement(int element)
{
  // check validity of the element
  if (element > 0 && element < maximumZ) {
    // seach whether this element already exists
    bool found = false;
    for (int i = 0; i < fNumberOfElements; ++i) {
      if (fElementArray[i] == element) {
        found = true;
        break;
      }
    }

    // add a new element to the array
    if (!found) {
      fElementArray[fNumberOfElements] = element;
      fNumberOfElements++;
    }
  }
}

VECCORE_CUDA_HOST
void MaterialHandler::PrepareTargetElements(int *targetElements, int ntracks, int elementMode)
{
  // only two modes for now based on John's original method
  static int noCalls = 0;
  noCalls++;

  bool report = (noCalls == 1);

  if (elementMode == 0) { // all elements in the material table
    if (report)
      printf(" Generating Target Elements from table of %d elements - mode # =  %d\n", fNumberOfElements, elementMode);
    int indEl;
    for (int i = 0; i < ntracks; ++i) {
      indEl = (i % fNumberOfElements);
      targetElements[i] = fElementArray[indEl];
    }
  }
  else if (elementMode == 1) { // using a single element
    if (report)
      printf(" Using *Constant* Target Element Z = %d - mode # = %d\n", fElementArray[0], elementMode);

    for (int i = 0; i < ntracks; ++i) {
      targetElements[i] = fElementArray[0];
    }
  }
  else {
    printf(" Illeagal - mode # = %d\n", elementMode);
    assert(0);
  }
}

VECCORE_CUDA_HOST
void MaterialHandler::BuildMaterialTable()
{
  // build a material table - 3 materials for testing

  using CLHEP::g;
  using CLHEP::mole;
  using CLHEP::cm3;

  // should have an interface to add elements to Material
  double density = 1.0 * g / cm3;

  // PbWO4

  double PbWO4_a[3] = {207.2 * g / mole, 183.84 * g / mole, 16.00 * g / mole};
  double PbWO4_z[3] = {82, 74, 8};
  double PbWO4_w[3] = {0.455, 0.404, 0.141};

  density = 8.28 * g / cm3;
  new vecgeom::Material("PbWO4", PbWO4_a, PbWO4_z, PbWO4_w, 3, density);

  // SiO2
  double SiO2_a[2] = {28.09 * g / mole, 16.00 * g / mole};
  double SiO2_z[2] = {14, 8};
  double SiO2_w[2] = {0.467, 0.533};

  density = 2.200 * g / cm3;
  new vecgeom::Material("Quartz", SiO2_a, SiO2_z, SiO2_w, 2, density);

  // Brass
  double Brass_a[2] = {63.54 * g / mole, 65.41 * g / mole};
  double Brass_z[2] = {29, 30};
  double Brass_w[2] = {0.493, 0.507};

  density = 8.6 * g / cm3;
  new vecgeom::Material("Brass", Brass_a, Brass_z, Brass_w, 2, density);
}

VECCORE_CUDA_HOST
void MaterialHandler::PrepareMaterialIndex(int *materialIndex, int ntracks, int materialMode)
{
  // temporary material index array (this should be provided the nagivator)
  static int nCalls = 0;
  nCalls++;

  bool report = (nCalls == 1);

  std::vector<vecgeom::Material *> &mtable = vecgeom::Material::GetMaterials();
  int numberOfMaterials = mtable.size();

  if (materialMode == 0) { // all elements in the current material table
    if (report)
      printf(" Material Mode =  %d ( Generating material index from a table of %d materials )\n", materialMode,
             numberOfMaterials);
    for (int i = 0; i < ntracks; ++i) {
      materialIndex[i] = (i % numberOfMaterials);
    }
  }
  else if (materialMode == 1) { // using a single material
    if (report)
      printf(" Material Mode =  %d ( Using *Constant* target material mame = %s )\n", materialMode,
             (mtable[0]->GetName()));
    for (int i = 0; i < ntracks; ++i) {
      materialIndex[i] = (mtable[0])->GetIndex();
    }
  }
  else {
    printf(" Illeagal - mode # = %d\n", materialMode);
    assert(0);
  }
}

} // end namespace impl
} // end namespace vecphys
