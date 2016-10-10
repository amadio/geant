

#ifndef REGION_H
#define REGION_H

#include <vector>
#include <string>

namespace geantphysics {

class Material;
/**
 * @brief   Dummy class for regions.
 * @class   Region.
 * @author  M Novak, A Ribon
 * @date    june 2016
 */
class Region {
public:
  Region(const std::string &name, bool iscutinlength, double gcut, double emcut, double epcut);
  Region(const std::string &name, bool iscutinlength = true);

 ~Region(){}

  void AddMaterial(const Material *mat);
  std::vector<const Material*>& GetMaterialList() { return fMaterialList; }

  int     GetIndex() const { return fIndex; }
  const std::string& GetName() const {return fName;}

  bool    IsProductionCutInLength() const { return fIsProductionCutInLength; }
  double  GetGammaCut()    const { return fGammaCutValue; }
  double  GetElectronCut() const { return fElectronCutValue; }
  double  GetPositronCut() const { return fPositronCutValue; }

  // get the global Region table
  static const std::vector<Region*>& GetTheRegionTable() { return gTheRegionTable; }

private:
  std::string fName;
  int    fIndex; // index of this region
  bool   fIsProductionCutInLength;
  double fGammaCutValue;
  double fElectronCutValue;
  double fPositronCutValue;

  // list of materials added so far to this Region by AddMaterial method;
  // can store duplicated material pointers !
  std::vector<const Material*> fMaterialList;
  // the global Region table
  static std::vector<Region*>  gTheRegionTable;
};

} // namespace geantphysics

#endif // REGION_H
