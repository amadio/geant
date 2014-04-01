#ifndef cmsExpMagneticField_H
#define cmsExpMagneticField_H 1

#include "cmsExpNameSpace.hh"

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4TwoVector.hh"

#include <fstream>
using namespace std;

//class cmsExpTrackerBfield;

class cmsExpMagneticField : public G4MagneticField
{
public:
  cmsExpMagneticField();
  ~cmsExpMagneticField();

  void GetFieldValue(const double point[3], double *bField ) const;
  void SetFlux(G4double mag) { flux = mag ; }
  void SetFieldType(G4String newType);

  bool IsDefined(const double& z, const double& rho) const;
  void ReadFieldMap(const char* filename);
  void GetVolumeBaseBfield(G4double const *point, G4double *bField) const;
  void GetParametricBfield(G4double const *point, G4double *bField,
			   G4double rho) const;
  void GetUniformBfield(G4double const *point, G4double *bField,
			G4double rho) const;
  G4TwoVector** GetFieldMapPointer() const { return fieldMap; }

private:

  G4double flux;
  cmsExp::BFieldType fieldType;

  G4TwoVector **fieldMap;
  // cmsExpTrackerBfield *theTrakcerBfield; 

};

#endif

