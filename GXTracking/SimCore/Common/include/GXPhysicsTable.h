#ifndef GXPhysicsTable_h
#define GXPhysicsTable_h

#include "GPConstants.h"

class GXPhysicsVector {

 public:
  FQUALIFIER GXPhysicsVector();
  FQUALIFIER GXPhysicsVector(GXPhysicsVector& right);
  FQUALIFIER ~GXPhysicsVector() {}

  FQUALIFIER G4double Value(double energy);
  FQUALIFIER G4double Energy(size_t index);
  FQUALIFIER bool SplinePossible();
  FQUALIFIER void ComputeSecDerivatives();
  FQUALIFIER void ComputeSecondDerivatives(G4double firstPointDerivative, double endPointDerivative);
  FQUALIFIER void FillSecondDerivatives();
  FQUALIFIER int FindBinLocation(G4double energy);
  FQUALIFIER G4double SplineInterpolation(double energy, int lastBin);
  FQUALIFIER G4double LinearInterpolation(double energy, int lastBin);
  FQUALIFIER G4double Interpolation(double energy, int lastBin);
  FQUALIFIER void SetSpline(bool val);

  //special treatment for the Inverse Range Table (G4LPhysicsFreeVector)
  FQUALIFIER G4double InvRangeValue(double energy);
  FQUALIFIER int InvRangeFindBinLocation(G4double energy);

  bool useSpline;
  bool isSecondDerivativeFilled;
  int type;   // The type of PhysicsVector (enumerator)
  int numberOfNodes;
  G4double edgeMin;           // Energy of first point
  G4double edgeMax;           // Energy of the last point
  
  G4double dataVector[maxPVDataVector];
  G4double binVector[maxPVDataVector];
  G4double secDerivative[maxPVDataVector];
  G4double dBin;          // Bin width - useful only for fixed binning
  G4double baseBin;       // Set this in constructor for performance
};


class GXPhysicsTable {

 public:
  FQUALIFIER GXPhysicsTable();
  FQUALIFIER ~GXPhysicsTable() {}

  FQUALIFIER void Print();

  int nPhysicsVector;
  GXPhysicsVector physicsVectors[maxPhysicsVector];
};

#endif
