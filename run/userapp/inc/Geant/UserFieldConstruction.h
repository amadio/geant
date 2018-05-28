//===--- UserFieldConstruction.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file         UserFieldConstruction.h
 * @brief        Base class for user's initialization of field & integrator
 *                 of tracks for this field.
 *
 * @description  Provides templated method CreateSolverAndField that does
 *                 the hard work - creating the concrete solver, using 
 *                 knowledge of the field's type.
 *               This version can be used create a uniform magnetic field.
 *
 * @author       John Apostolakis
 * @date         June 2017 - May 2018
 */
//===----------------------------------------------------------------------===//

#ifndef UserFieldConstruction_H
#define UserFieldConstruction_H 1

#include "base/Vector3D.h"
// #include "Geant/Typedefs.h"
// #include "base/Global.h"
#include "Geant/Error.h"

#include "Geant/FieldLookup.h"
#include "Geant/FieldPropagatorFactory.h"
#include "Geant/SystemOfUnits.h"

#include "Geant/UniformMagField.h" // For use in scalar and vector/flexible driver/stepper etc

// GEANT_DEVICE_DECLARE_CONV(class,UserFieldConstruction);

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class UserFieldConstruction {
public:
  // UserFieldConstruction();  // RootDL - i.e. moved below, where it is defined
  virtual ~UserFieldConstruction(){};
  // @brief   key virtual method - must create the field for detector & the integration driver ('solver')
  virtual bool CreateFieldAndSolver(bool /*useRungeKutta*/, VVectorField **fieldPP = nullptr); // RootDL
     
  /** Register a constanct B-field */
  // void UseConstantMagField( float value[3], const char* Unit= 0 ); // Default unit is kilogauss // RootDL

  /** Register a B-field, and create integrator for it. */
  // template <class Field_t>                   // RootDL
  // bool CreateSolverForField(Field_t* field); // RootDL

  void SetEpsilonRK(double val) { fEpsilonRK = val; }
  void SetMinimumStep(double v) { fMinimumStepInField = v; }

  // Inquiry methods
  bool IsFieldUniform() { return fUseUniformField; }
  double GetEpsilonRK() { return fEpsilonRK; }
  double GetMinimumStep() { return fMinimumStepInField; }
  bool IsFieldCreated() { return fCreatedField; }
  VVectorField *GetField() { return fpField; }

private:
  double fEpsilonRK;
  double fMinimumStepInField;
  vecgeom::Vector3D<float> fMagFieldValue;
  bool fUseUniformField;
  bool fZeroField;

protected:
  bool fCreatedField;
  bool fCalled;

protected:
  VVectorField *fpField;

public:
  static constexpr double fEpsilonDefault        = 3.0e-5;
  static constexpr double fMinimumStepInFieldDef = 1.0e-4; // GV units = cm
                                                           // vecgeom::Vector3D<float>  fMagFieldValueVec;

  // };   // RootComm

  // --> Changed to accomodate Root needs (problem in dictionary creation - Root version 6.10 )
public: // RootAdded
  // UserFieldConstruction::
  UserFieldConstruction()
      : fEpsilonRK(fEpsilonDefault), fMinimumStepInField(fMinimumStepInFieldDef), fUseUniformField(false),
        fZeroField(true), fCreatedField(false), fCalled(false), fpField(nullptr)
  {
  }

  template <class Field_t>
  bool
      // UserFieldConstruction::
         CreateSolverForField(Field_t *ptrField, int stepperTypeId= 0 )         
  {
    // printf(" -UserFieldConstruction::CreateSolverForField() called.\n");
     FieldPropagatorFactory::CreatePropagator<Field_t>(*ptrField, fEpsilonRK, fMinimumStepInField, stepperTypeId);
    fCreatedField = true;
    return true;
  }

  void
      // UserFieldConstruction:: // RootComm
      UseConstantMagField(float fieldVal[3], const char *Units = 0)
  {
    const char *methodName = "UserFieldConstruction::UseConstantMagField";
    bool defaultUsed       = false;
    double unit            = 1;

    if (Units == 0 || strcmp(Units, "kilogauss") == 0) {
      unit        = geant::units::kilogauss;
      defaultUsed = (Units == 0);
    } else if ((strcmp(Units, "gauss") == 0) || (strcmp(Units, "Gauss") == 0)) {
      unit = geant::units::gauss;
    } else if ((strcmp(Units, "tesla") == 0) || (strcmp(Units, "Tesla") == 0)) {
      unit = geant::units::gauss;
    } else {
      unit        = geant::units::kilogauss;
      defaultUsed = (Units == 0);
    }

    if (defaultUsed) printf("%s - WARNING: No units provided - using kilogauss as default unit", methodName);

    fMagFieldValue = vecgeom::Vector3D<float>(fieldVal[0] * unit, fieldVal[1] * unit, fieldVal[2] * unit);

    /*
    printf("%s called. Field value = %9.3g , %9.3g  %9.3g  kiloGauss\n",
           methodName,
           fMagFieldValue[0] / geant::units::kilogauss,
           fMagFieldValue[1] / geant::units::kilogauss,
           fMagFieldValue[2] / geant::units::kilogauss );
     */

    fUseUniformField = true;
    fZeroField       = (fMagFieldValue.Mag2() == 0.0);
  }
};

  /** @brief Create the global magnetic field and classes to integrate it. Register field. */
  /** @description  Must call the templated CreateSolverForField method.                   */
  // virtual
  inline bool
  UserFieldConstruction:: // RootComm     
  CreateFieldAndSolver(bool /*useRungeKutta*/, VVectorField **fieldPP  /* = nullptr */ )
  {
    static const char *method = "UserFieldConstruction::CreateFieldAndSolver";
    bool rtv                  = false;
    if (fieldPP) *fieldPP     = nullptr;

    geant::Print(method, "%s - method called.  Use uniform= %d  Value= %f %f %f - kiloggauss.  Zero-Flag= %d", method,
                 fUseUniformField, fMagFieldValue[0] / geant::units::kilogauss,
                 fMagFieldValue[1] / geant::units::kilogauss, fMagFieldValue[2] / geant::units::kilogauss, fZeroField);

    if (fUseUniformField) {
      auto gvUniformField = new UniformMagField(fMagFieldValue);
      bool isUniform      = true;
      auto fieldConfig    = new FieldConfig(gvUniformField, isUniform);
      FieldLookup::SetFieldConfig(fieldConfig);

      // printf("   Field class created - address= %p \n", gvUniformField );
      fpField = gvUniformField;

      // Check that field was correctedly created ...
      ThreeVector Position(0.0, 0.0, 0.1);
      ThreeVector fieldVal(0.0, 0.0, 0.13579);
      gvUniformField->GetFieldValue(Position, fieldVal);

      rtv = CreateSolverForField<UniformMagField>(gvUniformField);

      if (fieldPP) *fieldPP = gvUniformField; // Return it ??

      if (fZeroField) {
        geant::Print(method, " Zero Magnetic Field configured.");
      }
      fCalled = true;
    } else {
      fCalled = true;
      geant::Error(method, "No user Magnetic Field is registered.");
    }
    return rtv;
  }
// -- };  // RootComm

} // GEANT_IMPL_NAMESPACE
} // Geant
#endif
