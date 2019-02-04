#ifndef GEANT_CMS_Detector_Construction
#define GEANT_CMS_Detector_Construction

#include "Geant/UserFieldConstruction.h"

#include <string>
#include "Geant/Error.h"
#include "Geant/UserFieldConstruction.h"

class CMSmagField;

class CMSFieldConstruction : public geant::cxx::UserFieldConstruction
{
public:
  /** @brief Default constructor - field with fixed name */
  CMSFieldConstruction() : fFieldFilename(std::string("cmsmagfield2015.txt")), fCMSfield(nullptr) {}

  /** @brief Construct using field map in named file */
  CMSFieldConstruction(const char *fieldFilename) : fFieldFilename(fieldFilename), fCMSfield(nullptr) {}
  CMSFieldConstruction(std::string fieldFilename) : fFieldFilename(fieldFilename), fCMSfield(nullptr) {}

  /** @brief Construct for uniform field */
  CMSFieldConstruction(bool useUniform, double fieldValue=3.8*geant::units::tesla)
   : fUseUniformField( useUniform), fMagFieldValue(fieldValue) {};
  ~CMSFieldConstruction();

  /** @brief Destructor */
  void SetFileForField(const char *filename) { fFieldFilename = filename; }
  void SetFileForField(std::string filename) { fFieldFilename = filename; }

  /** @brief Method to register a B-field, and create integrator for it. */
  bool CreateFieldAndSolver( bool useRungeKutta,
                             VVectorField **fieldPP = nullptr
                             ) override final;


private:
  std::string fFieldFilename;
  bool             fUseUniformField = false;  //  If true, use uniform field 
  float            fMagFieldValue   = 3.8*geant::units::tesla;
  CMSmagField*     fCMSfield        = nullptr;
  UniformMagField* fUniformField    = nullptr;
  
  // ScalarUniformMagField*  fUniformField; // Alternative - for debugging only
  /** Field is created and owned by this class */

public:
};

#endif
