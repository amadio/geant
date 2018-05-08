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
  /** @brief Destructor */
  CMSFieldConstruction() : fFieldFilename(std::string("cmsmagfield2015.txt")), fCMSfield(nullptr) {}
  // CMSFieldConstruction(const char* fieldFilename);
  // CMSFieldConstruction(std::string fieldFilename);
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
  CMSmagField *fCMSfield;
  // ScalarUniformMagField*  fUniformField; // Alternative - for debugging only
  /** Field is created and owned by this class */

public:
  // CMSFieldConstruction::
  CMSFieldConstruction(const char *fieldFilename) : fFieldFilename(fieldFilename), fCMSfield(nullptr) {}

  // CMSFieldConstruction::
  CMSFieldConstruction(std::string fieldFilename) : fFieldFilename(fieldFilename), fCMSfield(nullptr) {}
};

#endif
