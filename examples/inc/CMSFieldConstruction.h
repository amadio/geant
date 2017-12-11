#ifndef GEANT_CMS_Detector_Construction
#define GEANT_CMS_Detector_Construction

#include "UserFieldConstruction.h"

// Delete ASAP - if possible.    JA 2017.09.14
#ifdef  USE_ROOT_TObject
#ifndef ROOT_TObject
#include "TObject.h"
#endif
#endif

#include <string>
#include "Geant/Error.h"
#include "UserFieldConstruction.h"

class CMSmagField;
class GUVMagneticField;

class CMSFieldConstruction : public Geant::cxx::UserFieldConstruction 
#ifdef  USE_ROOT_TObject
                               , public TObject
#endif
{
  public:
    /** @brief Destructor */
    CMSFieldConstruction():  fFieldFilename(std::string("")), fCMSfield(nullptr) {}
    // CMSFieldConstruction(const char* fieldFilename);
    // CMSFieldConstruction(std::string fieldFilename);
    inline ~CMSFieldConstruction();

    /** @brief Destructor */
    void SetFileForField(const char *filename){ fFieldFilename= filename; }
    void SetFileForField(std::string filename){ fFieldFilename= filename; }

    /** @brief Method to register a B-field, and create integrator for it. */
    bool CreateFieldAndSolver(bool useRungeKutta= true,
                              GUVField** fieldPP= nullptr )  override  final;
  private:
    std::string   fFieldFilename;
    CMSmagField*  fCMSfield;
    // ScalarUniformMagField*  fUniformField; // Alternative - for debugging only
    /** Field is created and owned by this class */

    // ClassDef(CMSFieldConstruction, 1) // User application

// };

//  Implementations made inline, in order to cope with need to load dynamically,
//   using ROOT v6.
public:

// CMSFieldConstruction::
CMSFieldConstruction(const char* fieldFilename) :
  fFieldFilename(fieldFilename),
  fCMSfield(nullptr)
{}

// CMSFieldConstruction::
CMSFieldConstruction(std::string fieldFilename) :
  fFieldFilename(fieldFilename),
  fCMSfield(nullptr)
{}

// ClassImp(CMSFieldConstruction);
};

#include "CMSmagField.h"

CMSFieldConstruction::~CMSFieldConstruction() { delete fCMSfield; }    

bool
CMSFieldConstruction::
CreateFieldAndSolver(bool useRungeKutta, GUVField** fieldPP )
{
  using FieldType = CMSmagField;  

  Geant::Print("CMSFieldConstruction::CreateFieldAndSolver", " Called with Arg: useRungeKutta=");
  if(useRungeKutta )   { printf("on"); }  else { printf("Off"); }

  if( fieldPP ) {
     *fieldPP= nullptr;
  }

  std::cout << "    Calling CMSmagField constructor with filename= " << fFieldFilename << std::endl;
  fCMSfield= new CMSmagField(fFieldFilename);
  // fUniformField= nullptr;
  useRungeKutta= true;  // Must initialize it always -- 
  printf("CMSFieldConstruction::CratedFieldAndSolver> useRungeKutta - forced ON, until 'general helix' is available ");
  
  auto fieldPtr = fCMSfield;  

  if( fieldPP && fieldPtr ) *fieldPP= fieldPtr;

  fpField= fieldPtr;  // UserFieldConstruction::SetField( fieldPtr ); 
  
  Geant::Print("CMSFieldConstruction::CreateFieldAndSolver", "CMSmagfield created.");
  
  if(useRungeKutta){
    CreateSolverForField<FieldType>(fieldPtr);
    printf("%s", "CMSFieldConstruction - Configured field propagation for Runge Kutta.");    
  } else {
    printf("%s", "CMSFieldConstruction - NOT configuring field propagation with Runge Kutta.");
  }  

  return true;
}

// };  // Added for ROOT
#endif
