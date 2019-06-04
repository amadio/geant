#include "CMSFieldConstruction.h"

#include "CMSmagField.h"

CMSFieldConstruction::~CMSFieldConstruction()
{
  delete fCMSfield;
}

bool CMSFieldConstruction::CreateFieldAndSolver( bool           useRungeKutta,
                                                 VVectorField **fieldPP)
{
  const char* methodName= "CMSFieldConstruction::CratedFieldAndSolver";
  // template using geant::UserFieldConstruction::CreateSolverForField<>;
  
  int stepperTypeId =   10;
  // Choices:  ( see FieldPrpagationFactory )
  // a) Using default Integration Driver  ( for now Simple Integration Driver )
  //  4 = Bogacki Shampine (3rd order, 4 stage FSAL method)
  //  6    = Cash Karp 4/5 order, 6 stage method  
  //  7    = Dormand Prince 5th order 7 stage method   ( = DoPri5 ) without FSAL (for now)
  // b) Using other Integration Driver - together with Dormand Prince stepper
  // 10    = Rolling Driver + Dormand Prince 5th order 7 stage method   ( = DoPri5 ) without FSAL (for now)  
  // 
  geant::Print("CMSFieldConstruction::CreateFieldAndSolver", " Called with Arg: useRungeKutta= ");
  if (useRungeKutta) {
    printf("On");
  } else {
    printf("Off");
  }
  printf("\n");
  
  if( fCMSfield ) {
     // Field is already initialised ... do not repeat !
     geant::Print("CMSFieldConstruction::CreateFieldAndSolver", " Second call - doing nothing.");
     if (fieldPP)
        *fieldPP = fCMSfield;
     return false;
  }
  
  if (fieldPP) {
     *fieldPP = nullptr;
  }

  // fUniformField= nullptr;
  geant::FieldConfig *fieldConfig = nullptr;
  bool rtv= false;
  VVectorField *fieldPtr; 

  if( fUseUniformField ) {
     auto gvUniformField = new UniformMagField(fMagFieldValue);
     fUniformField = gvUniformField;
     fCMSfield = nullptr;
     fieldConfig = new geant::FieldConfig(gvUniformField, fUseUniformField);
     printf("Calling Create Solver For Uniform Field\n");
     rtv = geant::UserFieldConstruction::
        CreateSolverForField<UniformMagField>(gvUniformField);
     fieldPtr= gvUniformField;
  }
  else
  {
     useRungeKutta = true; // Must initialize it always
     geant::Print(methodName, "useRungeKutta is ON (only available method for non-uniform field)" );
     //  "until 'general helix' is available ");

     std::cout << "    Calling CMSmagField constructor with filename= " << fFieldFilename << std::endl;
     fCMSfield = new CMSmagField(fFieldFilename);
     printf("Calling Create Solver For CMS Field\n");     
     fieldConfig = new geant::FieldConfig(fCMSfield,        false );
     rtv= geant::UserFieldConstruction::
          CreateSolverForField<CMSmagField>(fCMSfield, stepperTypeId);
     fieldPtr= fCMSfield;
  }
  geant::FieldLookup::SetFieldConfig(fieldConfig);

  if (fieldPP ) *fieldPP = fieldPtr;

  // fpField = fieldPtr; // UserFieldConstruction::SetField( fieldPtr );

  geant::Print("CMSFieldConstruction::CreateFieldAndSolver", "CMSmagfield created.");

  if (useRungeKutta) {
    printf("%s", "CMSFieldConstruction - Configured field propagation for Runge Kutta.");
  } else {
    printf("%s", "CMSFieldConstruction - NOT configuring field propagation with Runge Kutta.");
  }

  return true;
}
