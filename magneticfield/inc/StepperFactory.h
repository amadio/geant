#ifndef STEPPER_FACTORY_H
#define STEPPER_FACTORY_H 1

#include <ostream>

// Base types
#include "GUVEquationOfMotion.h"
#include "GUVIntegrationStepper.h"

// Concrete Types being created
#include "TSimpleRunge.h"
#include "TClassicalRK4.h"
#include "GUTCashKarpRKF45.h"

// namespace vecFieldPropagation {

class StepperFactory
{
   public:   
     static const unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec
     static const int DefaultStepperType= 5;  // Cash Karp

     template<typename EquationType>
        static
        GUVIntegrationStepper *          
            CreateStepper(EquationType *equation, int StepperType = DefaultStepperType);

     // static StepperFactory* Instance();
};

template<typename EquationType>
  GUVIntegrationStepper *          
StepperFactory::CreateStepper(EquationType *equation, int StepperType )
{
    GUVIntegrationStepper *stepper; // , *exactStepper;

    const char *stepperName=0;
    const char * const NameSimpleRunge = "TSimpleRunge";
    const char * const NameClassicalRK4 = "TClassicalRK4";
    const char * const NameCashKarpRKF45 = "TCashKarpRKF45";

    int MaxStepperType= 5;
    
    if( (StepperType <= 0)
        || (StepperType> MaxStepperType)
        || (StepperType == 2)  // Missing in range  min - max
        || (StepperType == 3)
       )
       StepperType= DefaultStepperType;
    
    switch(StepperType)
    {
      case 1: stepper = new TSimpleRunge<EquationType,Nposmom>(equation);
         stepperName= NameSimpleRunge;
         break;         
         // case 2: stepper = new G4SimpleHeum(equation);   break;
         // case 3: stepper = new BogackiShampine23(equation); break;
      case 4:
         stepper = new TClassicalRK4<EquationType,Nposmom>(equation);
         stepperName= NameClassicalRK4;
         break;
      case 5:
         stepper = new GUTCashKarpRKF45<EquationType,Nposmom>(equation);
         stepperName= NameCashKarpRKF45;
         break;
      default : stepper = (GUVIntegrationStepper*) 0 ;
         std::cerr << " ERROR> StepperFactory: No stepper selected. " << endl;
         // exit(1); 
    }
    if( stepperName )
       std::cout << "StepperFactory: Chosen the  " << stepperName << " stepper." << std::endl;

    return stepper;
}

// } // end of namespace vecFieldPropagation
#endif
