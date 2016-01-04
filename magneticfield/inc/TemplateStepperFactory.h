//===----------------------------------------------------------------------===//
/**
 * @file StepperFactory.h
 * @brief  Abstract field class for Geant-V prototype
 * @author John Apostolakis
 */
//===----------------------------------------------------------------------===//

#ifndef TEMPLATESTEPPER_FACTORY_H
#define TEMPLATESTEPPER_FACTORY_H 1

#include <ostream>

// Base types
#include "TemplateGUVEquationOfMotion.h"
#include "TemplateGUVIntegrationStepper.h"

// Concrete Types being created
#include "TemplateTSimpleRunge.h"
#include "TemplateTClassicalRK4.h"
#include "TemplateGUTCashKarpRKF45.h"

// namespace vecFieldPropagation {

/**
 * @brief Class StepperFactory
 */

template <class Backend>
class TemplateStepperFactory
{
   public:   
     static const unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec
     static const int DefaultStepperCode= 5;  // Cash Karp

      /**
       * @brief GeantTrack parametrized constructor
       *
       * @param EquationType - Type of Equation of Motion
       * @param equation     - Instance of Equaiton of Motion (type: EquationType)
       * @param StepperCode  - Integer Code to identify type of Stepper
       */
     template<typename EquationType>
        static
        TemplateGUVIntegrationStepper<Backend> *          
            CreateStepper(EquationType *equation, int StepperCode = DefaultStepperCode);

     // static TemplateStepperFactory* Instance();
};

template <class Backend>
template<typename EquationType>
TemplateGUVIntegrationStepper<Backend> *          
TemplateStepperFactory<Backend>
  ::CreateStepper(EquationType *equation, int StepperCode ) //Ananya: discuss : StepperCode : same for all or different? Should be same .... 
{
    TemplateGUVIntegrationStepper<Backend> *stepper; // , *exactStepper;

    const char *stepperName=0;
    const char * const NameSimpleRunge = "TSimpleRunge";
    const char * const NameClassicalRK4 = "TClassicalRK4";
    const char * const NameCashKarpRKF45 = "TCashKarpRKF45";

    int MaxStepperCode= 5;
    
    if( (StepperCode <= 0)
        || (StepperCode> MaxStepperCode)
        || (StepperCode == 2)  // Missing in range  min - max
        || (StepperCode == 3)
       )
       StepperCode= DefaultStepperCode;
    
    switch(StepperCode)
    {
      case 1: stepper = new TemplateTSimpleRunge<Backend,EquationType,Nposmom>(equation);
         stepperName= NameSimpleRunge;
         break;         
         // case 2: stepper = new G4SimpleHeum(equation);   break;
         // case 3: stepper = new BogackiShampine23(equation); break;
      case 4:
         stepper = new TemplateTClassicalRK4<Backend,EquationType,Nposmom>(equation);
         stepperName= NameClassicalRK4;
         break;
      case 5:
         stepper = new TemplateGUTCashKarpRKF45<Backend,EquationType,Nposmom>(equation);
         stepperName= NameCashKarpRKF45;
         break;
      default : stepper = (TemplateGUVIntegrationStepper<Backend>*) 0 ;
         std::cerr << " ERROR> TemplateStepperFactory: No stepper selected. " << endl;
         // exit(1); 
    }
    if( stepperName )
       std::cout << "TemplateStepperFactory: Chosen the  " << stepperName << " stepper." << std::endl;

    return stepper;
}

// } // end of namespace vecFieldPropagation
#endif
