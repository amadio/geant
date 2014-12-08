//
// $Id: NTSTFieldSetup.cc,v 1.2 2007-10-26 09:51:29 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//   User Field setup class implementation.
//

#include "NTSTFieldSetup.hh"
#include "NTSTFieldMessenger.hh"

//template version 
#include "TNTSTField.hh"
//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

bool fieldFlag = true;

NTSTFieldSetup::NTSTFieldSetup(G4MagneticField *pCommonField)
: fChordFinder(0),fEquation(0), fMagneticField(0),
    pAField1(0),pAField2(0),
    ffield(0),
    fStepper(0),fStepperType(4),
    fFieldName(0),fMinStep(0.01),fGradofField(0.000001)
{
    fMagneticField = pCommonField; 
    fFieldMessenger = new NTSTFieldMessenger(this) ; 

}

NTSTFieldSetup::NTSTFieldSetup()
:  fChordFinder(0),fEquation(0), fMagneticField(0),
    pAField1(0),pAField2(0),
    ffield(0),
    fStepper(0),fStepperType(4),
    fFieldName(0),fMinStep(0.01),fGradofField(0.000001)
{
    fMagneticField = new TUniformMagField( ThreeVector(0.0, 0.0, 0.0 ) );
    std::cout << " NTSTFieldSetup: magnetic field set to Uniform( 0.0, 0, 0 ) " << G4endl;
    InitialiseAll();
    fieldFlag = false;
}

    void
NTSTFieldSetup::InitialiseAll()
{
    fMinStep = 1.0*mm ; // minimal step of 1 mm is default
    fMaxEpsilon= 0.00001;
    fMinEpsilon= 0.001;
    fFieldManager = G4TransportationManager::GetTransportationManager()
        ->GetFieldManager();
    CreateStepperAndChordFinder();
}

////////////////////////////////////////////////////////////////////////////////

NTSTFieldSetup::~NTSTFieldSetup()
{
    // GetGlobalFieldManager()->SetDetectorField(0);
    //std::cout<<" Deleting NTSTFieldSetup"<<G4endl;
    //if(fMagneticField) delete fMagneticField;
    if(fChordFinder)   delete fChordFinder;
    if(fStepper)       delete fStepper;
    //std::cout<<"End of Deleting NTSTFieldSetup"<<G4endl;
}



/////////////////////////////////////////////////////////////////////////////
//
// Update field
//
#include "G4MagIntegratorDriver.hh"
void NTSTFieldSetup::CreateStepperAndChordFinder()
{
    SetField();  
    if(fEquation) delete fEquation;
    fEquation = new G4Mag_UsualEqRhs(fMagneticField);

    SetStepper();
    std::cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

    fFieldManager->SetDetectorField(fMagneticField );

    //    if(fChordFinder) delete fChordFinder;

    if(!fChordFinder)
      fChordFinder = new G4ChordFinder( fMagneticField, fMinStep,fStepper);

    fFieldManager->SetChordFinder( fChordFinder );

    fFieldManager->SetMinimumEpsilonStep( fMinEpsilon );    // Old value
    fFieldManager->SetMaximumEpsilonStep( fMaxEpsilon );      // FIX - old value
    fFieldManager->SetDeltaOneStep( 0.25 * mm );       // original value
    fFieldManager->SetDeltaIntersection( 0.10 * mm );  // original value

    std::cout << "Field Manager's parameters are " 
        << " minEpsilonStep= " << fFieldManager->GetMinimumEpsilonStep() << " "
        << " maxEpsilonStep= " << fFieldManager->GetMaximumEpsilonStep() << " " 
        << " deltaOneStep=   " << fFieldManager->GetDeltaOneStep() << " "
        << " deltaIntersection= " << fFieldManager->GetDeltaIntersection() 
        << G4endl;
    //G4MagInt_Driver *pDriver;
    //To have verbose from MagInt_Driver
    //fChordFinder->SetVerbose(1);
    //  pDriver=fpChordFinder->GetIntegratorDriver();
    //pDriver=new G4MagInt_Driver(fMinStep, 
    //                                   fStepper, 
    //                                   fStepper->GetNumberOfVariables(),2 );

    //fChordFinder->SetIntegrationDriver(pDriver);
    //fFieldManager->SetChordFinder( fpChordFinder );
    return;
}


/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

//=============test template mode================

#include "TUniformMagField.hh"
#include "TMagFieldEquation.hh"
#include "TCashKarpRKF45.hh"
#include "TClassicalRK4.hh"
#include "TExplicitEuler.hh"
#include "TSimpleRunge.hh"
#include "TSimpleHeum.hh"
//#include "TMagIntegratorDriver.hh"
#include "TNTSTTabulatedField3d.hh"

// this is ugly at moment 
#ifdef TABFIELD
typedef TNTSTTabulatedField3d Field_t;
#else
typedef TNTSTField            Field_t;
#endif
typedef TMagFieldEquation<Field_t, 6UL> Equation_t;
typedef TCashKarpRKF45<Equation_t, 6UL> StepperCashKarp_t;
typedef TClassicalRK4<Equation_t, 6UL>  StepperRK4_t;
typedef TSimpleHeum<Equation_t, 6UL>    StepperHeum_t;
typedef TSimpleRunge<Equation_t, 6UL>   StepperRunge_t;
typedef TExplicitEuler<Equation_t, 6UL> StepperExEuler_t;

//===============================================

void NTSTFieldSetup::SetStepper()
{
    //if(fChordFinder) delete fChordFinder;
  Equation_t* tEquation = new Equation_t( static_cast<Field_t*>(fMagneticField) );
    switch ( fStepperType ) 
    {
        case 0:  
            assert(fieldFlag);
            StepperExEuler_t *tStepperEuler;
            tStepperEuler = new StepperExEuler_t(tEquation);
            fStepper= tStepperEuler;

            // fChordFinder = 
            //     new TChordFinder
		    //        <Field_t, Equation_t, SteTMagInt_Driver<StepperExEuler_t> >
            //     (static_cast<Field_t*>(fMagneticField),
            //      fMinStep,
            //      tStepper);

            std::cout<<"TExplicitEuler is used."<<G4endl;     
            break;

        case 2:  
            StepperRunge_t *tStepperRunge;
            tStepperRunge = new StepperRunge_t(tEquation);
            // fChordFinder = 
            //     new TChordFinder
            //     <TMagInt_Driver<StepperRunge_t> >            
            //     (static_cast<Field_t*>(fMagneticField),
            //      fMinStep,
            //      tStepper);
            fStepper= tStepperRunge;
            std::cout<<"TSimpleRunge is used."<<G4endl;     
            break;
        case 3:  
            assert(fieldFlag); 
            StepperHeum_t *tStepper;
            tStepper = new StepperHeum_t(tEquation);
            // fChordFinder = 
            //     new TChordFinder
            //     <TMagInt_Driver<StepperHeum_t> >
            //     (static_cast<Field_t*>(fMagneticField),
            //      fMinStep,
            //      tStepper);
            fStepper= tStepper;
            std::cout<<"TSimpleHeum is called"<<G4endl;     
            break;
        case 4:  
            //assert(fieldFlag); 
            //StepperRK4_t *tStepper;
            //tStepper 
            fStepper = new StepperRK4_t(tEquation);
            //                fChordFinder = 
            //                    new TChordFinder
            //		//		  <Field_t, Equation_t, StepperRK4_t, TMagInt_Driver<StepperRK4_t> >
            //		  ( static_cast<Field_t*>(fMagneticField),
            //                     fMinStep,
            //                     tStepper);
    
            std::cout<<"TClassicalRK4 (default) is used."<<G4endl;     
            break;
        // case 5:  
        //    fStepper = new G4HelixExplicitEuler( fEquation ); 
        //    std::cout<<"G4HelixExplicitEuler is called"<<G4endl;     
        //    break;

        case 8:  
            fStepper  = new StepperCashKarp_t(tEquation);
            std::cout<<"TCashKarpRKF45 is used."<<G4endl;     
            break;

        default: 
            fStepper = 0;
            G4cerr<<"ERROR> no Stepper is set."<<G4endl;  
    }
    return; 
}

#include "NTSTGradientField.hh"
#include "TNTSTTabulatedField3d.hh"
void NTSTFieldSetup::SetField()
{
    switch(fFieldName)
    { 
        case 0: std::cout<<"Field set to  UniformField(default)"<<G4endl;break;

        case 1: if(pAField1)delete pAField1;
                    pAField1= new NTSTGradientField(fGradofField);
                fMagneticField=pAField1;
                std::cout<<"Field set to Gradient Field"<<G4endl;break;
        case 2: pAField2= new TNTSTTabulatedField3d("TableST5.dat", 0.0);
                fMagneticField=pAField2;
                std::cout<<"Field set to Tabulated Solenoid Field"<<G4endl;break;

        case 3: //ffield= new Field_t(1.5*tesla,0.,0.);
	  //       fMagneticField = new Field_t(1.5*tesla,0.,0.);
                //fMagneticField=ffield; 
                std::cout<<"Field set to  UniformField with new value"<<G4endl;break;
        default: std::cout<<"Field set to  UniformField(default)"<<G4endl;break;
    };

}


void NTSTFieldSetup::GetChordFinderStats()
{
    if(fChordFinder){
        fChordFinder->PrintStatistics();
    }
    else{
        std::cout<<"Copy of ChordFinder doesn't exist"<<G4endl;
    }

}

void NTSTFieldSetup::GetFieldCallStats()
{
    if(fChordFinder){
        switch(fFieldName)
        { 
	        case 0: std::cout<< ffield->GetCount() << G4endl;break;
     	    case 1: std::cout<< pAField1->GetCount()<<G4endl;pAField1->ClearCount();break;      
            case 2: std::cout<< pAField2->GetCount()<<G4endl;pAField2->ClearCount();break;
            case 3: std::cout<< ffield->GetCount()<<G4endl;ffield->ClearCount();break;
            default:std::cout<<"no field"<<G4endl;break;
        };
    } 
    else{
        std::cout<<"Field doesn't exist"<<G4endl;
    }  
}



////////////////////////////////////////////////////////////////////////////////
//
//  Utility method

G4FieldManager*  NTSTFieldSetup::GetGlobalFieldManager()
{
    return G4TransportationManager::GetTransportationManager()
        ->GetFieldManager();
}
