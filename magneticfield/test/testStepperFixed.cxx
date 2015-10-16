//
//  Compare the output of different steppers
// 
//  Based on testStepperFixed.cc
//    was the work of Somnath Banerjee in GSoC 2015
//

// #include "G4UniformMagField.hh"
// #include "G4SystemOfUnits.hh"
#include "Units.h"

// using fieldUnits::meter;
using fieldUnits::millimeter;   
using fieldUnits::second;  
using fieldUnits::c_light;
using fieldUnits::eplus;  
using fieldUnits::tesla;
using fieldUnits::degree;

// #include "G4ios.hh"

#include "TUniformMagField.h"

#include "TMagFieldEquation.h"

#include "GUVIntegrationStepper.h"
#include "TClassicalRK4.h"
// #include "GUCashKarpRKF45.h"

#include "GUExactHelixStepper.h"
// #include "GULineSection.hh"

// #include "BogackiShampine23.hh"
// #include "DormandPrince745.hh"
// #include "BogackiShampine45.h"
#include <iomanip>
// #include "G4SimpleHeum.hh"

//#include
//#include <system.h>
//#include "G4Types.h"

using namespace std;
// using namespace CLHEP;

//Version 2.0 - Includes direct comparison with ExactHelix

/* Stepper No.
 0. GUExactHelix
 4. TClassicalRK4
 5. TCashKarpRKF45 (to do)

Potential expansion:
 2. G4SimpleHeum
 3. BogackiShampine23
 6. BogackiShampine45
 7. DormandPrince745
 */

const unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

int main(int argc, char *args[])
{
    using  EquationType=  TMagFieldEquation<TUniformMagField, Nposmom>;
   
    /* -----------------------------SETTINGS-------------------------------- */
    
    /* Parameters of test
     - Modify values  */
    
    int no_of_steps = 3;         // No. of Steps for the stepper
    int stepper_no =  4;         // Choose stepper no., for refernce see above
    double step_len = 5.0 * millimeter;  //Step length 
    
    //Set coordinates here
    double
    x_pos = 0.,                 //pos - position
    y_pos = 0.,
    z_pos = 0.,
    
    x_mom = 0.,                 //mom - momentum
    y_mom = 10.,
    z_mom = 0.,
    
    x_field = 0.,               //Uniform Magnetic Field (x,y,z)
    y_field = 0.,
    z_field = -0.1*tesla ;
    
    //Set Charge etc.
    double particleCharge = +1.0 * eplus,     // in e+ units
    spin=0.0,                                   // ignore the spin
    magneticMoment= 0.0,                        // ignore the magnetic moment
    mass = 1;
    
    //Choice of output coordinates
    int
    columns[] =
    {
        1 ,  //x_pos
        1 ,  //y_pos
        0 ,  //z_pos
        0 ,  //x_mom
        0 ,  //y_mom
        0    //z_mom
    }; //Variables in yOut[] you want to display - 0 for No, 1 for yes
    
    /*----------------------------END-SETTINGS-------------------------------*/
    
    /************************************XXXXXX*****************************************/
    
    
    /*-------------------------PREPARING STEPPER-----------------------------*/
    
    /* CODER SPACE
     - don't modify values here */
    
    //Checking for command line values :
    if(argc>1)
        stepper_no = atoi(args[1]);
    if(argc > 2)
        step_len = (float)(stof(args[2])*millimeter);
    if(argc > 3)
        no_of_steps = atoi(args[3]);
    
    //Initialising coordinates
    double yIn[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
    double yInX[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
    
    //Empty buckets for results
    double dydx[8] = {0.,0.,0.,0.,0.,0.,0.,0.};   // Extra 2 as 'safety buffer'
    double dydxRef[8]= {0.,0.,0.,0.,0.,0.,0.,0.};  
    double yout[8]  = {0.,0.,0.,0.,0.,0.,0.,0.};
    double youtX[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
    double yerr[8]  = {0.,0.,0.,0.,0.,0.,0.,0.};
    double yerrX[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
    
    //1. Create a field :
    TUniformMagField*  gvField= new TUniformMagField( ThreeVector(x_field, y_field, z_field) );
    
    //Create an Equation :
    auto fEquation = new TMagFieldEquation<TUniformMagField, Nposmom>(gvField);

    fEquation->InitializeCharge( particleCharge );



#if COMPARE_TO_G4
    G4UniformMagField myField(G4ThreeVector(x_field, y_field, z_field));
    
    G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                              spin=0.0,
                              magneticMoment=0.0);
    
    g4Equation->SetChargeMomentumMass( chargeState,
                                     G4ThreeVector(x_mom, y_mom, z_mom).mag(), //momentum magnitude
                                     mass);//No place fo mass in fEquation though
#endif
    
    //Create a stepper :
    GUVIntegrationStepper *myStepper, *exactStepper;
    // G4MagIntegrationStepper *g4refStepper;    
    
    //Choose the stepper based on the command line argument
    switch(stepper_no){
      case 0: myStepper = new GUExactHelixStepper(fEquation); break;
       // case 2: myStepper = new G4SimpleHeum(fEquation);   break;
       // case 3: myStepper = new BogackiShampine23(fEquation); break;
      case 4: myStepper = new TClassicalRK4<EquationType,Nposmom>(fEquation); break;
       // case 5: myStepper = new TCashKarpRKF45(fEquation); break;         
       // case 6: myStepper = new BogackiShampine45(fEquation); break;
       // case 7: myStepper = new DormandPrince745(fEquation);  break;
      default : myStepper = 0 ;
    }
    
    //Creating the baseline stepper
    auto fEquation2 = new TMagFieldEquation<TUniformMagField, Nposmom>(gvField);
    // Should be able to share the Equation -- eventually
    // For now, it checks that it was Done() -- and fails an assert
    
    exactStepper = new GUExactHelixStepper(fEquation2);

    // Configure Stepper for current particle
    exactStepper->InitializeCharge( particleCharge ); // Passes to Equation, is cached by stepper

    /*-----------------------END PREPARING STEPPER---------------------------*/


    /*---------------------------------------------------*/
    //        -> First Print the (commented) title header
    cout<<"\n#";
    cout<<setw(7)<<"StepNo";
    for (int i=0; i<6;i++)
        if (columns[i])
            cout << setw(13)<< "yOut[" << i << "]"
            << setw(13) << "yErr[" << i << "]"
            << setw(13) << "yOut-yOutX[" << i << "]";
    cout<<setw(13)<<"tan-1(y/x)";
    
    
    //-> Then print the data
    cout<<"\n";
    
    cout.setf (ios_base::scientific);
    cout.precision(3);
    
    /*----------------NOW STEPPING-----------------*/
    
    for(int j=0; j<no_of_steps; j++)
    {
        cout<<setw(8)<<j + 1;           //Printing Step number
        
        myStepper->RightHandSide(yIn, dydx);               //compute dydx - to supply the stepper
        myStepper->StepWithErrorEstimate(yIn,dydx,step_len,yout,yerr);   //Call the 'trial' stepper
        
        exactStepper->RightHandSide(yInX, dydxRef);        //compute the value of dydx for the exact stepper
        
        exactStepper->StepWithErrorEstimate(yInX,dydxRef,step_len,youtX,yerrX); //call the reference stepper
        // g4exactStepper->Stepper(yInX,dydxRef,step_len,youtX,yerrX); //call the reference stepper
        
        //-> Then print the data
        cout.setf (ios_base::scientific);
        cout.precision(3);
        for(int i=0; i<6;i++)
            if(columns[i]){
                cout<<setw(15)<<yout[i]<<setw(15);
                cout<<setw(15)<<yerr[i];
                cout<<setw(15)<<yout[i] - youtX[i];
            }
        cout.unsetf(ios_base::scientific);
        cout.precision(6);
        cout<<setw(13)<<atan(yout[1]/yout[0])/degree;
        
        //Copy yout into yIn
        for(int i=0; i<6;i++){
            yIn[i] = yout[i];
            yInX[i] = youtX[i];
        }
        
        
        cout<<"\n";
    }
    
    /*-----------------END-STEPPING------------------*/

    /*------ Clean up ------*/
    fEquation->InformDone();    
    
    cout<<"\n\n#-------------End of output-----------------\n";
    
}
