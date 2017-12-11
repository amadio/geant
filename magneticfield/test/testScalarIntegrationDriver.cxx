//
//  Test ScalarIntegrationDriver 
//   * compares with the output of a reference stepper (high accuracy) - ok for small steps
// 
//  Based on testStepperFixed.cc
//    which was started from the work of Somnath Banerjee in GSoC 2015
//

#include <cfloat>
#include <iomanip>

// #include "G4UniformMagField.hh"
// #include "G4SystemOfUnits.hh"
#include "Units.h"
#include <cfloat>

// using fieldUnits::meter;
using fieldUnits::millimeter;   
using fieldUnits::second;  
using fieldUnits::eplus;  
using fieldUnits::tesla;
using fieldUnits::degree;

// #include "ThreeVector.h"
#include "base/Vector3D.h"
typedef vecgeom::Vector3D<double>  ThreeVector;

#include "ScalarUniformMagField.h"

#include "ScalarMagFieldEquation.h"

#include "VScalarIntegrationStepper.h"
#include "StepperFactory.h"

#include "ScalarFieldTrack.h"
#include "ScalarIntegrationDriver.h"


// #define  COMPARE_TO_G4  1

#ifdef COMPARE_TO_G4
//#include <system.h>
//#include "G4Types.h"
#include "G4UniformMagField.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ClassicalRK4.hh"
#endif

using namespace std;
// using namespace CLHEP;

int main(int argc, char *args[])
{
    constexpr unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

    using  GvEquationType=  ScalarMagFieldEquation<ScalarUniformMagField, Nposmom>;
 
    /* -----------------------------SETTINGS-------------------------------- */
    
    /* Parameters of test
     - Modify values  */
    
    int no_of_steps = 20;         // No. of Steps for the stepper
    int stepper_no =  5;         // Choose stepper no., for refernce see above
    double step_len_mm = 200.;    // meant as millimeter;  //Step length 
    double z_field_in = DBL_MAX;
    double epsTolInp =  -1.0; 
    
    //Checking for command line values :
    if(argc>1)
        stepper_no = atoi(args[1]);
    if(argc > 2)
       step_len_mm = (float)(stof(args[2]));   // *mm);
    if(argc > 3)
        no_of_steps = atoi(args[3]);
    if(argc > 4)
       z_field_in = (float) (stof(args[4]));     // tesla
    if(argc > 5)
       epsTolInp  = (float) (stof(args[5]));     // number 

    double step_len = step_len_mm * fieldUnits::millimeter;
    
    //Set Charge etc.
    double particleCharge = +1.0;      // in e+ units

    //Choice of output coordinates
    int
    columns[] =
    {
       1 , 1 , 0 ,  // position  x, y, z 
       1 , 1 , 0 ,  // momentum  x, y, z
       0 , 0 , 0 ,  // dydx pos  x, y, z
       0 , 0 , 0    // dydx mom  x, y, z
    }; //Variables in yOut[] & dydx[] we want to display - 0 for No, 1 for yes

    bool printErr= 1,   // Print the predicted Error
         printRef= 0,   // Print the reference Solution 
         printDiff= 1;  // Print the difference 
    bool printSep = 0;  // separator  '|'
    int  printOk  = 3;  // "Ok" or "Bad" - for Vec (val>0) & for each component (val>2)
                        //   Format:   1 Vec-Terse:    "Ok" Or "Bad"
                        //   Format:  >1 Vec-Verbose:  "0 Ok" "1 Bad"    
                        //            >3 [i]-Verbose: 
    bool printRatio    = 0;  //  yDiff / yAverage - per component
    bool printRatioVec = 1;  //                   - per vector  (ie Pos or Vector )
    bool printMagShift = 1;  //  Print  | p | / | p_0 | - 1.0
    bool printInp = 0;  // print the input values
    bool printInpX= 0;  // print the input values for Ref 

    const int    epsDigits = 3;       //  Digits - for printing the different
    const double errAmplif = 1.0e5;   //  Amplification for printing of error

    const unsigned int wd= 12;  // Width for float/double
    const unsigned int wdf = wd + epsDigits - 3;
    const unsigned int wok = 5; //  2 for (0/1," ") + 3 for "Bad"/"Ok "

    //Set coordinates here
    double
       x_pos = 0.,                 //pos - position  : input unit = mm
       y_pos = 0.,
       z_pos = 0.;
    double   
       x_mom = 0.0,                 //mom - momentum  : input unit = GeV / c
       y_mom = 1.0,
       z_mom = 1.0;
    double
       x_field = 0.,               //Uniform Magnetic Field (x,y,z)
       y_field = 0.;               //  Tesla // *tesla ;
    double z_field;
    if( z_field_in < DBL_MAX )
       z_field = z_field_in;
    else
       z_field = -1.0;  //  Tesla // *tesla ;

    // Field
    auto gvField= new ScalarUniformMagField( fieldUnits::tesla * ThreeVector(x_field, y_field, z_field) );

    cout << "#  Initial  Field strength (GeantV) = "
         << x_field << " , " << y_field << " , " << z_field 
       // << (1.0/fieldUnits::tesla) * gvField->GetValue()->X() << ",  "
         << " Tesla " << endl;
    cout << "#  Initial  momentum * c = " << x_mom << " , " << y_mom << " , " << z_mom << " GeV " << endl;
    //Create an Equation :
    auto gvEquation =
       new GvEquationType(gvField);
       // new ScalarMagFieldEquation<ScalarUniformMagField, Nposmom>(gvField);
    // gvEquation->InitializeCharge( particleCharge );

    /*-------------------------PREPARING STEPPER-----------------------------*/

    //Create a stepper :
    VScalarIntegrationStepper *myStepper; // , *exactStepper;
    // G4MagIntegrationStepper *g4refStepper;    

    const int cloneBump= 10;
    bool useClonedStepper= (stepper_no > cloneBump);
    if(  useClonedStepper )
       stepper_no -= cloneBump;

    myStepper= StepperFactory::CreateStepper<GvEquationType>(gvEquation, stepper_no);
    //         *****************************
    if( useClonedStepper ){
       auto cloneStepper = myStepper->Clone();
       delete myStepper;
       myStepper = cloneStepper;
    }

    /*-------------------------PREPARING DRIVER -----------------------------*/
    const double hminimum  = 0.0001 * millimeter;  // Minimum step = 0.1 microns
    const double epsTolDef = 1.0e-5;              // Relative error tolerance
    int   statisticsVerbosity= 1;

    const double epsTol =  ( epsTolInp < 0.0 ) ? epsTolDef : epsTolInp; 
    cout << "#  Driver parameters:  eps_tol= "  << epsTol << "  h_min= " << hminimum << endl;
 
    auto integrDriver= new ScalarIntegrationDriver( hminimum,
                                                    myStepper,
                                                    Nposmom,
                                                    statisticsVerbosity); 
    // myStepper->InitializeCharge( particleCharge );
    //integrDriver->InitializeCharge( particleCharge );
 
    //Initialising coordinates
    const double mmGVf = fieldUnits::millimeter;
    const double ppGVf = fieldUnits::GeV ;  //   it is really  momentum * c_light
                                         //   Else it must be divided by fieldUnits::c_light;
    // const double ppGVf = fieldUnits::GeV / Constants::c_light;     // OLD

    // double yIn[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
    double yIn[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                    x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};
    
#if COMPARE_TO_G4
    double mass = 1;

    const double mmG4 = CLHEP::millimeter;
    const double ppG4 = CLHEP::GeV ;  //  In G4 too 'p' means p*c -- so no division  / CLHEP::c_light;

    const double mmRef = mmG4;  // Unit for reference of lenght   - milli-meter
    const double ppRef = ppG4;  // Unit for reference of momentum - GeV / c^2
    
    // G4double yInX[] = {x_pos * mmG4, y_pos * mmG4 ,z_pos * mmG4,
    //                   x_mom * ppG4 ,y_mom * ppG4 ,z_mom * ppG4};

    G4UniformMagField myField( CLHEP::tesla * G4ThreeVector(x_field, y_field, z_field));
    cout << " Field strength = " << (1.0/CLHEP::tesla) * myField << " Tesla " 
         << endl;
    
    G4Mag_UsualEqRhs *g4Equation = new G4Mag_UsualEqRhs(&myField);
    
    G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                              spin=0.0,
                              magneticMoment=0.0);
    
    g4Equation->SetChargeMomentumMass( chargeState,
                                       G4ThreeVector(x_mom, y_mom, z_mom).mag(), //momentum magnitude
                                       mass);  // unused
    auto g4exactStepper =
                           new G4ClassicalRK4(g4Equation);
                        // new G4ExactHelixStepper(g4Equation);
    auto exactStepper = g4ExactStepper;
#else
    // double yInX[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
    //                 x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};    
    const double mmRef = mmGVf; // Unit for reference of lenght   - milli-meter
    const double ppRef = ppGVf; // Unit for reference of momentum - GeV / c^2
    
    auto gvEquation2 = new GvEquationType(gvField);
                   // new ScalarMagFieldEquation<ScalarUniformMagField, Nposmom>(gvField);
    // gvEquation2->InitializeCharge( particleCharge ); // Let's make sure
    
    // Should be able to share the Equation -- eventually
    // For now, it checks that it was Done() -- and fails an assert

    //Creating the baseline stepper
    auto exactStepperGV =
        new GUTCashKarpRKF45<GvEquationType,Nposmom>(gvEquation2);
        // new TClassicalRK4<GvEquationType,Nposmom>(gvEquation2);
    cout << "#  Reference stepper is: GUTCashKarpRKF45 <GvEquationType,Nposmom>(gvEquation2);" << endl;
    // cout << "#  Reference stepper is: TClassicalRK4 <GvEquationType,Nposmom>(gvEquation2);" << endl;    
       // new TSimpleRunge<GvEquationType,Nposmom>(gvEquation2);    
       // new GUExactHelixStepper(gvEquation2);

    // Configure Stepper for current particle
    //exactStepperGV->InitializeCharge( particleCharge ); // Passes to Equation, is cached by stepper
    // gvEquation2->InitializeCharge( particleCharge ); //  Different way - in case this works
    
    auto exactStepper = exactStepperGV;
#endif
    std::cout << "# step_len_mm = " << step_len_mm;
    std::cout << " mmRef= " << mmRef << "   ppRef= " << ppRef << std::endl;
    
    double yInX[] = {x_pos * mmRef, y_pos * mmRef ,z_pos * mmRef,
                     x_mom * ppRef ,y_mom * ppRef ,z_mom * ppRef};

    double stepLengthRef = step_len_mm * mmRef;
    
    //Empty buckets for results
    double dydx[8] = {0.,0.,0.,0.,0.,0.,0.,0.},  // 2 extra safety buffer
        dydxRef[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
           yout[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
          youtX[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
      //  yerr [8] = {0.,0.,0.,0.,0.,0.,0.,0.},
          yerrX[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
          yDiff[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
          yAver[8] = {0.,0.,0.,0.,0.,0.,0.,0.};          
    /*-----------------------END PREPARING STEPPER---------------------------*/


    /*---------------------------------------------------*/
    //        -> First Print the (commented) title header
    cout<<"\n#";
    cout<<setw(5)<<"Step";
    cout<<setw(9)<<" StepLen";    
    for (int i=0; i<6;i++) {
        if (columns[i])
        {
           if( printSep ) cout << " | " ;  // Separator
           if( printInp )
              cout << setw(wd-2)<< "yIn[" << i << "]";
           if( printInp )
              cout << setw(wd-2)<< "yInX[" << i << "]";           
           cout << setw(wd-2)<< "yOut[" << i << "]";
           if( printRef )
              cout << setw(wd-2) << "yOutX[" << i << "]";
           if( printErr )
              cout << setw(wd-1) << "yErrX[" << i << "]";
           if( printDiff )
              cout << setw(wd-2) << " Diff[" << i << "]";
           if( printOk>1  ) cout << " Cmp";
           if( printRatio ) cout << setw(wd) << " Ratio ";
        }
        if( (i+1)%3 == 0 ){
           if( printOk ){
              if( i == 2 ) cout << " PosCmp";
              if( i == 5 ) cout << " MomCmp";              
           }
           if( printRatioVec ) cout << setw(wd) << " Ratio ";
           if( (i == 5) && printMagShift ){
              cout << setw(wd) << " MomMagShift";
           }
        }
    }
    
    for (int i=0; i<6;i++){
        if (columns[i+6])
        {
           if( printSep ) cout << " | " ;  // Separator
           if( printInp )
              cout << setw(wd-2)<< "yIn[" << i << "]";           
           cout << setw(wd-2)<< "dydx[" << i << "]";
           if( printRef )
              cout << setw(wd-2) << "dydxRef[" << i << "]";
        }
    }
    cout<<setw(wd)<<"tan-1(y/x)";
    cout<<"\n";
    
    // Units
    const char *nameUnitLength= "mm";
    const char *nameUnitMomentum= "GeV/c";
    cout<<setw(6)<<"#Numbr";
    cout<<setw(9)<<nameUnitLength;
    for (int i=0; i<6;i++){
        if (columns[i])
        {
           const char* nameUnit = ( i<3 ) ? nameUnitLength : nameUnitMomentum ; 
           if( printSep ) cout << " | " ;  // Separator
           if( printInp ) cout << setw(wd-2)<< nameUnit;
           if( printInpX) cout << setw(wd-2)<< nameUnit;
           cout << setw(wd)<< nameUnit;
           if( printRef )  cout << setw(wd) << nameUnit;  // Reference output
           if( printErr ){
              cout << setw(4) << nameUnit;  // Estim. error
              cout << " * " << setw(wd-7) << errAmplif << " ";
           }
           if( printDiff ) cout << " " << setw(wdf-1) << nameUnit;    // Diff  new - ref
           if( printOk > 1  )  cout << "    ";
           if( printRatio ) cout << setw(wd) << "/epsTol ";                      
        }
        if( (i+1)%3 == 0 ){
           if( printOk )       cout << "  Ok?  ";
           // if( printRatioVec ) cout << setw(wd+4) << " Ratio ";
           if( printRatioVec ) cout << setw(wd) << "/epsTol ";                    
           if( (i == 5) && printMagShift ){
              // cout << setw(wd) << "  % ";
              cout << setw(wd) << "/epsTol ";                                    
           }
        }        
    }
    cout<<"\n";
    
    //-> Then print the data
    // cout.setf (ios_base::scientific);
    // cout.setf (ios_base::scientific);    
    cout.precision(3);

    const ThreeVector  startPosition( yIn[0], yIn[1], yIn[2]);
    const ThreeVector  startMomentum( yIn[3], yIn[4], yIn[5]);
    const double momentumMagInit = startMomentum.Mag();
    cout << "# momentumMagInit = " << momentumMagInit << endl;
    
    ScalarFieldTrack yStart( startPosition, startMomentum, particleCharge); 
    double total_step =0;
    /*----------------NOW STEPPING-----------------*/
    
    for(int j=0; j<no_of_steps; j++)
    {
        bool goodAdvance= true;

        cout<<setw(4)<<j ;           //Printing Step number

        // myStepper->RightHandSideVIS(yIn, dydx);               //compute dydx - to supply the stepper
        exactStepper->RightHandSideVIS(yInX, particleCharge, dydxRef);   //compute the value of dydx for the exact stepper

        // Driver begins at the start !
        ScalarFieldTrack  yTrackIn(  startPosition, startMomentum, particleCharge );  // yStart
        ScalarFieldTrack  yTrackOut( startPosition, startMomentum, particleCharge );  // yStart
        
        if( j > 0 )  // Let's print the initial points!
        {
           total_step += step_len;
           
           goodAdvance=
              integrDriver->AccurateAdvance( yTrackIn, total_step, epsTol, yTrackOut ); // , hInitial );
           // *****************************
           
           ThreeVector PositionOut = yTrackOut.GetPosition();
           yout[0]= PositionOut.x();
           yout[1]= PositionOut.y();
           yout[2]= PositionOut.z();
           ThreeVector MomentumOut = yTrackOut.GetMomentum();
           yout[3]= MomentumOut.x();
           yout[4]= MomentumOut.y();
           yout[5]= MomentumOut.z();
  
           // myStepper->StepWithErrorEstimate(yIn,dydx,step_len,yout,yerr);   //Call the 'trial' stepper

           // Compare with a high-quality stepper -- expect it works well for this step size (check!)
           //   This builds on its previous step to create the solution !

           // Call the reference stepper
#ifdef COMPARE_TO_G4        
           g4ExactStepper->Stepper(yInX,dydxRef,stepLengthRef,youtX,yerrX);
#else
           exactStepperGV->StepWithErrorEstimate(yInX,dydxRef,particleCharge,stepLengthRef,youtX,yerrX); 
#endif
        }
        if( goodAdvance ) cout << " o";  // OK
        else  cout << " x"; // Bad advance;

        // Check the results
        double sumDiffPos2= 0.0, sumDiffMom2= 0.0;  // Mag^2 |X-Xref|, |P-Pref|
        double maxDiffPos = 0.0, maxDiffMom = 0.0;
        for(int i=0; i<3;i++){
           yDiff[i] =        yout[i] / mmGVf - youtX[i] / mmRef ;
           yAver[i] = 0.5 * (yout[i] / mmGVf + youtX[i] / mmRef) ;

           sumDiffPos2 += ( yDiff[i] * yDiff[i] ); 
           maxDiffPos  = std::max( maxDiffPos, fabs( yDiff[i] ) );
        }
        double magDiffPos = sqrt( sumDiffPos2 );
        double relDiffPos =  
           (total_step > 0.0 ) ? magDiffPos / (total_step/mmGVf) : 0.0;
        bool     goodPos= relDiffPos < epsTol;
        double  ratioPos= relDiffPos / epsTol;

        // Check the momentum 
        ThreeVector  MomentumNew( yout[3], yout[4], yout[5] ); MomentumNew *= 1.0 / ppGVf;
        ThreeVector  MomentumRef( youtX[3], youtX[4], youtX[5] ); MomentumRef *= 1.0 / ppRef;
        double sumMom2new =0;
        for(int i=3; i<6;i++){
           yDiff[i] =        yout[i] / ppGVf - youtX[i] / ppRef ;
           yAver[i] = 0.5 * (yout[i] / ppGVf + youtX[i] / ppRef) ;

           sumMom2new  += yout[i] * yout[i]; 
           sumDiffMom2 += (yDiff[i] * yDiff[i]); 
           maxDiffMom  = std::max( maxDiffMom, fabs(yDiff[i]) );
        }
        ThreeVector  MomentumDiff( yDiff[3], yDiff[4], yDiff[5] ); MomentumRef *= 1.0 / ppRef;                 
        
        double magDiffMom = MomentumDiff.Mag(); // sqrt( sumDiffMom2 );
        double relDiffMom =
           (momentumMagInit > 0.0 ) ? magDiffMom / momentumMagInit : 0.0;
        bool     goodMom= relDiffMom < epsTol;
        double  ratioMom= relDiffMom / epsTol;
        double MomMagShift= sqrt(sumMom2new) / momentumMagInit - 1.0; // MomentumNew.Mag(); 

        static ThreeVector totErrMomRefV= ThreeVector( 0., 0., 0. );        
        ThreeVector errMomRefV= ThreeVector( yerrX[3], yerrX[4], yerrX[5] );
        totErrMomRefV += errMomRefV.Abs();
        
        //-> Then print the data
        cout.setf (ios_base::fixed);

        cout.precision(1);
        cout << setw(9) << total_step / mmGVf << " ";

        // Report Position        
        cout.precision(4);
        for(int i=0; i<3;i++){
           if(columns[i]){
              if( printSep ) cout << " | " ;  // Separator
              if( printInp ) cout << setw(wd-2)<< yIn[i] / mmGVf;
              if( printInpX ) cout << setw(wd-2)<< yInX[i] / mmRef;
              cout << setw(wd) << yout[i] / mmGVf ;
              // cout.precision(3);
              if( printRef ) cout<<setw(wd)<< youtX[i] / mmRef; // Reference Solution
              if( printErr ) cout<<setw(wd)<< errAmplif * (yerrX[i] / mmGVf) ;
              if( printDiff ){
                 cout.precision(epsDigits);
                 cout<< setw(wdf)<< yDiff[i];
              }               
              bool goodI =  fabs(yDiff[i]) <  epsTol * fabs(yAver[i]);
              if( printOk>2 ){
                 if( printOk > 3 ) cout << " " << setw(wok)<< (goodI ? "1 ok " : "0 bad" );
                 else              cout << " " << setw(3)<< (goodI ? "ok " : "bad" );
              }
              double ratioPi = (yAver[i] != 0.0) ? yDiff[i] / yAver[i] : 0.0;
              if( printRatio ) cout << " " << setw(wd-2) << ratioPi/epsTol << " ";
           }
        }
        if( printOk ) {
           // cout<< " Vec: ";
           if( printOk > 1 ) cout << " " << setw(wok)<< (goodPos ? "1 Ok " : "0 Bad" );
           else              cout << " " << setw(3)<< (goodPos ? "Ok " : "Bad" );
        }
        if( printRatioVec ){
           // if( printSep && ! ( columns[0] || columns[1] || columns[2] ) ) cout << " || "; 
           cout << " " << setw(wd) << ratioPos << " ";
        }
        cout.unsetf (ios_base::fixed);

        // Report momentum 
        for(int i=3; i<6;i++){
            if(columns[i]){
               cout.setf (ios_base::scientific);
               if( printSep ) cout << " | " ;  // Separator
               if( printInp ) cout << setw(wd-1)<< yIn[i] / ppGVf << " ";
               if( printInpX ) cout << setw(wd-1)<< yInX[i] / ppRef << " ";
               cout<<setw(wd) << yout[i] / ppGVf ;
               if( printRef ) cout<<setw(wd)<< youtX[i] / ppRef; // Reference Solution
               if( printErr ) cout<< " " << setw(wd)<< errAmplif * yerrX[i] / ppGVf ;
               if( printDiff ){
                  cout.precision(epsDigits);
                  cout<<setw(wdf)<< yDiff[i];
               }
               cout.unsetf (ios_base::scientific);
               if( printOk || printRatio ){
                  bool goodI =  fabs(yDiff[i]) <  epsTol * fabs(yAver[i]);
                  // if( printOk>1 ) cout<< " " << setw(3)<< (goodI ? "ok " : "bad" );
                  if( printOk>2 ){
                     if( printOk > 3 ) cout << " " << setw(wok)<< (goodI ? "1 ok " : "0 bad" );
                     else              cout << " " << setw(3)<< (goodI ? "ok " : "bad" );
                  }
                  double ratioMi = (yAver[i] != 0.0) ? yDiff[i] / yAver[i] : 0.0;
                  if( printRatio ){
                     // cout.unsetf (ios_base::scientific);
                     cout.setf (ios_base::fixed);                     
                     cout << " " << setw(wd) << ratioMi/epsTol << " ";
                     cout.unsetf (ios_base::fixed);
                     // cout.setf (ios_base::scientific);                     
                  }
               }
            }
        }
        // if( printOk )  cout << " vMem: " << setw(3)<< (goodMom ? "Ok " : "Bad" );
        if( printOk ) {
           // cout<< " Vec: ";
           cout << " ";
           if( printOk > 1 ) cout << goodMom; // << setw(wok)<< (goodMom ? "Ok " : "Bad" );
           cout  << " " << setw(3)<< (goodMom ? "Ok " : "Bad" );
           // else              cout << setw(3)<< (goodMom ? "Ok " : "Bad" );           
        }
        if( printRatioVec ){
           cout.setf (ios_base::fixed);
           cout << " " << setw(wd) << ratioMom << " ";
           cout.unsetf (ios_base::fixed);           
        }
        // cout << " m-m= " << setw(wd) << magDiffMom ;
        // cout << " max= " << setw(wd) << maxDiffMom ;
        // cout << " rel= " << setw(wd) << relDiffMom ; // Un-normalised to epsilon
        cout // << " MomMagShift% = " << setw(wd) << 100.0 * MomMagShift
             << setw(wd) << MomMagShift / epsTol;

        
        for(int i=0; i<6;i++)   // Print auxiliary components
        {
           double unitGVf=1;
           double unitRef=1;
           // if( i < 3 )             // length / length

           if( i >= 3 ){
              unitGVf = ppGVf / mmGVf; //  dp / ds
              unitRef = ppRef / mmRef; //  dp / ds
           }

           if( i == 0 ){            // length / length
              // cout.unsetf (ios_base::scientific);
              cout.setf (ios_base::fixed);
           }else if( i == 3 ){
              cout.unsetf (ios_base::fixed);              
              cout.setf (ios_base::scientific);
           }
           
           if(columns[i+6])
           {
               // cout << " dy/dx["<<i<<"] = ";
               if( printSep ) cout << " | " ;  // Separator
               // cout<<setw(wd)<< dydx[i] / unitGVf ; -- Not available from Driver
               // if( printRef ) 
               cout<<setw(wd)<< dydxRef[i] / unitRef; // Reference Solution
               cout.precision(epsDigits);               
               if( 0 ) {  //  ( printDiff ){  -- Not meaningful
                  cout<<setw(wdf)<< ( dydx[i] / unitGVf )
                                   - ( dydxRef[i] / unitRef );
               }
           }
           // if( i == 2 )     { cout.unsetf(ios_base::fixed);      }
           // else if ( i==5 ) { cout.unsetf(ios_base::scientific); }
        }
        cout.unsetf(ios_base::scientific);
        if( j > 0 )  // Step 0 did not move -- printed the starting values
        {
           // cout.unsetf(ios_base::scientific);
           cout.setf (ios_base::fixed);                         
           cout.precision(2);
           cout<<setw(wd) << atan2(yout[1],yout[0])/degree;
                             // atan(yout[1]/yout[0])/degree;
           if( printRef ) 
             cout<<setw(wd) << atan2(youtX[1],youtX[0])/degree;
           
           //Copy yout into yIn
           for(int i=0; i<6;i++){
              // yIn[i] = yout[i];
              yInX[i] = youtX[i];
           }
        }
        
        cout<<"\n";
    }
    
    /*-----------------END-STEPPING------------------*/

    /*------ Clean up ------*/
    
    delete myStepper;
    delete exactStepper;
    // delete gvEquation;  // The stepper now takes ownership of the equation
    // delete gvEquation2;    
    delete gvField;
    
    cout<<"\n\n#-------------End of output-----------------\n";
    
}
