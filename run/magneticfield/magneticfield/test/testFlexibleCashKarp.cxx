//
//  Compare the output of a selected steppers (scalar and vector.)
// 
//  Based on testStepperFixed.cc
//    which was part of the work of Somnath Banerjee in GSoC 2015
//
#include <iomanip>

#include "base/Vector3D.h"
#include <Geant/VectorTypes.h>

using namespace vecCore::math;

#include "UniformMagField.h"          // New type (universal) class
#include "MagFieldEquation.h"

// #include "VIntegrationStepper.h"   // Supressed ?
#include "CashKarp.h"

// ---- Future 'general' steppers (i.e. scalar + vector )
// #include "SimpleRunge.h"
// #include "ExactHelixStepper.h"    

#include "SystemOfUnits.h"

// using geant::meter;
using geant::millimeter;   
using geant::second;  
using geant::eplus;  
using geant::tesla;
using geant::degree;

// using namespace std;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;

template <typename Real_v, typename Stepper_t> // , typename Equation_t>
   bool TestFlexibleStepper(Stepper_t *stepper) ; // , Equation_t *equation);

/* -----------------------------SETTINGS-------------------------------- */
    
// Parameters of test - values can be modified here (or in arguments to main)
int    no_of_steps = 250;        // No. of Steps for the stepper
int    stepper_no =  5;          // Choose stepper no., for refernce see above
double step_len_mm = 200.;    // meant as millimeter;  //Step length 
double z_field_in = DBL_MAX;

bool debug= true;

/* ------------------ Parameters for PRINTING --------------------------- */
//Choice of output coordinates
int
columns[] =
    {
       1 , 1 , 1 ,  // position  x, y, z 
       1 , 1 , 1 ,  // momentum  x, y, z
       0 , 0 , 0 ,  // dydx pos  x, y, z
       1 , 1 , 0    // dydx mom  x, y, z
    }; //Variables in yOut[] & dydx[] we want to display - 0 for No, 1 for yes

bool printErr= 0,   // Print the predicted Error
   printRef= 0,   // Print the reference Solution 
   printDiff= 0;  // Print the diffrence 
bool printSep = 0;  // separator  '|'
bool printInp = 0;  // print the input values
// bool printInpX= 0;  // print the input values for Ref 

const unsigned int nwdf= 12;  // Width for float/double
// -------   End of Parameters for PRINTING ------------------------------ */

void processArguments(int argc, char *args[] )
{
    //Checking for command line values :
    if(argc>1)
        stepper_no = atoi(args[1]);
    if(argc > 2)
       step_len_mm = (float)(atof(args[2]));   // *mm);   //  stof( std::string, size_t ... )
    if(argc > 3)
        no_of_steps = atoi(args[3]);
    if(argc > 4)
       z_field_in = (float) (atof(args[4]));     // tesla
}

void printBanner();

constexpr unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

// ---------------  Main ------------------------------------------------ */

int main(int argc, char *args[])
{

    // using Backend = vecgeom::kVc ;
    // typedef typename Backend::precision_v Double_v;

    using Double_v  = Geant::Double_v;
    // template <typename Tpod>
    // using Vector3D  = vecgeom::Vector3D<Tpod>;
    using ThreeVector_f   = vecgeom::Vector3D<float>;
    using ThreeVector_d   = vecgeom::Vector3D<double>;
    // using ThreeVectorSimd = vecgeom::Vector3D<Double_v>;

    using  FieldType      =  UniformMagField;
    using  GvEquationType =  MagFieldEquation<FieldType>;
   

    if( debug ) cout << "Running with debug (progress) messages." << endl;

    processArguments( argc, args );
    
    double  x_field = 0., y_field = 0., z_field = 0.;  //Uniform Magnetic Field (x,y,z)
    z_field= ( z_field_in < DBL_MAX ) ? z_field_in :  -1.0;  //  Tesla // *tesla ;

    if( debug ) cout<<"---- Creating UniformMagField object"<<endl;
    // Vector3D<float>
    ThreeVector_f fieldValInput( x_field, y_field, z_field ); // ( 0.0, 0.0, -1.0);
    fieldValInput *= geant::tesla;
    
    // Field
    auto gvUniformField= new UniformMagField( fieldValInput ); // New class
    //  UniformMagField( geant::tesla * ThreeVector_f(x_field, y_field, z_field) );

    if( debug ) { 
       cout << "#  Initial  Field strength (GeantV) = "
            << x_field << " , " << y_field << " , " << z_field
            << " Tesla " << endl;
       
       ThreeVector_d origin(0.0, 0.0, 0.0), fieldValue;
       gvUniformField->GetFieldValue( origin, fieldValue );
       cout << "#    Values in object created       = "
            << (1.0/geant::tesla) * fieldValue.x() << ",  "
            << (1.0/geant::tesla) * fieldValue.y() << ",  "
            << (1.0/geant::tesla) * fieldValue.z() << endl;
       cout << endl;
    }
    
    //Create Equation :
    if( debug ) cout << "Create Equation" << endl;

    // 1. Original way of creating an equation
    auto   magEquation = new GvEquationType(gvUniformField);
       // new MagFieldEquation<UniformMagField>(gvUniformField);

    if( debug ) cout<<"----Equation instantiated. "<<endl;    

    //Create a stepper :
    if( debug ) cout<<"---- Preparing to create (Vector) CashKarpRKF45 Stepper "<<endl;

    // 2. Create stepper
    // using StepperType = VectorCashKarp< GvEquationType, Nposmom>;    
    using StepperType = CashKarp< GvEquationType, Nposmom>;

    StepperType  myStepper2(magEquation);
    auto myStepper = &myStepper2;
    // myStepper = new VectorCashKarp<GvEquationType,Nposmom>(gvEquation);
                                                                                  
    if( debug ) {
       cout << "---- constructed (flexible) CashKarp" << endl;
       cout << " Stepper  information: ptr= " << myStepper << endl;
       // cout << "       full " << *myStepper << endl;
    }
                                                              
    // Phase 1 - get it to work without cloning

    cout << " Testing scalar.      " << endl;
    bool okScalar    = TestFlexibleStepper<double, StepperType >(myStepper);
    // cout << " Testing Vec Float.   " << endl;  
    // bool okVecFloat  = TestFlexibleStepper<Float_v, StepperType >(myStepper); // , magEquation);
    cout << " Testing Vec Double . " << endl;    
    bool okVecDouble = TestFlexibleStepper<Double_v, StepperType >(myStepper); 

    bool good =  okScalar &&  // okVecFloat &&
       okVecDouble;

    // delete myStepper;  // Only if object was new'ed !!

    if( debug )  cout<<"----deletion of stepper done "<<endl;

#ifdef BASELINESTEPPER
    delete exactStepper;
#endif 

    // delete gvEquation;  // The stepper now takes ownership of the equation
    // delete gvEquation2;    
    delete gvUniformField;

    if( debug ) 
       cout << "\n\n#-------------End of all tests -----------------\n";
    
    return good;
}

const double mmGVf = geant::millimeter;
const double ppGVf = geant::GeV ;  //   it is really  momentum * c_light
                                        //   Else it must be divided by geant::c_light;

template < typename Real_v, typename scalar_t, int Ncomp >
void SelectOutput( const Real_v arrVec[Ncomp] , scalar_t arrLane[Ncomp], int lane )
{
   for(unsigned int i=0; i<Ncomp;i++) {   
      arrLane[i]   = vecCore::Get( arrVec[i], lane );
   }
}

// template < typename scalar_t, typename scalar_t, int Ncomp >
// void SelectOutput( const double arrVec[Ncomp] , scalar_t arrLane[Ncomp], int lane )
// template < typename double, typename double, int Ncomp >
template < int Ncomp >
void SelectOutput( const double arrVec[Ncomp] , double arrLane[Ncomp], int /*lane*/ )
{
   for(unsigned int i=0; i<Ncomp;i++) {
      arrLane[i]   =  arrVec[i];
   }
}

template <typename Real_v, typename Stepper_t> // , typename Equation_t>
bool TestFlexibleStepper(Stepper_t *stepper) // , Equation_t *equation)
{
    // Initialising 'parameters'
    Real_v  chargeVec(-1.);
    double stepLengthValue = step_len_mm * geant::millimeter;

    // auto scratch = stepper->ObtainScratchSpace<Double_v>();
        // Scratch space needed for this thread / task ...
    
    // Initial coordinates , momentum
    double  x_pos = 0., y_pos = 0., z_pos = 0.; //pos - position  : input unit = mm
    double  x_mom = 0., y_mom = 1., z_mom = 1.;    //mom - momentum  : input unit = GeV / c

    if( debug ) {
       cout << "TestFlexibleStepper called again ----------  " << endl;
       cout << "#  Initial  momentum * c = " << x_mom << " , " << y_mom << " , " << z_mom << " GeV " << endl;
    }
    // Double_v yIn[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
    //                  x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};
    Real_v yInVec[]= { Real_v(x_pos * mmGVf), Real_v(y_pos * mmGVf) ,
                      Real_v(z_pos * mmGVf),
                      Real_v(x_mom * ppGVf), Real_v(y_mom * ppGVf) ,
                      Real_v(z_mom * ppGVf) };
    if( debug ) cout << "Initialized yIn: values [0]= " << yInVec[0] << endl;

    // double yInX[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
    //                 x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};    
    
    // auto gvEquation2 = new GvEquationType(gvUniformField);
                   // new TMagFieldEquation<ScalarUniformMagField, Nposmom>(gvUniformField);

    
    // Should be able to share the Equation -- eventually
    // For now, it checks that it was Done() -- and fails an assert

    //Empty buckets for results
    Real_v dydxVec[8] = { Real_v(0.),0.,0.,0.,0.,0.,0.,0.},  // 2 extra safety buffer
           yOutVec[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
           yErrVec[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
    
#ifdef BASELINESTEPPER
    const double mmRef = mmGVf; // Unit for reference of lenght   - milli-meter
    const double ppRef = ppGVf; // Unit for reference of momentum - GeV / c^2
    
    //Creating the baseline stepper
    auto exactStepperGV =
        new TClassicalRK4<GvEquationType,Nposmom>(gvEquation2);
    cout << "#  Reference stepper is: TClassicalRK4<GvEquationType,Nposmom>(gvEquation2);" << endl;

       // new TSimpleRunge<GvEquationType,Nposmom>(gvEquation2);    
       // new GUExactHelixStepper(gvEquation2);

    auto exactStepper = exactStepperGV;

    Real_v dydxVecRef[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
             yOutVecX[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
             yErrVecX[8] = {0.,0.,0.,0.,0.,0.,0.,0.};    
    cout << " mmRef= " << mmRef << "   ppRef= " << ppRef << endl;

    Real_v yInX[] = {x_pos * mmRef, y_pos * mmRef ,z_pos * mmRef,
                      x_mom * ppRef ,y_mom * ppRef ,z_mom * ppRef};

    // Simple arrays for outputing and/or checking values
    double yinX[Nposmom], youtX[Nposmom], yerrX[Nposmom];
#endif        
    cout << "# step_len_mm = " << step_len_mm;

    // double stepLengthRef = step_len_mm * mmRef;
    /*-----------------------END PREPARING STEPPER---------------------------*/

    printBanner();
    
    // Units
    const char *nameUnitLength= "mm";
    const char *nameUnitMomentum= "GeV/c";
    cout<<setw(6)<<"#Numbr";
    for (int i=0; i<6;i++)
        if (columns[i])
        {
           if( printSep ) cout << " | " ;  // Separator                      
           const char* nameUnit = ( i<3 ) ? nameUnitLength : nameUnitMomentum ; 
           cout << setw(nwdf)<< nameUnit;
           if( printRef )  cout << setw(nwdf) << nameUnit;  // Reference output
           if( printErr )  cout << setw(nwdf) << nameUnit;  // Estim. error
           if( printDiff ) cout << setw(15) << nameUnit;    // Diff  new - ref
        }    
    cout<<"\n";
    
    //-> Then print the data

    
    // cout.setf (ios_base::scientific);
    // cout.setf (ios_base::scientific);    
    cout.precision(3);
    
    /*----------------NOW STEPPING-----------------*/
    no_of_steps = 25;
    for(int j=0; j<no_of_steps; j++)
    {
        cout<<setw(6)<<j ;           //Printing Step number
        Real_v  step_len( stepLengthValue ); 
        stepper->RightHandSideInl(yInVec, chargeVec, dydxVec);  //compute dydx - to supply the stepper
#ifdef  BASELINESTEPPER
        exactStepper->RightHandSideVIS(yInVecX, dydxVecRef);   //compute the value of dydx for the exact stepper
#endif

        int lane=0; // Lane for printing

        using scalar_t  = double;
        scalar_t  yin[Nposmom] = { -999., -99.9, -9.99, 10.0, 1.0, 0.1 };
        scalar_t  yout[Nposmom], dydx[Nposmom], yerr[Nposmom];                  
        if( j > 0 )  // Do nothing for j=0, so we can print the initial points!
        {
           stepper->StepWithErrorEstimate( yInVec, dydxVec, chargeVec, step_len, yOutVec, yErrVec /*, scratch */ );
           //       *********************
#ifdef  BASELINESTEPPER
           exactStepperGV->
              StepWithErrorEstimate(yInVecX,dydxVecRef,chargeVec,stepLengthRef,yOutVecX,yErrVecX /* ,scratchRef*/ );
           // *********************
#endif
        }

        // Try1 : 
        // SelectoOutput( yInVec,  yin,  lane);        
        // Try2 : 
        // void CopyToArr(Real_v realArr[], scalar_t scalArr[], int lane ) {
        // SelectOutput<Real_v, scalar_t, N>(realArr, scalArr, lane); }
        // CopyToArr( yInVec,  yin,  lane);

        // Select one lane for printing.

        // General code - for vecCore types or built-in types
        SelectOutput<Real_v, scalar_t, Nposmom> ( yInVec,  yin,  lane);
        if( j > 0 )
           SelectOutput<Real_v, scalar_t, Nposmom>  ( yOutVec, yout, lane);
        else
           SelectOutput<Real_v, scalar_t, Nposmom>  ( yInVec,  yout, lane);
        SelectOutput<Real_v, scalar_t, Nposmom>  ( yErrVec, yerr, lane);
        SelectOutput<Real_v, scalar_t, Nposmom>  ( dydxVec, dydx, lane);

#ifdef BASELINESTEPPER
        SelectOutput<Real_v, scalar_t, Nposmom> ( yInVecX,  yinX,  lane);
        if( j > 0 )
           SelectOutput<Real_v, scalar_t, Nposmom>  ( yOutVecX, youtX, lane);
        else
           SelectOutput<Real_v, scalar_t, Nposmom>  ( yInVecX, youtX, lane);
        SelectOutput<Real_v, scalar_t, Nposmom>  ( yErrVecX,   yerrX, lane);
        SelectOutput<Real_v, scalar_t, Nposmom>  ( dydxVecRef, dydxX, lane);        
#endif 
        // -->> End of Selecting 'lane' for printing
        
        //-> Then print the data
        cout.setf (std::ios_base::fixed);
        cout.precision(4);
        for(int i=0; i<3;i++)
            if(columns[i]){
               if( printSep ) cout << " | " ;  // Separator
               if( printInp ) cout << setw(nwdf-2)<< yin[i] / mmGVf;
              #ifdef BASELINESTEPPER
               if( printInpX ) cout << setw(nwdf-2)<< yinX[i] / mmRef;
               if( printRef ) cout<<setw(nwdf)<< youtX[i] / mmRef; // Reference Solution
               if( printDiff )                
                  cout<<setw(nwdf)<< yout[i] /  mmGVf - youtX[i] / mmRef ;
              #endif
               cout<<setw(nwdf)<< yout[i] / mmGVf ;
               
               if( printErr ) cout<<setw(nwdf)<< yerr[i] / mmGVf ;
               cout.precision(3);
               
            }

        cout.unsetf (std::ios_base::fixed);        
        cout.setf (std::ios_base::scientific);
        for(int i=3; i<6;i++)
            if(columns[i]){
               if( printSep ) cout << " | " ;  // Separator
               if( printInp ) cout << setw(nwdf-1)<< yin[i] / ppGVf << " ";
              #ifdef BASELINESTEPPER
               if( printInpX ) cout << setw(nwdf-1)<< yinX[i] / ppRef << " ";
               if( printRef ) cout<<setw(nwdf)<< youtX[i] / ppRef; // Reference Solution
               if( printDiff ) 
                  cout<<setw(nwdf+2)<< ( yout[i] /  ppGVf )
                                     - ( youtX[i] / ppRef );
              #endif
               cout<<setw(nwdf) << yout[i] / ppGVf ;
               if( printErr ) cout<<setw(nwdf)<< yerr[i] / ppGVf ;

            }
        cout.unsetf (std::ios_base::scientific);
        
        for(int i=0; i<6;i++)   // Print auxiliary components
        {
           double unitGVf=1;  
           // double unitRef=1;
           // if( i < 3 )             // length / length

           if( i >= 3 ){
              unitGVf = ppGVf / mmGVf; //  dp / ds
            #ifdef BASELINESTEPPER
              unitRef = ppRef / mmRef; //  dp / ds
            #endif
           }
           
           if( i == 0 ){            // length / length                      
              // cout.unsetf (ios_base::scientific);
              cout.setf (std::ios_base::fixed);
           }else if( i == 3 ){
              cout.unsetf (std::ios_base::fixed);              
              cout.setf (std::ios_base::scientific);
           }
           
           if(columns[i+6])
           {
               // cout << " dy/dx["<<i<<"] = ";
               if( printSep ) cout << " | " ;  // Separator
               cout<<setw(nwdf)<< dydx[i] / unitGVf ;
               #ifdef BASELINESTEPPER
               if( printRef )
                 cout<<setw(nwdf)<< dydxRef[i] / unitRef; // Reference Solution
                 if( printDiff ) // Deriv ) 
                  cout<<setw(nwdf)<< ( dydx[i] / unitGVf )
                                   - ( dydxRef[i] / unitRef );
                #endif 
               // bool printDiffDeriv = true;

           }
           // if( i == 2 )     { cout.unsetf(std::ios_base::fixed);      }
           // else if ( i==5 ) { cout.unsetf(std::ios_base::scientific); }
        }
        cout.unsetf(std::ios_base::scientific);
        if( j > 0 )  // Step 0 did not move -- printed the starting values
        {
           // cout.unsetf(std::ios_base::scientific);
           cout.setf (std::ios_base::fixed);                         
           cout.precision(2);
           cout<<setw(nwdf) << atan2(yout[1],yout[0])/degree;
           
          #ifdef BASELINESTEPPER                   // atan(yout[1]/yout[0])/degree;
           if( printRef ) 
             cout<<setw(nwdf) << atan2(youtX[1],youtX[0])/degree;
           #endif

           // Prepare the state of yIn(Vec)
           //Copy yout into yIn
           for(unsigned int i=0; i<Nposmom;i++){
              yInVec[i]  = yOutVec[i];
              #ifdef BASELINESTEPPER
              yInVecX[i] = yOutVecX[i];
              #endif
           }
        }
        
        cout<<"\n";
    }
    
    /*-----------------END-STEPPING------------------*/

    if( debug ) cout<<"----Stepping done "<<endl;

    /*------ Clean up ------*/
    #ifdef BASELINESTEPPER
    exactStepper->InformDone();
    #endif

    return true; // Should perform some checks -- for not just PASS
}

void printBanner()
{
    /*---------------------------------------------------*/
    //        -> First Print the (commented) title header
    cout<<"\n#";
    cout<<setw(5)<<"Step";
    for (int i=0; i<6;i++)
        if (columns[i])
        {
           if( printSep ) cout << " | " ;  // Separator
           if( printInp )
              cout << setw(nwdf-2)<< "yIn[" << i << "]";
           if( printInp )
              cout << setw(nwdf-2)<< "yInX[" << i << "]";           
           cout << setw(nwdf-2)<< "yOut[" << i << "]";
           if( printRef )
              cout << setw(nwdf-2) << "yOutX[" << i << "]";
           if( printErr )
              cout << setw(nwdf-1) << "yErr[" << i << "]";
           if( printDiff )
              cout << setw(nwdf) << " yOut-yOutX[" << i << "]";
        }
    for (int i=0; i<6;i++)
        if (columns[i+6])
        {
           if( printSep ) cout << " | " ;  // Separator
           if( printInp )                               
              cout << setw(nwdf-2)<< "yIn[" << i << "]";           
           cout << setw(nwdf-2)<< "dydx[" << i << "]";
           if( printRef )                    
              cout << setw(nwdf-2) << "dydxRef[" << i << "]";
           // if( printErr ) cout << setw(nwdf-2) << "yErr[" << i << "]";
           if( printDiff ) //  printDiffDeriv )
              cout << setw(nwdf-2) << "dydx-Ref[" << i << "]";
        }    
    cout<<setw(nwdf)<<"tan-1(y/x)";
    cout<<"\n";
}
