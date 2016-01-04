//
//
#include "base/Vector3D.h"

#include "GUVEquationOfMotion.h"

#include "GUVVectorEquationOfMotion.h"
#include "GUVVectorField.h"
#include "TVectorMagFieldEquation.h"
#include "GUVMagneticField.h"
#include "GUVVectorMagneticField.h"
#include "TVectorUniformMagField.h"

#include "GUVField.h"
#include "TMagFieldEquation.h"
#include "FieldEquationFactory.h"

#include "TUniformMagField.h"

#define CMS_FIELD 1

#ifdef CMS_FIELD
// #include "VecMagFieldRoutine/CMSmagField.h"
#include "CMSmagField.h"
#endif


using ThreeVector_f = vecgeom::Vector3D<float>;
using ThreeVector_d = vecgeom::Vector3D<double>;

using std::cout;
using std::cerr;
using std::endl;

GUVEquationOfMotion* CreateUniformFieldAndEquation(ThreeVector_f constFieldValue);
GUVEquationOfMotion* CreateFieldAndEquation(const char* filename);

bool  TestEquation(GUVEquationOfMotion* );

constexpr unsigned int gNposmom= 6; // Position 3-vec + Momentum 3-vec

const char *defaultFieldFileName= "cmsmagfield2015.txt";

int
main( int, char** )
{
  ThreeVector_f  FieldValue(0.0, 0.0, 1.0);
   
  GUVEquationOfMotion* eq = CreateUniformFieldAndEquation( FieldValue );
  TestEquation(eq);

#ifdef CMS_FIELD
  GUVEquationOfMotion* eq2 = CreateFieldAndEquation( defaultFieldFileName ); // ("cmsMagneticField2015.txt");
  TestEquation(eq2);
#endif
  
  return 1;
}

GUVEquationOfMotion* CreateUniformFieldAndEquation(ThreeVector_f FieldValue)
{
  // using Field_t = TUniformMagField;

  TUniformMagField*   pConstBfield = new TUniformMagField( FieldValue );

  // 1. Original way of creating an equation
  // using EquationType = TMagFieldEquation<TUniformMagField, gNposmom>;
  // GUVEquationOfMotion*  magEquation = new EquationType(pConstBfield);
  // return magEquation;

  //  2. Different method of creating equation:  Factory
  return FieldEquationFactory::CreateMagEquation<TUniformMagField>(pConstBfield);
}

#ifdef CMS_FIELD
GUVEquationOfMotion* CreateFieldAndEquation(const char* filename)
{
  // const char *defaultFieldFileName= "cmsMagneticField2015.txt";
   
  //  3. Equation for CMS field
  auto cmsField = new CMSmagField( filename ? filename : defaultFieldFileName ); 
  auto equationCMS = FieldEquationFactory::CreateMagEquation<CMSmagField>(cmsField);

  return equationCMS;
}
#endif

int gVerbose = 1;

bool TestEquation(GUVEquationOfMotion *equation)
{
  constexpr double perMillion = 1e-6;
  bool   hasError = false;  // Return value
  
  ThreeVector_d PositionVec( 1., 2.,  3.);  // initial
  ThreeVector_d MomentumVec( 0., 0.1, 1.);
  ThreeVector_f FieldVec( 0., 0., 1.);  // Magnetic field value (constant)

  // double PositionTime[4] = { PositionVec.x(), PositionVec.y(), PositionVec.z(), 0.0};
  double PositionMomentum[gNposmom];
  double dydx[gNposmom];
  double charge= -1;

  PositionMomentum[0] = PositionVec[0];
  PositionMomentum[1] = PositionVec[1];
  PositionMomentum[2] = PositionVec[2];
  PositionMomentum[3] = MomentumVec[0];
  PositionMomentum[4] = MomentumVec[1];
  PositionMomentum[5] = MomentumVec[2];

  // double FieldArr[3]= { FieldVec.x(), FieldVec.y(), FieldVec.z() };
  
  equation->InitializeCharge( charge );
  equation->EvaluateRhsGivenB( PositionMomentum, FieldVec, /* charge, */ dydx );

  ThreeVector_d  ForceVec( dydx[3], dydx[4], dydx[5]);

  // Check result
  double MdotF = MomentumVec.Dot(ForceVec);
  double BdotF = FieldVec.Dot(ForceVec);

  double momentumMag = MomentumVec.Mag();
  double fieldMag =    FieldVec.Mag();
  double sineAngle =   FieldVec.Cross( MomentumVec ).Mag() / ( momentumMag  * fieldMag );

  double ForceMag =   ForceVec.Mag();
  const double c = Constants::c_light;

  std::cout << "Test output:  "  << std::endl;     
  if( ForceMag != c * std::fabs(charge) * fieldMag * sineAngle ) {
     cerr << "ERROR: Force magnitude is not equal to   c * |charge| * |field| * sin( p, B )."  << endl;     
     cerr << "       Force magnitude = " << ForceMag << endl;
     cerr << "         other side =    " <<  c * std::fabs(charge) * fieldMag * sineAngle ; 
     cerr << " charge = " << charge 
               << " field-Mag= " << fieldMag  << std::endl;     
     cerr << "       Force = " << ForceVec[0] << " " << ForceVec[1] << " " << ForceVec[2] << " "  << endl;
  }
     
  assert( ForceMag != momentumMag * fieldMag );  // Must add coefficient !!
  
  if( std::fabs(MdotF) > perMillion * MomentumVec.Mag() * ForceVec.Mag() )
  { 
     cerr << "ERROR: Force due to magnetic field is not perpendicular to momentum!!"  << endl;
     hasError= true;
  }
  else if ( gVerbose )
  {
     cout << " Success:  Good (near zero) dot product momentum . force " << endl;
  }
  if( std::fabs(BdotF) > perMillion * FieldVec.Mag() * ForceVec.Mag() )
  { 
    cerr << "ERROR: Force due to magnetic field is not perpendicular to B field!"
              << std::endl; 
    cerr << " Vectors:  BField   Force " << std::endl;
    for ( int i = 0; i < 3; i ++ )
       cerr << "   [" << i << "] " << FieldVec[i] << " " << ForceVec[i] << std::endl;

    hasError = true;
  }
  else if ( gVerbose )
  {
    cout << " Success:  Good (near zero) dot product magnetic-field . force " << std::endl;
  }

  return hasError;
}
