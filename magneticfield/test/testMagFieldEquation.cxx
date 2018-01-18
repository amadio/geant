//
//
#include "base/Vector3D.h"

#include "VScalarEquationOfMotion.h"

#include "GUVVectorEquationOfMotion.h"
#include "TVectorMagFieldEquation.h"

#include "VVectorField.h"
#include "MagFieldEquation.h"
#include "FieldEquationFactory.h"

#include "UniformMagField.h"

//#define CMS_FIELD 1

#ifdef CMS_FIELD
// #include "VecMagFieldRoutine/CMSmagField.h"
#include "CMSmagField.h"
#endif


using Double_v = Geant::Double_v;
using Float_v = Geant::Float_v;

template <typename T>
using Vector3D = vecgeom::Vector3D<T>;

using std::cout;
using std::cerr;
using std::endl;
using namespace vecCore::math;

MagFieldEquation<UniformMagField>* 
CreateUniformFieldAndEquation(Vector3D<float> const &fieldValue);

VScalarEquationOfMotion* CreateFieldAndEquation(const char* filename);

template <typename Real_v, typename Equation_t>
bool  TestEquation(Equation_t *);

constexpr unsigned int gNposmom= 6; // Position 3-vec + Momentum 3-vec

const char *defaultFieldFileName= "cmsmagfield2015.txt";

int
main( int, char** )
{
  Vector3D<float>  fieldValue(0.0, 0.0, 1.0);

  using EquationConstField_t = MagFieldEquation<UniformMagField>;
  
  auto eq = CreateUniformFieldAndEquation(fieldValue);

  cout << " Testing scalar.      "; // << endl;
  bool okUniformScalar = TestEquation<double, EquationConstField_t>(eq);
  cout << " Testing Vec Float.   "; // << endl;  
  bool okUniformVecFloat = TestEquation<Float_v, EquationConstField_t>(eq);
  cout << " Testing Vec Double . "; // << endl;    
  bool okUniformVecDouble = TestEquation<Double_v, EquationConstField_t>(eq);

  bool good = okUniformScalar && okUniformVecFloat && okUniformVecDouble;
  
#ifdef CMS_FIELD
  VScalarEquationOfMotion* eq2 = CreateFieldAndEquation( defaultFieldFileName ); // ("cmsMagneticField2015.txt");
  bool okCMSfield = TestEquation(eq2);

  good = good && okCMSfield;
#endif
  
  return good;
}

MagFieldEquation<UniformMagField>* 
CreateUniformFieldAndEquation(Vector3D<float> const &fieldValue)
{
  //  1. Simple method of creating equation
  UniformMagField*   pConstBfield = new UniformMagField( fieldValue );
  return new MagFieldEquation<UniformMagField>(pConstBfield);
}

#ifdef CMS_FIELD
VScalarEquationOfMotion* CreateFieldAndEquation(const char* filename)
{
  //  2. Equation for CMS field
  auto cmsField = new CMSmagField( filename ? filename : defaultFieldFileName ); 
  auto equationCMS = FieldEquationFactory::CreateMagEquation<CMSmagField>(cmsField);

  return equationCMS;
}
#endif

int gVerbose = 1;

template <typename Real_v, typename Equation_t>
bool TestEquation(Equation_t *equation)
{
  const Real_v perMillion = Real_v(1e-6);
  bool   hasError = false;  // Return value
  using Bool_v = vecCore::Mask<Real_v>;
  
  Vector3D<Real_v> PositionVec( 1., 2., 3.);  // initial
  Vector3D<Real_v> MomentumVec( 0., 0.1, 1.);
  Vector3D<Real_v> FieldVec( 0., 0., 1.);  // Magnetic field value (constant)

  // double PositionTime[4] = { PositionVec.x(), PositionVec.y(), PositionVec.z(), 0.0};
  Real_v PositionMomentum[gNposmom];
  Real_v dydx[gNposmom];
  int chargeVec[8] = { -1, 1, 2, -2, -3, -3, -3, -3 };

  // Revise the values, so that they are no longer equal
  Real_v charge = Real_v(0.0); 
  for (size_t lane = 0; lane < vecCore::VectorSize<Real_v>(); ++lane)
    vecCore::Set( charge, lane, chargeVec[lane % 8] );

  PositionMomentum[0] = PositionVec[0];
  PositionMomentum[1] = PositionVec[1];
  PositionMomentum[2] = PositionVec[2];
  PositionMomentum[3] = MomentumVec[0];
  PositionMomentum[4] = MomentumVec[1];
  PositionMomentum[5] = MomentumVec[2];

  equation->EvaluateRhsGivenB( PositionMomentum, FieldVec, charge, dydx );

  Vector3D<Real_v>  ForceVec( dydx[3], dydx[4], dydx[5]);

  // Check result
  Real_v MdotF = MomentumVec.Dot(ForceVec);
  Real_v BdotF = FieldVec.Dot(ForceVec);

  Real_v momentumMag = MomentumVec.Mag();
  Real_v fieldMag =    FieldVec.Mag();
  Real_v sineAngle =   FieldVec.Cross( MomentumVec ).Mag() / ( momentumMag  * fieldMag );

  Real_v ForceMag =   ForceVec.Mag();
  const Real_v c = Real_v(Constants::c_light);

  // Tolerance of difference in values (used below)
  Real_v tolerance = perMillion;
  
  if ( gVerbose ) { std::cout << "Test output:  "  << std::endl; }
  Bool_v error = (Abs(ForceMag - c * Abs(charge) * fieldMag * sineAngle)) > (tolerance * ForceMag);
  if( !vecCore::MaskEmpty(error) ) {
     cerr << "ERROR: Force magnitude is not equal to   c * |charge| * |field| * sin( p, B )."  << endl;     
     cerr << "       Force magnitude = " << ForceMag << endl;
     cerr << "         other side =    " <<  c * Abs(charge) * fieldMag * sineAngle ; 
     cerr << " charge = " << charge 
               << " field-Mag= " << fieldMag  << std::endl;     
     cerr << "       Force = " << ForceVec[0] << " " << ForceVec[1] << " " << ForceVec[2] << " "  << endl;
     hasError = true;
  }
     
  error = (ForceMag == momentumMag * fieldMag);
  assert( vecCore::MaskEmpty(error) );  // Must add coefficient !!

  error = Abs(MdotF) > tolerance * MomentumVec.Mag() * ForceVec.Mag();
  if( !vecCore::MaskEmpty(error) )
  { 
     cerr << "ERROR: Force due to magnetic field is not perpendicular to momentum!!"  << endl;
     hasError= true;
  }
  else if ( gVerbose )
  {
     cout << " Success:  Good (near zero) dot product momentum . force " << endl;
  }
  
  error = Abs(BdotF) > tolerance * FieldVec.Mag() * ForceVec.Mag();
  if( !vecCore::MaskEmpty(error) )
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
