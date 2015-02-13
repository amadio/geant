//
//
#include "TMagFieldEquation.h"
#include "TUniformMagField.h"

GUVEquationOfMotion*  CreateFieldAndEquation();
bool  TestEquation(TGUVEquationOfMotion* );

const unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

main()
{
  CreateFieldAndEquation();
  TestEquation();  
}

ThreeVector    FieldValue(0.0, 0.0, 1.0);
GUField*       pField;

GUVEquationOfMotion* CreateFieldAndEquation(ThreeVector FieldValue)
{
  TUniformMagField*     ConstBfield= new TUniformMagField( FieldValue );

  pField = ConstBField; 
  GUVEquationOfMotion*  magEquation= TMagFieldEquation<TUniformMagField, Nposmom>;
}

G4int gVerbose= 1;

bool TestEquation(GUVEquationOfMotion* equation)
{
  ThreeVector PositionVec( 1., 2., 3.);
  ThreeVector MomentumVec( 0., 0., 1.);

  double PositionTime[4] = { PositionVec.x(), PositionVec.y(), PositionVec.z(), 0.0};
  // double PositionTime[4] = { 1., 2., 3., 4.};
  PositionTime[0]= PositionVec.x();
  PositionTime[1]= PositionVec.y();
  PositionTime[2]= PositionVec.z();

  double magField[3];

  double dydx[Nposmem]
  double PositionMomentum[Nposmom];

  double charge= -1;

  PositionMomentum[0]= PositionVec[0];
  PositionMomentum[1]= PositionVec[1];  
  PositionMomentum[2]= PositionVec[2];
  PositionMomentum[3]= MomentumVec[0];
  PositionMomentum[4]= MomentumVec[1];  
  PositionMomentum[5]= MomentumVec[2];

  FieldArr[3]= { FieldVec.x(), FieldVec.y(), FieldVec.z() };

  equation->EvaluateRHSGivenB( PositionTime, FieldArr, charge, dydx );

  ThreeVector  Force( dydx[3], dydx[4], dydx[5]);

  // Check result
  MdotF= MomentumVec.dot(Force);
  BdotF= FieldVec.dot(Force); 

  if( std::fabs(MdotF) > perMillion * Momentum.Mag() * Force.Mag() )
  { 
    std::cerr << "ERROR: Force due to magnetic field is not perpendicular to momentum!!" 
  }
  else if ( gVerbose )
  {
    std::cout << " Success:  Good (near zero) dot product momentum . force " << G4endl;
  }
  if( std::fabs(BdotF) > perMillion * Field.Mag() * Force.Mag() )
  { 
    std::cerr << "ERROR: Force due to magnetic field is not perpendicular to B field!" 
    std::cerr << " Vectors:  BField   Force " << G4endl;
    for ( i=0; i<3; i++ )
      std::cerr << "   [" << setwidth() << i << "] << Bfield(i) << " " << Force(i) << G4endl;    
  }
  else if ( gVerbose )
  {
    std::cout << " Success:  Good (near zero) dot product magnetic-field . force " << G4endl;
  }
}