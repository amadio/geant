//
//
#include "GUVEquationOfMotion.h"

#include "TMagFieldEquation.h"
#include "TUniformMagField.h"

GUVEquationOfMotion *CreateFieldAndEquation(ThreeVector constFieldValue);
bool TestEquation(GUVEquationOfMotion *);

const unsigned int Nposmom = 6; // Position 3-vec + Momentum 3-vec

ThreeVector FieldValue(0.0, 0.0, 1.0);

int main(int, char **)
{
  GUVEquationOfMotion *eq = CreateFieldAndEquation(FieldValue);
  TestEquation(eq);

  return 1;
}

// GUVField*       pField;

GUVEquationOfMotion *CreateFieldAndEquation(ThreeVector FieldValue)
{
  TUniformMagField *ConstBfield = new TUniformMagField(FieldValue);
  // typedef typename TMagFieldEquation<TUniformMagField, Nposmom>  EquationType;
  using EquationType = TMagFieldEquation<TUniformMagField, Nposmom>;

  // GUVEquationOfMotion*  magEquation= new TMagFieldEquation<TUniformMagField, Nposmom>(ConstBfield);
  GUVEquationOfMotion *magEquation = new EquationType(ConstBfield);

  // pField = ConstBfield;
  return magEquation;
}

int gVerbose = 1;

bool TestEquation(GUVEquationOfMotion *equation)
{
  const double perMillion = 1e-6;
  bool hasError = false; // Return value

  ThreeVector PositionVec(1., 2., 3.); // initial
  ThreeVector MomentumVec(0., 0.1, 1.);
  ThreeVector FieldVec(0., 0., 1.); // Magnetic field value (constant)

  double PositionTime[4] = {PositionVec.x(), PositionVec.y(), PositionVec.z(), 0.0};
  // double PositionTime[4] = { 1., 2., 3., 4.};
  PositionTime[0] = PositionVec.x();
  PositionTime[1] = PositionVec.y();
  PositionTime[2] = PositionVec.z();

  // double magField[3];

  double dydx[Nposmom];
  //  double PositionMomentum[Nposmom];

  double charge = -1;

  //  PositionMomentum[0]= PositionVec[0];
  //  PositionMomentum[1]= PositionVec[1];
  //  PositionMomentum[2]= PositionVec[2];
  // PositionMomentum[3]= MomentumVec[0];
  // PositionMomentum[4]= MomentumVec[1];
  // PositionMomentum[5]= MomentumVec[2];

  double FieldArr[3] = {FieldVec.x(), FieldVec.y(), FieldVec.z()};

  equation->InitializeCharge(charge);
  equation->EvaluateRhsGivenB(PositionTime, FieldArr, /* charge, */ dydx);

  ThreeVector ForceVec(dydx[3], dydx[4], dydx[5]);

  // Check result
  double MdotF = MomentumVec.Dot(ForceVec);
  double BdotF = FieldVec.Dot(ForceVec);

  double momentumMag = MomentumVec.Mag();
  double fieldMag = FieldVec.Mag();
  double ForceMag = ForceVec.Mag();

  if (ForceMag != momentumMag * fieldMag) {
    std::cerr << "ERROR: Force magnitude is not equal to momentum * field." << std::endl;
  }

  assert(ForceMag != momentumMag * fieldMag); // Must add coefficient !!

  if (std::fabs(MdotF) > perMillion * MomentumVec.Mag() * ForceVec.Mag()) {
    std::cerr << "ERROR: Force due to magnetic field is not perpendicular to momentum!!" << std::endl;
    hasError = true;
  }
  else if (gVerbose) {
    std::cout << " Success:  Good (near zero) dot product momentum . force " << std::endl;
  }
  if (std::fabs(BdotF) > perMillion * FieldVec.Mag() * ForceVec.Mag()) {
    std::cerr << "ERROR: Force due to magnetic field is not perpendicular to B field!" << std::endl;
    std::cerr << " Vectors:  BField   Force " << std::endl;
    for (int i = 0; i < 3; i++)
      std::cerr << "   [" << i << "] " << FieldVec[i] << " " << ForceVec[i] << std::endl;

    hasError = true;
  }
  else if (gVerbose) {
    std::cout << " Success:  Good (near zero) dot product magnetic-field . force " << std::endl;
  }

  return hasError;
}
