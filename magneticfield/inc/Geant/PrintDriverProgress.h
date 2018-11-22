
#include "Geant/FormattedReporter.h"

#ifndef PRINTDRIVERPROGRESS_H
#define PRINTDRIVERPROGRESS_H

namespace ReportValuesOfVectors
{

// Auxiliary printing methods 

template <class Real_v>
void ReportConditionLanes(vecCore::Mask_v<Real_v> problemLane, Real_v x, Real_v xnew, Real_v h, Real_v htry); // const;

template <class Real_v>
void ReportStatus(const Real_v x, const Real_v charge, const Real_v hTry, const Real_v errmaxSqFinal,
                  const Real_v yValues[]); // const;

template <class Real_v>
void ReportResults(const Real_v hFinal, const Real_v errmaxSqFinal, const Real_v yOutput[]); // const;

// ------------------------------------------------------------------------------------
// Definitions



template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
void // SimpleIntegrationDriver</*Real_v,*/ T_Stepper, Nvar>::
   ReportConditionLanes(vecCore::Mask_v<Real_v> problemVec,
                        Real_v xVec, Real_v xnewVec,
                        Real_v hVec, Real_v htryVec) // const
{
  using std::cerr;
  using std::endl;
  using vecCore::Get;

  for (size_t i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
    if (Get(problemVec, i)) { // problemVec[i]
      double x = Get(xVec, i);
      double h = Get(hVec, i);

      double xnew = Get(xnewVec, i);
      double htry = Get(htryVec, i);

      double xnewCheck = x + h; // x[i] + h[i];
      double diffX     = xnew - x;
      cerr.precision(16);
      cerr << " WARNING (SimpleIntegrationDriver::OneGoodStep()> Underflow in lane " << i << endl
           << "   Step's start and end are equal !  Values: " << endl
           << "      Start x = " << std::setw(20) << x << endl
           << "      End   x = " << std::setw(20) << xnew << " check= " << xnewCheck << endl;
      cerr.precision(6);
      cerr << "      Diff    = " << diffX << endl;
      cerr.precision(9);
      cerr << "     Step-size = " << h << endl << "   Input step = " << htry << endl;
    }
  }
}

// ---------------------------------------------------------


template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
void SimpleIntegrationDriver<T_Stepper, Nvar>::ReportStatus(const Real_v x, const Real_v charge, const Real_v hTry,
                                                            const Real_v errmaxSqFinal, const Real_v yValues[]) const
{
  // Method to report intermediate state of integration, including
  //   - error ratios
  //   - status of finishing, underflow, storing, ..
  //   - integrated variables
  using FormattedReporter::ReportRowOfDoubles;

  ReportRowOfDoubles("x", x);
  ReportRowOfDoubles("charge", charge);
  // std::cout << " Check - charge: " << charge << std::endl;
  ReportRowOfDoubles("hTry", hTry);
  ReportRowOfDoubles("errMaxSq/F", errmaxSqFinal);

  // ReportRowOfDoubles( "yValues- 0/F", yValues[i] );
  std::string baseName = "yNow";
  for (unsigned int i = 0; i < Nvar; i++) {
    std::string varName = baseName + "[" + std::to_string(i) + "]/Rs";
    ReportRowOfDoubles(varName, yValues[i]);
    // ReportRowOfDoubles( "yNow- 0/F", yValues[i] );
  }
}



template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
void SimpleIntegrationDriver<T_Stepper, Nvar>::ReportResults(const Real_v hFinal, const Real_v errmaxSqFinal,
                                                             const Real_v yOutput[]) const
{
  // Method to report intermediate state of integration, including
  //   - error ratios
  //   - status of finishing, underflow, storing, ..
  //   - integrated variables
  using std::cout;
  using std::endl;
  using std::setw;
  using vecCore::Get;

  const int chName = 10, chValue = 10;
  std::cout << setw(chName) << "h"
            << " : ";
  for (size_t i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
    cout << " " << setw(chValue) << Get(hFinal, i) << " | ";
  }
  cout << endl;

  ReportRowOfDoubles("hFinal", hFinal);
  //
  cout << setw(chName) << "errMaxSq/M"
       << " : ";
  for (size_t i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
    cout << " " << setw(chValue) << Get(errmaxSqFinal, i) << " | ";
  }
  cout << endl;

  ReportRowOfDoubles("errMaxSq/F", errmaxSqFinal);

  ReportManyRowsOfDoubles("yOut/R", yOutput, Nvar);

  std::cout << "-- Same value with 'manual' loop (to check): " << std::endl;
  std::string baseName = "yOut";
  for (int i = 0; i < Nvar; i++) {
    std::string varName = baseName + "[" + std::to_string(i) + "]/M";
    ReportRowOfDoubles(varName, yOutput[i]);
    // ReportRowOfDoubles( "yOut- 0/F", yOut[i] );
  }
}

#endif

