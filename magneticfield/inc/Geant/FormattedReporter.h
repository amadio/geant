#include <fstream>
#include <cstdlib>
#include <iomanip>

// #include "Geant/IntegrationDriverConstants.h"  // For ReportOneLane

//  Auxiliary methods - should be encapsulated into a separate helper class

#ifndef FORMATTED_REPORTER_H
#define FORMATTED_REPORTER_H

namespace FormattedReporter // Was ReportValuesOfVectors
{

const int sDefaultNameLength = 14;
const int sDefaultPrecision  = 9;
const int sDefaultVarSize    = sDefaultPrecision + 7;  // Was 12

template <typename Real_v>
void ReportRowOfDoubles(std::string varName, const Real_v &varValue, int charName = -1, int widthVal = -1)
{
  using std::cout;
  using std::endl;

  const int minPrec= 9;
  
  // int argWV= widthVal;
  if (charName <= 0) {
    charName = sDefaultNameLength;
  }
  if (widthVal <= 0) {
    widthVal = sDefaultVarSize;
  }

  int prec = (widthVal - 7);
  prec     = std::min(minPrec, prec);
  widthVal = prec + 7;

  cout << std::setw(charName) << varName << " : ";
  for (size_t i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
    cout << " " << std::setw(widthVal) << std::setprecision(prec) << vecCore::Get(varValue, i) << " | ";
  }
  // Auxiliary information about width, precision
  // cout << " withVal: arg= " << argWV << " used= " << widthVal << " ( prec= " << prec << " ) ";
  // cout << " lenName = " << charName;
  cout << endl;
}

// ---------------------------------------------
template <typename Real_v>
void ReportRowOfSquareRoots(std::string varName,
                            const Real_v &valueSq, // Square of interesting value
                            int charName = -1, int widthVal = -1)
{
  if (vecCore::MaskEmpty(valueSq < 0.0)) {
    Real_v value = vecCore::math::Sqrt(valueSq);
    ReportRowOfDoubles(varName, value, charName, widthVal);
  } else {
    // There is an erroneous value !
    std::string varNameAndWarning = "WARNING> some values are Negative> " + varName;
    ReportRowOfDoubles(varNameAndWarning, valueSq, charName, widthVal);
  }
}
// ---------------------------------------------

template <typename Real_v>
void ReportManyRowsOfDoubles(std::string varArrName, const Real_v varArr[], int arrLen, int charName = -1,
                             int widthVal = -1)
{
  for (int i = 0; i < arrLen; i++) {
    // std::ostringstream nameAndIndex;
    // nameAndIndex << varArrName << "[" << i << "]/AF";
    // ReportRowOfDoubles( nameAndIndex.str(), yValues[i] );
    std::string nameAndIndex = varArrName + "[" + std::to_string(i) + "]"; // + "/AF";
    ReportRowOfDoubles<Real_v>(nameAndIndex, varArr[i], charName, widthVal);
  }
  std::cout << "##-------------------------------------------------------------------------------" << std::endl;
}

// -----------------------------------------------------------
template <typename Real_v>
Real_v GetMomentumMag(const Real_v varPositionsMomenta[6])
{
  Real_v px = varPositionsMomenta[3];
  Real_v py = varPositionsMomenta[4];
  Real_v pz = varPositionsMomenta[5];

  return vecCore::math::Sqrt(px * px + py * py + pz * pz);
}

// -----------------------------------------------------------

template <typename Real_v>
void ReportRowsOfPositionsMomenta(std::string varName, const Real_v varPositionsMomenta[], int arrLen,
                                  const Real_v &momentumMagStart, //
                                  int widthNm = -1, int widthVal = -1)
{
  using vecCore::MaskEmpty;
  using vecCore::MaskFull;
  using vecCore::math::Sqrt;
  if (widthVal < 0) {
    widthVal = sDefaultVarSize;
  }

  ReportManyRowsOfDoubles(varName, varPositionsMomenta, arrLen, widthNm, widthVal);
  assert(arrLen >= 6);

  Real_v momEnd = GetMomentumMag(varPositionsMomenta);

  Real_v diffMagP = momEnd - momentumMagStart;
  ReportRowOfDoubles("diff|p|", diffMagP);

  if (!MaskFull(momentumMagStart == Real_v(0.0))) {
    double tinyVal = 1.0e-80;
    Real_v relDiff = diffMagP / (momentumMagStart + Real_v(tinyVal));
    ReportRowOfDoubles("d|p|/|p|", relDiff);

    double thresholdRelativeDiff = 1.0e-5; //  Later: 3 * epsilon ??
    if (!MaskEmpty(vecCore::math::Abs(relDiff) > Real_v(thresholdRelativeDiff))) {
      int extraWidth = widthVal + 12;
      ReportRowOfDoubles("|momEnd|", momEnd, widthNm, extraWidth);
      ReportRowOfDoubles("|momStart|", momentumMagStart, widthNm, extraWidth);
    } else {
      ReportRowOfDoubles("|momEnd|", momEnd);
    }
    std::cout << "##-------------------------------------------------------------------------------" << std::endl;
  }
}

// ---------------------------------------------

template <typename Real_v>
inline void ReportRowOfBools(std::string varName, const vecCore::Mask_v<Real_v> &var, int widthName = -1,
                             int widthVal = -1)
{
  using std::cout;

  if (widthName < 0) {
    widthName = sDefaultNameLength;
  }
  if (widthVal < 0) {
    widthVal = sDefaultVarSize;
  }

  cout << std::setw(widthName) << varName << " : ";
  for (size_t i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
    cout << " " << std::setw(widthVal) << vecCore::Get(var, i) << " | ";
  }
  cout << std::endl;
}

// ---------------------------------------------

// ===============  Selective Reporting / Printing ==================

template <typename Real_v>
inline void ReportRowOfDoublesIf(std::string varName, const Real_v var, vecCore::Mask_v<Real_v> cond,
                                 int widthName = -1, int widthVal = -1)
{
  using std::cout;
  if (widthName < 0) {
    widthName = sDefaultNameLength;
  }
  if (widthVal < 0) {
    widthVal = sDefaultVarSize;
  }

  cout << std::setw(widthName) << varName << " : ";
  for (unsigned int i = 0; i < vecCore::VectorSize<Real_v>(); ++i)
  {
    if (vecCore::Get(cond, i))
       cout << " " << std::setw(widthVal) << vecCore::Get(var, i) << " | ";
    else
       cout << " " << std::setw(widthVal) << "-/NA "     << " | ";

    // if( i+1 << vecCore::VectorSize<Real_v>() ) { cout << " | "; }
       
  }
  cout << std::endl;
}

// ----------------------------------------------------------------------------------
inline void ReportArray(const char *methodName, const std::string &variableName, const double Arr[], int numTracks,
                        bool banner = false)
{
  using std::cout;
  using std::endl;

  const int precisionVal = 4;
  const int wdName       = 12;
  const int charWidth    = precisionVal + 2;

  if (banner) {

    cout << " **** Method " << std::setw(wdName) << methodName << " values of arrays: " << endl;
    cout << std::setw(wdName) << "Variable Name"
         << " :";
    for (int i = 0; i < numTracks; ++i) {
      cout << " [" << std::setw(charWidth - 3) << i << "] ";
    }
    cout << endl;
  }
  cout << std::setw(wdName) << variableName << " : ";
  int oldPrec = cout.precision(precisionVal);
  for (int i = 0; i < numTracks; ++i) {
    // cout << " [" << i << "]= ";
    cout << std::setw(charWidth) << Arr[i] << " ";
  }
  cout << std::endl;
  cout.precision(oldPrec);
}

// ----------------------------------------------------------------------------------

template<typename Double_v, typename Bool_v>
inline
void
   ReportOneLane( Double_v hStep,
                  Double_v xStepStart,
                  Double_v epsPosition, Double_v errPosSq,
                  Double_v errMomSq,    Double_v errmax_sq,
                  Bool_v   lanesDone,   int     allDone,
                  int      iter,        int      noCall,
                  int      lanePrint,   int      trackNum,
                  const char *methodName)
{
   using std::cout;
   using std::setw;   
   bool  laneIsDone = vecCore::Get( lanesDone , lanePrint );
   int   prec = 10; // precision
   int   wd = prec + 5;
   int   oldPrec= cout.precision(prec);
   bool  printSquares = false; // Old version - true
   bool  printValues  = true; 

   // const int trackToPrint = IntegrationDriverConstants::GetInstance()->GetTrackToCheck();   
   
   cout << std::setw(12) << methodName << " - ReportOneLane : "
        << " trk# "     << setw(3) << trackNum  /* trackToPrint */ << " "
        << " lane: "    << setw(3) << lanePrint << " > "
        << " iter = "   << setw(3) << iter << " #call= " << setw(5) << noCall;
   prec=6;
   wd = prec + 5;
   cout << std::setprecision( prec )
        << " h = "          << setw( wd ) << vecCore::Get( hStep ,       lanePrint )
        << " xStepStart = " << setw( wd ) << vecCore::Get( xStepStart ,  lanePrint )      
        << " Eps/x = "      << setw( wd ) << vecCore::Get( epsPosition , lanePrint );

   double errPosLane2 = vecCore::Get( errPosSq ,    lanePrint );
   double errMomLane2 = vecCore::Get( errMomSq ,    lanePrint );
   double errMax2     = vecCore::Get( errmax_sq ,   lanePrint );

   int   prec2 = 9; // precision
   int   wd2 = prec2 + 5;

   cout.precision( prec2 );
            
   if( printSquares ) 
      cout
             << " errSq: x,p = "  << setw( wd2) << errPosLane2  // vecCore::Get( errPosSq ,    lanePrint )
             << " "               << setw( wd2) << errMomLane2  // vecCore::Get( errMomSq ,    lanePrint ) // << " "
             << " errMax^2 = "    << setw( wd2) << errMax2;     // vecCore::Get( errmax_sq ,   lanePrint );
   if( printValues ) 
      cout
             << " error-x/p = "   << setw( wd2) << sqrt( errPosLane2 ) // vecCore::Get( errPosSq ,    lanePrint ) )
             << " "               << setw( wd2) << sqrt( errMomLane2 ) // vecCore::Get( errMomSq ,    lanePrint ) ) // << " "
             << " errMax = "      << setw( wd2) << sqrt( errMax2 )  ; // vecCore::Get( errmax_sq ,   lanePrint ) );

   cout << " lane done = "   << laneIsDone;

   if( allDone >= 0 )
      cout << " allDone = "     << allDone;
             
   cout << std::endl;

   if( laneIsDone  ) cout << std::endl;

   cout.precision(oldPrec);
   
   // if( allDone ) cout << std::endl;
   // **> Old mode - that printed all updates - not just ones in which this lane was active.
}

template <class Real_v> // , unsigned int Nvar>
void
FullReport(const Real_v yStepStart[],
           const Real_v charge,
           const Real_v dydx[],
           const Real_v hStep,
           const Real_v yStepEnd[],
           const Real_v yEstErr[],
           const Real_v errmax_sq,
           const vecCore::Mask_v<Real_v> Active
   )
{
  // using FormattedReporter::ReportRowOfDoubles;
  // using FormattedReporter::ReportManyRowsOfDoubles;
  // using FormattedReporter::ReportRowOfBools;
    
  // std::cout << "RID: Status after stepper call ----------------------------------------------" << std::endl;
  std::cout << "Status   -----------------------------------------------------------------------------------" << std::endl;  
  std::cout << "-- Argument values" << std::endl;
  ReportRowOfDoubles("hStep(ask)", hStep );
  std::cout << "-- Start  values" << std::endl;
  ReportManyRowsOfDoubles("StartX/P", yStepStart, 6 );
  ReportRowOfDoubles("Charge",  charge, 6 );              
  ReportManyRowsOfDoubles("d[X,P]/ds", dydx, 6 );

  std::cout << "-- Return values" << std::endl;
  ReportManyRowsOfDoubles("X,P_out", yStepEnd, 6 );

  // Report Estimated Errors       
  std::cout << "-- Estimated Errors" << std::endl;
  ReportManyRowsOfDoubles("errXP/xyz", yEstErr, 6 );
  
  // ReportManyRowsOfDoubles("err-p/xyz", &yerr[3], 3 );
  
  // ReportRowOfSquareRoots("|err-p|", yerr[3]*yerr[3] + yerr[4]*yerr[4] + yerr[5]*yerr[5] );       
  // ReportRowOfDoubles("up = SumErr^2", sumerr_sq );
  // ReportRowOfDoubles("dwn= magMom^2+e", magmom_sq + tinyValue );
  // ReportRowOfDoubles("mul:1/e_vel^2", invEpsilonRelSq );
  // ReportRowOfDoubles("ErrMom", errmom_sq );
  // ReportRowOfSquareRoots("ErrMom", errmom_sq );       
  ReportRowOfDoubles("ErrMaxSq", errmax_sq );
  ReportRowOfSquareRoots("ErrMax", errmax_sq );       
  ReportRowOfBools<Real_v>("Active(old)", Active);
  std::cout << "End Status -------------------------------------------------------------------------------------" << std::endl;
}



}; // End of namespace FormattedReporter

#endif // FORMATTED_REPORTER_H
