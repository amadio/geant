#include <fstream>
#include <cstdlib>
#include <iomanip>

//  Auxiliary methods - should be encapsulated into a separate helper class

#ifndef  FORMATTED_REPORTER_H
#define  FORMATTED_REPORTER_H

namespace FormattedReporter   // Was ReportValuesOfVectors
{

const int sDefaultNameLength= 14;
const int sDefaultVarSize =   12;

// using std::setw;
// using std::cout;
// using std::endl;

template< typename Real_v >
   void
   ReportRowOfDoubles( std::string    varName,
                       const Real_v & varValue, 
                       int            charName = -1,
                       int            widthVal = -1 )
{
   using std::cout;
   using std::endl;

   // int argWV= widthVal;
   if( charName <=  0 ) { charName  = sDefaultNameLength; }
   if( widthVal <= 0 ) { widthVal = sDefaultVarSize; }

   int prec = (widthVal - 7);
   prec = std::min( 5, prec );
   widthVal = prec + 7;
   
   cout << std::setw( charName ) << varName << " : ";
   for (int i = 0; i < vecgeom::kVectorSize; ++i) {   
      cout << " " << std::setw( widthVal )
                << std::setprecision( prec )
                << vecCore::Get( varValue, i ) << " | ";      
   }
   // Auxiliary information about width, precision
   // cout << " withVal: arg= " << argWV << " used= " << widthVal << " ( prec= " << prec << " ) ";
   // cout << " lenName = " << charName; 
   cout << endl;
}

// ---------------------------------------------
template< typename Real_v >
   void
   ReportRowOfSquareRoots( std::string    varName,
                           const Real_v & valueSq,      // Square of interesting value
                           int            charName = -1,
                           int            widthVal = -1 )
{
   if( vecCore::MaskEmpty( valueSq < 0.0)  ) {
      Real_v value= vecCore::math::Sqrt( valueSq );
      ReportRowOfDoubles( varName, value, charName, widthVal );
   } else {
      // There is an erroneous value !
      std::string varNameAndWarning= "WARNING> some values are Negative> " + varName;
      ReportRowOfDoubles( varNameAndWarning, valueSq, charName, widthVal );
   }   
}
// ---------------------------------------------

template< typename Real_v >
   void
   ReportManyRowsOfDoubles( std::string    varArrName,
                            const Real_v   varArr[],
                            int            arrLen, 
                            int            charName = -1,
                            int            widthVal = -1 )
{
   for( int i=0; i< arrLen; i++ ) {
      // std::ostringstream nameAndIndex;
      // nameAndIndex << varArrName << "[" << i << "]/AF";
      // ReportRowOfDoubles( nameAndIndex.str(), yValues[i] );
      std::string nameAndIndex = varArrName + "[" + std::to_string(i) + "]"; // + "/AF";
      ReportRowOfDoubles<Real_v>( nameAndIndex, varArr[i], charName, widthVal );
   }
   std::cout << "##-------------------------------------------------------------------------------" << std::endl;
}

// ---------------------------------------------
template< typename Real_v >
   Real_v GetMomentumMag( const Real_v varPositionsMomenta[6] )
{
   Real_v px = varPositionsMomenta[3];
   Real_v py = varPositionsMomenta[4];
   Real_v pz = varPositionsMomenta[5];   

   return vecCore::math::Sqrt( px * px + py * py + pz * pz );
}

template< typename Real_v >
   void
   ReportRowsOfPositionsMomenta( std::string    varName, 
                                 const Real_v   varPositionsMomenta[],
                                 int            arrLen,
                                 const Real_v & momentumMagStart, // 
                                 int            widthNm  = -1,
                                 int            widthVal = -1 )
{
   using vecCore::math::Sqrt;
   using vecCore::MaskEmpty;
   using vecCore::MaskFull;   
   if( widthVal < 0 ) { widthVal = sDefaultVarSize; }
   
   ReportManyRowsOfDoubles( varName,
                            varPositionsMomenta,
                            arrLen,
                            widthNm, 
                            widthVal );
   assert( arrLen >= 6 ); 
   
   Real_v momEnd= GetMomentumMag( varPositionsMomenta );
   
   Real_v diffMagP = momEnd - momentumMagStart;
   ReportRowOfDoubles( "diff|p|", diffMagP );

   if( ! MaskFull( momentumMagStart == Real_v( 0.0 ) ) )
   {
      double tinyVal = 1.0e-80;
      Real_v relDiff = diffMagP / ( momentumMagStart + Real_v(tinyVal) ) ;
      ReportRowOfDoubles( "d|p|/|p|",   relDiff  );
      
      double thresholdRelativeDiff = 1.0e-5;   //  Later: 3 * epsilon ??         
      if( ! MaskEmpty( abs(relDiff) > thresholdRelativeDiff ) )
      {
         int extraWidth= widthVal + 12;
         ReportRowOfDoubles( "|momEnd|",   momEnd,           widthNm, extraWidth  );
         ReportRowOfDoubles( "|momStart|", momentumMagStart, widthNm, extraWidth  );
      }else{
         ReportRowOfDoubles( "|momEnd|",   momEnd  );
      }
      std::cout << "##-------------------------------------------------------------------------------"
                << std::endl;
   }
}

// ---------------------------------------------


template< typename Real_v >
   void
   ReportRowOfBools( std::string                     varName,
                     const vecCore::Mask_v<Real_v> & var,
                     int          widthName= -1,
                     int          widthVal= -1)
{
   using std::cout;
   
   if( widthName < 0 ) { widthName  =  sDefaultNameLength; }
   if( widthVal < 0 ) { widthVal = sDefaultVarSize; }
   
   cout << std::setw( widthName ) << varName << " : ";
   for (int i = 0; i < vecgeom::kVectorSize; ++i) {
      cout << " " << std::setw( widthVal ) << vecCore::Get( var, i ) << " | ";
   }
   cout << std::endl;
}

// ---------------------------------------------


// ===============  Selective Reporting / Printing ==================

template< typename Real_v >
void
ReportRowOfDoublesIf( std::string    varName,
                      const Real_v   var, 
                      vecCore::Mask_v<Real_v> cond,
                      int            widthName = -1,
                      int            widthVal = -1                 
   )
{
   using std::cout;
   if( widthName < 0 ) { widthName  =  sDefaultNameLength; }
   if( widthVal < 0 )  { widthVal = sDefaultVarSize; }
   
   cout << std::setw( widthName ) << varName << " : ";
   for (int i = 0; i < vecgeom::kVectorSize; ++i) {
      if( Get( cond, i ) )
         cout << " " << std::setw( widthVal ) << Get(var, i) << " | ";      
   }
   cout << std::endl;
}


// ----------------------------------------------------------------------------------

void ReportArray( const std::string& methodName,
                  const std::string& variableName,                  
                  const double       Arr[],
                  int                numTracks,
                  bool               banner= false)
{
   using std::cout;
   using std::endl;

   const int precisionVal = 4;
   const int wdName= 12;
   const int charWidth    = precisionVal + 2;

   if( banner ) {
      
      cout << " **** Method " 
         << std::setw(wdName) << methodName << " values of arrays: " << endl;
      cout << std::setw(wdName) << "Variable Name" << " :";      
      for (int i = 0; i < numTracks; ++i)
      {
         cout << " [" << std::setw(charWidth-3) << i << "] ";
      }
      cout << endl;
   }
   cout << std::setw(wdName) << variableName << " : ";
   int oldPrec= cout.precision( precisionVal );
   for (int i = 0; i < numTracks; ++i)
   {
      // cout << " [" << i << "]= ";
      cout << std::setw(charWidth) << Arr[i] << " ";
   }
   cout<< std::endl;
   cout.precision( oldPrec ); 
}

};  // End of namespace

#endif // FORMATTED_REPORTER_H
