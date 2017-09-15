#ifndef _FIELDTRACK_H_
#define _FIELDTRACK_H_

/*---------------------
Data structure in place of GUFieldTrack to be used 
for input and output stream arrays of AccurateAdvance 
in IntegrationDriver. 

Functions DumpToArray and LoadFromArray can be removed 
          if this is changed into a struct - i.e. PosMomVector is made public data member.
Same goes for SetCurveLength and GetCurveLength functions.
----------------*/
#include <iostream>

struct FieldTrack{

public: 
  static constexpr int NumCompFT = 6;  // Number of components

  // Constructors
  FieldTrack() : fDistanceAlongCurve(0.0) { LoadZeroes(); } 
  FieldTrack( double PositionMomentum[NumCompFT], double length= 0.0) : fDistanceAlongCurve(length)
                        { LoadFromArray( PositionMomentum ); }
  FieldTrack( std::vector<double> PositionMomentumVec, double length= 0.0) : fDistanceAlongCurve(length)
                        { LoadFromVector( PositionMomentumVec ); }   
  ~FieldTrack(){};

  double GetComponent( int i ) const {
     // assert( 0 <= i && i < NumCompFT );
     return fPosMomArr[i];
  }

  void   SetComponent( int i, double val ) {
     // assert( 0 <= i && i < NumCompFT );
     fPosMomArr[i] = val;
  }   
  
  // Access & set methods
  void DumpToArray(double valArr[]) const {
    for (int i = 0; i < NumCompFT; ++i)        //and GUFieldTrack function
    {
      valArr[i] = fPosMomArr[i];
    }
  }

  void LoadFromArray(const double valArr[], int noVarsIntegrated = -1 )
  {
    if( noVarsIntegrated == -1 )
       noVarsIntegrated= 6; // NumCompFT;
    int top= std::min( noVarsIntegrated, 6 );  // NumCompFT ); 
    for (int i = 0; i < top; ++i)
    {
      fPosMomArr[i] = valArr[i];
    }
  }

  void LoadFromVector( const std::vector<double> valVec, double valRest = 0.0 )
  {
    int top= std::min( (int)(valVec.size()) , NumCompFT ); 
    for (int i = 0; i < top; ++i)
    {
      fPosMomArr[i] = valVec[i];
    }
    for (int i= top; i < NumCompFT; ++i)
    {
      fPosMomArr[i] = valRest;    //  Fill the rest, if any
    }
  }  

  void LoadZeroes()
  {
    for (int i = NumCompFT; i >=0; --i)
       fPosMomArr[i] = 0.0;
  }

  void   SetCurveLength(double len){ fDistanceAlongCurve = len; }
  double GetCurveLength() const { return fDistanceAlongCurve; }

private: 
  //data members   
  double fDistanceAlongCurve = 0.0;
  double fPosMomArr[NumCompFT];

public:
  
  friend std::ostream&
          operator<<( std::ostream& os, const FieldTrack& fieldTrack)
          {
            os<< " ( ";
            os<< " X= "<< fieldTrack.fPosMomArr[0]<<" "
                       << fieldTrack.fPosMomArr[1]<<" "
                       << fieldTrack.fPosMomArr[2]<<" "; //Position
            os<< " P= "<< fieldTrack.fPosMomArr[3]<<" "
                       << fieldTrack.fPosMomArr[4]<<" "
                       << fieldTrack.fPosMomArr[5]<<" "; //Momentum
            os<< " ) ";

            return os;
          }
};
#endif 
