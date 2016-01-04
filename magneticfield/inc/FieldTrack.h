#ifndef _FIELDTRACK_H_
#define _FIELDTRACK_H_

/*---------------------
Data structure in place of GUFieldTrack to be used 
for input and output stream arrays of AccurateAdvance 
in IntegrationDriver. Functions DumpToArray and LoadFromArray
can be removed if PosMomVector is made public data member.
Same goes for SetCurveLength and GetCurveLength functions.
----------------*/
#include <iostream>

struct FieldTrack{

public: 
  static constexpr int NumComp = 6;  // Number of components

  // Constructors
  FieldTrack() : fDistanceAlongCurve(0.0) { LoadZeroes(); } 
  FieldTrack( double PositionMomentum[NumComp], double length= 0.0) : fDistanceAlongCurve(length)
                        { LoadFromArray( PositionMomentum ); }
  FieldTrack( std::vector<double> PositionMomentumVec, double length= 0.0) : fDistanceAlongCurve(length)
                        { LoadFromVector( PositionMomentumVec ); }   
  ~FieldTrack(){};

  double GetComponent( int i ) {
     // assert( 0 <= i && i < NumComp );
     return fPosMomArr[i];
  }

  void   SetComponent( int i, double val ) {
     // assert( 0 <= i && i < NumComp );
     fPosMomArr[i] = val;
  }   
  
  // Access & set methods
  void DumpToArray(double valArr[]) { //12 from ncompSVEC as in both TemplateGUIntegrationDriver
    for (int i = 0; i < NumComp; ++i)        //and GUFieldTrack function
    {
      valArr[i] = fPosMomArr[i];
    }
  }

  void LoadFromArray(const double valArr[], int noVarsIntegrated = -1 )
  {
    if( noVarsIntegrated == -1 )
       noVarsIntegrated= 6; // NumComp;
    int top= std::min( noVarsIntegrated, 6 );  // NumComp ); 
    for (int i = 0; i < top; ++i)
    {
      fPosMomArr[i] = valArr[i];
    }
  }

  void LoadFromVector( const std::vector<double> valVec, double valRest = 0.0 )
  {
    int top= std::min( (int)(valVec.size()) , NumComp ); 
    for (int i = 0; i < top; ++i)
    {
      fPosMomArr[i] = valVec[i];
    }
    for (int i= top; i < NumComp; ++i)
    {
      fPosMomArr[i] = valRest;    //  Fill the rest, if any
    }
  }  

  void LoadZeroes()
  {
    for (int i = NumComp; i >=0; --i)
       fPosMomArr[i] = 0.0;
  }

  void   SetCurveLength(double len){ fDistanceAlongCurve = len; }
  double GetCurveLength(){ return fDistanceAlongCurve; }

private: 
  //data members   
  double fDistanceAlongCurve = 0.0;
  double fPosMomArr[NumComp];

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
