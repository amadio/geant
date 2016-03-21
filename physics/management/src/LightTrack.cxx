#include "LightTrack.h"

using namespace geant;

//-------------------------------------
//--- LightTrack non-inline methods ---
//-------------------------------------


LightTrack::LightTrack() :
  fGVcode( -1 ), fGTrackIndex( -1 ), fMaterialCutCoupleIndex( -1 ), fTrackStatus( -1 ),
  fProcessIndex( 0 ), fTargetZ( 0 ), fTargetN( 0 ), 
  fXdir( 0.0 ), fYdir( 0.0 ), fZdir( 0.0 ), fKinE( -1.0 ), fMass( -1.0 ), 
  fTime( -1.0 ), fWeight( 1.0 ), fStepLength( 0.0 ), fEdep ( 0.0 ),
  fExtraInfo( nullptr ) 
{}


LightTrack::LightTrack( const int aGVcode, const int aGTrackIndex, 
                        const int aMaterialCutCoupleIndex,
                        const int aTrackStatus, const int aProcessIndex, 
                        const int aTargetZ, const int aTargetN, 
                        const double aXdir, const double aYdir, const double aZdir,
                        const double aKinE, const double aMass,  const double aTime,
                        const double aWeight, const double aStepLength, const double aEdep,
                        ExtraInfo *aExtraInfo ) :
  fGVcode( aGVcode ), fGTrackIndex( aGTrackIndex ), 
  fMaterialCutCoupleIndex( aMaterialCutCoupleIndex ), fTrackStatus( aTrackStatus ),
  fProcessIndex( aProcessIndex ), fTargetZ( aTargetZ ), fTargetN( aTargetN ), 
  fXdir( aXdir ), fYdir( aYdir ), fZdir( aZdir ), 
  fKinE( aKinE ), fMass( aMass ), fTime( aTime ), fWeight( aWeight ), 
  fStepLength( aStepLength ), fEdep ( aEdep ), fExtraInfo( aExtraInfo )
{}


LightTrack::LightTrack( const LightTrack &other ) :
  fGVcode( other.fGVcode ), fGTrackIndex( other.fGTrackIndex ), 
  fMaterialCutCoupleIndex( other.fMaterialCutCoupleIndex ), fTrackStatus( other.fTrackStatus ),
  fProcessIndex( other.fProcessIndex ), fTargetZ( other.fTargetZ ), fTargetN( other.fTargetN ),
  fXdir( other.fXdir ), fYdir( other.fYdir ), fZdir( other.fZdir ),  
  fKinE( other.fKinE ), fMass( other.fMass ), fTime( other.fTime ), fWeight( other.fWeight ),
  fStepLength( other.fStepLength ), fEdep( other.fEdep ), fExtraInfo( other.fExtraInfo )
{}


LightTrack& LightTrack::operator=( const LightTrack &other ) {
  if ( this != &other ) {
    fGVcode = other.fGVcode;
    fGTrackIndex = other.fGTrackIndex; 
    fMaterialCutCoupleIndex = other.fMaterialCutCoupleIndex;
    fTrackStatus = other.fTrackStatus;
    fProcessIndex = other.fProcessIndex;
    fTargetZ = other.fTargetZ;
    fTargetN = other.fTargetN;
    fXdir = other.fXdir;
    fYdir = other.fYdir;
    fZdir = other.fZdir;  
    fKinE = other.fKinE;
    fMass = other.fMass;
    fTime = other.fTime;
    fWeight = other.fWeight;
    fStepLength = other.fStepLength;
    fEdep = other.fEdep;
    SetExtraInfo( other.fExtraInfo );
  }
  return *this;
}


LightTrack::~LightTrack() {
  // Note: we assume that the light track has the ownership of the extra information.
  if ( fExtraInfo != nullptr ) delete fExtraInfo;
}


//---------------------------------------
//--- LightTrack_v non-inline methods ---
//---------------------------------------


LightTrack_v::LightTrack_v() : 
  fNtracks( 0 ), fGVcodeV( nullptr ), fGTrackIndexV( nullptr ),
  fMaterialCutCoupleIndexV( nullptr ), fTrackStatusV( nullptr ),
  fProcessIndexV( nullptr ), fTargetZV( nullptr ), fTargetNV( nullptr ),
  fXdirV( nullptr ), fYdirV( nullptr ), fZdirV( nullptr ), fKinEV( nullptr ),
  fMassV( nullptr ), fTimeV( nullptr ), fWeightV( nullptr ), 
  fStepLengthV( nullptr), fEdepV( nullptr), fExtraInfoV( nullptr )
{}


void LightTrack_v::GetTrack( const int i, LightTrack &aLightTrack ) const {
  if ( i >= fNtracks ) return;
  aLightTrack.SetGVcode( fGVcodeV[ i ] );
  aLightTrack.SetTrackIndex( fGTrackIndexV[ i ] ); 
  aLightTrack.SetMaterialCutCoupleIndex( fMaterialCutCoupleIndexV[ i ] );
  aLightTrack.SetTrackStatus( fTrackStatusV[ i ] );
  aLightTrack.SetProcessIndex( fProcessIndexV[ i ] );
  aLightTrack.SetTargetZ( fTargetZV[ i ] );
  aLightTrack.SetTargetN( fTargetNV[ i ] );
  aLightTrack.SetDirX( fXdirV[ i ] );
  aLightTrack.SetDirY( fYdirV[ i ] );
  aLightTrack.SetDirZ( fZdirV[ i ] );  
  aLightTrack.SetKinE( fKinEV[ i ] );
  aLightTrack.SetMass( fMassV[ i ] );
  aLightTrack.SetTime( fTimeV[ i ] );
  aLightTrack.SetWeight( fWeightV[ i ] );
  aLightTrack.SetStepLength( fStepLengthV[ i ] );
  aLightTrack.SetEnergyDeposit( fEdepV[ i ] );
  aLightTrack.SetExtraInfo( fExtraInfoV[ i ] );
}


void LightTrack_v::AddTrack( LightTrack &aLightTrack ) {
  int itrack = fNtracks;
  fGVcodeV[ itrack ] = aLightTrack.GetGVcode();
  fGTrackIndexV[ itrack ] = aLightTrack.GetTrackIndex(); 
  fMaterialCutCoupleIndexV[ itrack ] = aLightTrack.GetMaterialCutCoupleIndex();
  fTrackStatusV[ itrack ] = aLightTrack.GetTrackStatus();
  fProcessIndexV[ itrack ] = aLightTrack.GetProcessIndex();
  fTargetZV[ itrack ] = aLightTrack.GetTargetZ();
  fTargetNV[ itrack ] = aLightTrack.GetTargetN();
  fXdirV[ itrack ] = aLightTrack.GetDirX();
  fYdirV[ itrack ] = aLightTrack.GetDirY();
  fZdirV[ itrack ] = aLightTrack.GetDirZ();  
  fKinEV[ itrack ] = aLightTrack.GetKinE();
  fMassV[ itrack ] = aLightTrack.GetMass();
  fTimeV[ itrack ] = aLightTrack.GetTime();
  fWeightV[ itrack ] = aLightTrack.GetWeight();
  fStepLengthV[ itrack ] = aLightTrack.GetStepLength();
  fEdepV[ itrack ] = aLightTrack.GetEnergyDeposit();
  fExtraInfoV[ itrack ] = aLightTrack.GetExtraInfo();
  fNtracks++;
}

