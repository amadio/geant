#include "GXTrackHandler.h"
#include "GPConstants.h"
#include "GPRandom.h"
#include <cmath>

GXTrackHandler::GXTrackHandler()
  : theNumberOfTracks(0), theNumberOfElectrons(0), theNumberOfPhotons(0),
    track_h(0), track_e(0), track_g(0)
{
  thePhotonFraction = 0.0;
}

GXTrackHandler::GXTrackHandler(size_t nTracks, double photonFrac)
{
  Allocate(nTracks, photonFrac);
}

GXTrackHandler::~GXTrackHandler()
{
  Deallocate();
}

void GXTrackHandler::SetNumberOfTracks(size_t nTracks, double photonFrac)
{
  theNumberOfTracks = nTracks;
  thePhotonFraction = photonFrac;
  theNumberOfPhotons = theNumberOfTracks*thePhotonFraction;
  theNumberOfElectrons = theNumberOfTracks - theNumberOfPhotons;  
}

void GXTrackHandler::Deallocate()
{
  if(theNumberOfTracks > 0)    free(track_h);  
  if(theNumberOfElectrons > 0) free(track_e);  
  if(theNumberOfPhotons > 0)   free(track_g);  
}

void GXTrackHandler::Allocate(size_t nTracks, double photonFrac)
{
  SetNumberOfTracks(nTracks, photonFrac);

  if(theNumberOfTracks > 0) {
    track_h = (GXTrack *) malloc (theNumberOfTracks*sizeof(GXTrack));
  }
  if(theNumberOfElectrons > 0) {
    track_e = (GXTrack *) malloc (theNumberOfElectrons*sizeof(GXTrack));
  }
  if(theNumberOfPhotons > 0) {
    track_g = (GXTrack *) malloc (theNumberOfPhotons*sizeof(GXTrack));
  }
}

void GXTrackHandler::Reallocate(size_t nTracks, double photonFrac)
{
  Deallocate();
  Allocate(nTracks,photonFrac);
}

void GXTrackHandler::FillOneTrack(GXTrack* aTrack) 
{
  //@@@G4FWP - add an implementation in conjunction with TaskManager
  ;
}

void GXTrackHandler::GenerateRandomTracks(size_t nTracks, double photonFrac,
					  double minP, double maxP)
{
  Reallocate(nTracks,photonFrac);

  G4double rho, z, p, mass;
  G4double theta, phi;
  G4double cosphi, sinphi;
  G4double sintheta;

  for(size_t i = 0 ; i < theNumberOfTracks ; ++i){
    
    rho = ecalRmim + (ecalRmax-ecalRmim)*GPUniformRand(0,-1);
    p = minP + (maxP - minP)*GPUniformRand(0,-1);
    z = ecalZmax*(2*GPUniformRand(0,-1)-1.0);
    phi = twopi*GPUniformRand(0,-1);
    theta = std::atan(rho/z);
    cosphi = std::cos(phi);
    sintheta = std::sin(phi);
    
    track_h[i].status = 0;
    track_h[i].proc = -1;
    track_h[i].q = ( i < theNumberOfPhotons ) ? 0.0 : -1.0;
    track_h[i].x = rho*cosphi;
    track_h[i].y = rho*sinphi; 
    track_h[i].z = z; 
    track_h[i].s = 1.0; 
    
    track_h[i].px = p*sintheta*cosphi;
    track_h[i].py = p*sintheta*sinphi;
    track_h[i].pz = p*std::cos(theta);
    
    mass = electron_mass_c2*track_h[i].q*track_h[i].q;
    track_h[i].E  = p*p/(sqrt(p*p + mass*mass) + mass);

    if(i < theNumberOfPhotons) CopyTrack(&track_h[i], &track_g[i]);
    else                       CopyTrack(&track_h[i], &track_e[i]);

  }
}

void GXTrackHandler::CopyTrack(GXTrack *This, GXTrack *track)
{ 
  track->status = This->status;
  track->id     = This->id;
  track->proc   = This->proc;
  track->q      = This->q;
  track->x      = This->x;
  track->y      = This->y;
  track->z      = This->z;
  track->px     = This->px;
  track->py     = This->py;
  track->pz     = This->pz;
  track->E      = This->E;
  track->s      = This->s;
}
