#include "GUTrackHandler.h"
#include "GURandom.h"
#include <cmath>
#include <iostream>
#include "mm_malloc.h"

GUTrackHandler::GUTrackHandler()
  : fNumberOfTracks(0),
    fTrack_aos(0),
    fBuffer(0)
{
}

GUTrackHandler::GUTrackHandler(size_t nTracks)
{
  Allocate(nTracks);
}

GUTrackHandler::~GUTrackHandler()
{
  Deallocate();
}

void GUTrackHandler::SetNumberOfTracks(size_t nTracks)
{
  fNumberOfTracks = nTracks;
}

void GUTrackHandler::Deallocate()
{
  if(fNumberOfTracks > 0) {
    _mm_free(fTrack_aos);  
    _mm_free(fBuffer);  
  }
}

void GUTrackHandler::Allocate(size_t nTracks)
{
  SetNumberOfTracks(nTracks);

  if(fNumberOfTracks > 0) {
    fTrack_aos = (GUTrack *)_mm_malloc (fNumberOfTracks*sizeof(GUTrack),32);

    char* fBuffer = (char *) _mm_malloc (sizeof(int)+
					fNumberOfTracks*sizeof(GUTrack),
					32);

    const int offset_int = fNumberOfTracks*sizeof(int);
    const int offset_double = fNumberOfTracks*sizeof(double);

    fBuffer += sizeof(int);
    
    fTrack_soa.status       = (int*)fBuffer; 
    fBuffer += offset_int;

    fTrack_soa.particleType = (int*) fBuffer;
    fBuffer += offset_int;

    fTrack_soa.id           = (int*) fBuffer;
    fBuffer += offset_int;

    fTrack_soa.parentId     = (int*) fBuffer;
    fBuffer += offset_int;

    fTrack_soa.proc         = (int*) fBuffer;
    fBuffer += offset_double;

    fTrack_soa.x            = (double*) fBuffer;
    fBuffer += offset_double;

    fTrack_soa.y            = (double*) fBuffer;
    fBuffer += offset_double;

    fTrack_soa.z            = (double*) fBuffer;
    fBuffer += offset_double;

    fTrack_soa.px           = (double*) fBuffer;
    fBuffer += offset_double;

    fTrack_soa.py           = (double*) fBuffer;
    fBuffer += offset_double;

    fTrack_soa.pz           = (double*) fBuffer;
    fBuffer += offset_double;

    fTrack_soa.E            = (double*) fBuffer;
    fBuffer += offset_double;

    fTrack_soa.q            = (double*) fBuffer;
    fBuffer += offset_double;

    fTrack_soa.s            = (double*) fBuffer;
  }
}

void GUTrackHandler::Reallocate(size_t nTracks)
{
  Deallocate();
  Allocate(nTracks);
}

void GUTrackHandler::FillOneTrack(GUTrack* aTrack) 
{
  //@@@G4FWP - add an implementation in conjunction with TaskManager
  ;
}

void GUTrackHandler::GenerateRandomTracks(size_t nTracks, 
					  double minP, double maxP)
{
  Reallocate(nTracks);

  double rho, z, p, mass;
  double theta, phi;
  double cosphi, sinphi;
  double sintheta;

  //constants - move to a header file
  double const electron_mass_c2 = 0.510998910;
  double const pi = 3.14159265358979323846;
  double const ecalRmim = 1290.; 
  double const ecalRmax = 1520.;
  double const ecalZmax = 3000.;

  fTrack_soa.numTracks = fNumberOfTracks;

  for(size_t i = 0 ; i < fNumberOfTracks ; ++i){
    
    rho = ecalRmim + (ecalRmax-ecalRmim)*GUUniformRand(0,-1);
    do {
      p = minP - 0.2*(maxP - minP)*log(GUUniformRand(0,-1));
    }
    while (p>maxP);

    z = ecalZmax*(2*GUUniformRand(0,-1)-1.0);
    phi = 2*pi*GUUniformRand(0,-1);
    theta = std::atan(rho/z);
    cosphi = std::cos(phi);
    sinphi = std::sin(phi);
    sintheta = std::sin(theta);
    
    (fTrack_soa.status)[i]       = fTrack_aos[i].status = 0;
    (fTrack_soa.proc)[i]         = fTrack_aos[i].proc = -1;
    (fTrack_soa.particleType)[i] = fTrack_aos[i].particleType = -1;
    (fTrack_soa.id)[i]           = fTrack_aos[i].id = i+1;
    (fTrack_soa.parentId)[i]     = fTrack_aos[i].parentId = 0;
    (fTrack_soa.x)[i]            = fTrack_aos[i].x = rho*cosphi;
    (fTrack_soa.y)[i]            = fTrack_aos[i].y = rho*sinphi; 
    (fTrack_soa.z)[i]            = fTrack_aos[i].z = z; 

    (fTrack_soa.q)[i]            = fTrack_aos[i].q = -11;
    (fTrack_soa.s)[i]            = fTrack_aos[i].s = 1.0; 
    		    		    
    (fTrack_soa.px)[i] = fTrack_aos[i].px = p*sintheta*cosphi;
    (fTrack_soa.py)[i] = fTrack_aos[i].py = p*sintheta*sinphi;
    (fTrack_soa.pz)[i] = fTrack_aos[i].pz = p*std::cos(theta);
    
    mass = electron_mass_c2*fTrack_aos[i].q*fTrack_aos[i].q;
    (fTrack_soa.E)[i]  = fTrack_aos[i].E  = p*p/(sqrt(p*p + mass*mass) + mass);

  }
}
