#include "GUTrackHandler.h"
#include "GUConstants.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "mm_malloc.h"

#include "base/Global.h"
#include "GUAuxFunctions.h"      // Define sincos on Apple, etc

namespace vecphys {

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

#undef NDEBUG
// Makes sure that assert are used, even in Release mode

void GUTrackHandler::Allocate(size_t nTracks)
{
  const int blockSize = 64;  // Bytes

  SetNumberOfTracks(nTracks);

  if(fNumberOfTracks > 0) {

    unsigned int numAlloc = ((fNumberOfTracks-1) / blockSize + 1 ) * blockSize;

    // std::cout << " Allocating " << numAlloc << " tracks - vs actual number= " << fNumberOfTracks << std::endl;
    unsigned int memSizeAlloc= numAlloc*sizeof(GUTrack);

    //allocation for aos
    fTrack_aos = (GUTrack *)_mm_malloc (memSizeAlloc, blockSize);
    assert( (long) fTrack_aos % blockSize == 0 );  // Double check

    memSizeAlloc += blockSize;
    //allocation for soa
    char* soaBuffer = (char *) _mm_malloc (
            memSizeAlloc,  // sizeof(int)+fNumberOfTracks*sizeof(GUTrack),
					  blockSize);
    assert( (long) soaBuffer % blockSize == 0 );  // Double check
    fBuffer = soaBuffer;

    // Single integer
    assert( blockSize >= sizeof(int) );
    assert( blockSize % sizeof(int) == 0);

    //stride for GUTrack_v.numTracks
    soaBuffer += blockSize;  // std::max( blockSize, sizeof(int) );
    assert( (long) soaBuffer % blockSize == 0 );  // Double check

    const int offset_int    = numAlloc*sizeof(int);
    const int offset_double = numAlloc*sizeof(double);

    // Ensures subsequent conditions hold
    assert( offset_int    >= blockSize );
    assert( offset_double >= blockSize );
    assert( offset_int    % blockSize == 0);
    assert( offset_double % blockSize == 0);

    //set ptr to each element of GUTrack_v
    fTrack_soa.status       = (int*)soaBuffer;
    soaBuffer += offset_int;
    assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.particleType = (int*) soaBuffer;
    soaBuffer += offset_int;
    // assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.id           = (int*) soaBuffer;
    soaBuffer += offset_int;
    // assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.parentId     = (int*) soaBuffer;
    soaBuffer += offset_int;
    // assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.proc         = (int*) soaBuffer;
    soaBuffer += offset_double;
    assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.x            = (double*) soaBuffer;
    soaBuffer += offset_double;
    // assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.y            = (double*) soaBuffer;
    soaBuffer += offset_double;
    // assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.z            = (double*) soaBuffer;
    soaBuffer += offset_double;
    // assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.px           = (double*) soaBuffer;
    soaBuffer += offset_double;
    // assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.py           = (double*) soaBuffer;
    soaBuffer += offset_double;
    // assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.pz           = (double*) soaBuffer;
    soaBuffer += offset_double;
    // assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.E            = (double*) soaBuffer;
    soaBuffer += offset_double;
    // assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.q            = (double*) soaBuffer;
    soaBuffer += offset_double;
    // assert( (long) soaBuffer % blockSize == 0 );

    fTrack_soa.s            = (double*) soaBuffer;

    // std::cout << " soaBuffer start = " << (long) soaBuffer  << std::endl;
    // std::cout << "   Trial end     = " << ((long) soaBuffer + offset_double ) << std::endl;
    // std::cout << "   Real  end     = " << ((long) fBuffer + memSizeAlloc ) << std::endl;

    assert( (long) soaBuffer + offset_double <= (long) fBuffer + memSizeAlloc );

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

  //constants - move to a header file
  double const pi = 3.14159265358979323846;
  double const ecalRmim = 1290.;
  double const ecalRmax = 1520.;
  double const ecalZmax = 3000.;

  fTrack_soa.numTracks = fNumberOfTracks;

  for(size_t i = 0 ; i < fNumberOfTracks ; ++i){
    double rho, z, p, mass;
    double theta, phi;
    double cosphi, sinphi;
    double sintheta, tantheta, costheta;

    rho = ecalRmim + (ecalRmax-ecalRmim)*UniformRandom<Real_t>(0,-1);

    if( minP == maxP )
    {
       p= maxP;
    }
    else
    {
       do {
	 p = minP - 0.2*(maxP - minP)*math::Log(UniformRandom<Real_t>(0,-1));
       }
       while (p>maxP);
    }

    z = ecalZmax*(2*UniformRandom<Real_t>(0,-1)-1.0);
    phi = 2*pi*UniformRandom<Real_t>(0,-1);
    tantheta = rho/z;
    theta = std::atan(tantheta); // (rho/z);

    // cosphi = std::cos(phi);
    // sinphi = std::sin(phi);
    sincos(phi,   &sinphi,   &cosphi);
    sincos(theta, &sintheta, &costheta);
    // sintheta = std::sin(theta);

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
    (fTrack_soa.pz)[i] = fTrack_aos[i].pz = p*costheta; // std::cos(theta);

    mass = 0; // electron_mass_c2*fTrack_aos[i].q*fTrack_aos[i].q;
    (fTrack_soa.E)[i]  = fTrack_aos[i].E  = p*p/(math::Sqrt(p*p + mass*mass) + mass);

  }
}

//utility functions - can be elsewhere

void GUTrackHandler::SortAoSTracksByEnergy(GUTrack* AoS, size_t numberOfTracks)
{
  //sort AoS tracks by energy in ascending order.
  std::sort(AoS, AoS+numberOfTracks, [](GUTrack const &a, GUTrack const &b){ return a.E < b.E; });
}

void GUTrackHandler::SortSoATracksByEnergy(GUTrack_v& SoA, size_t numberOfTracks)
{
  //temporary AoS tracks
  const int blockSize = 64;  // Bytes
  GUTrack* AoS = (GUTrack *)_mm_malloc (numberOfTracks*sizeof(GUTrack), blockSize);

  //sort SoA tracks by energy in ascending order.
  CopySoATracksToAoS(SoA,AoS,numberOfTracks);
  SortAoSTracksByEnergy(AoS ,numberOfTracks);
  CopyAoSTracksToSoA(AoS,SoA,numberOfTracks);

  _mm_free(AoS);
}

void GUTrackHandler::CopyAoSTracks(GUTrack* fromAoS, GUTrack* toAoS, size_t numberOfTracks)
{
  for(size_t i = 0 ; i < numberOfTracks ; ++i){
    toAoS[i].status       = fromAoS[i].status       ;
    toAoS[i].proc         = fromAoS[i].proc         ;
    toAoS[i].particleType = fromAoS[i].particleType ;
    toAoS[i].id           = fromAoS[i].id           ;
    toAoS[i].parentId     = fromAoS[i].parentId     ;
    toAoS[i].x            = fromAoS[i].x            ;
    toAoS[i].y            = fromAoS[i].y            ;
    toAoS[i].z            = fromAoS[i].z            ;
    toAoS[i].q            = fromAoS[i].q            ;
    toAoS[i].s            = fromAoS[i].s            ;
    toAoS[i].px           = fromAoS[i].px           ;
    toAoS[i].py           = fromAoS[i].py           ;
    toAoS[i].pz           = fromAoS[i].pz           ;
    toAoS[i].E            = fromAoS[i].E            ;
  }
}

void GUTrackHandler::CopySoATracks(GUTrack_v& fromSoA, GUTrack_v& toSoA, size_t numberOfTracks)
{
  toSoA.numTracks = fromSoA.numTracks;
  for(size_t i = 0 ; i < numberOfTracks ; ++i){
    (toSoA.status)[i]       = (fromSoA.status)[i]       ;
    (toSoA.proc)[i]         = (fromSoA.proc)[i]         ;
    (toSoA.particleType)[i] = (fromSoA.particleType)[i] ;
    (toSoA.id)[i]           = (fromSoA.id)[i]           ;
    (toSoA.parentId)[i]     = (fromSoA.parentId)[i]     ;
    (toSoA.x)[i]            = (fromSoA.x)[i]            ;
    (toSoA.y)[i]            = (fromSoA.y)[i]            ;
    (toSoA.z)[i]            = (fromSoA.z)[i]            ;
    (toSoA.q)[i]            = (fromSoA.q)[i]            ;
    (toSoA.s)[i]            = (fromSoA.s)[i]            ;
    (toSoA.px)[i]           = (fromSoA.px)[i]           ;
    (toSoA.py)[i]           = (fromSoA.py)[i]           ;
    (toSoA.pz)[i]           = (fromSoA.pz)[i]           ;
    (toSoA.E)[i]            = (fromSoA.E)[i]            ;
  }
}

void GUTrackHandler::CopyAoSTracksToSoA(GUTrack* fromAoS, GUTrack_v& toSoA, size_t numberOfTracks)
{
  toSoA.numTracks = numberOfTracks;
  for(size_t i = 0 ; i < numberOfTracks ; ++i){
    (toSoA.status)[i]        = fromAoS[i].status       ;
    (toSoA.proc)[i]          = fromAoS[i].proc         ;
    (toSoA.particleType)[i]  = fromAoS[i].particleType ;
    (toSoA.id)[i]            = fromAoS[i].id           ;
    (toSoA.parentId)[i]      = fromAoS[i].parentId     ;
    (toSoA.x)[i]             = fromAoS[i].x            ;
    (toSoA.y)[i]             = fromAoS[i].y            ;
    (toSoA.z)[i]             = fromAoS[i].z            ;
    (toSoA.q)[i]             = fromAoS[i].q            ;
    (toSoA.s)[i]             = fromAoS[i].s            ;
    (toSoA.px)[i]            = fromAoS[i].px           ;
    (toSoA.py)[i]            = fromAoS[i].py           ;
    (toSoA.pz)[i]            = fromAoS[i].pz           ;
    (toSoA.E)[i]             = fromAoS[i].E            ;
  }
}

void GUTrackHandler::CopySoATracksToAoS(GUTrack_v& fromSoA, GUTrack* toAoS, size_t numberOfTracks)
{
  for(size_t i = 0 ; i < numberOfTracks ; ++i){
    toAoS[i].status       = (fromSoA.status)[i]       ;
    toAoS[i].proc         = (fromSoA.proc)[i]         ;
    toAoS[i].particleType = (fromSoA.particleType)[i] ;
    toAoS[i].id           = (fromSoA.id)[i]           ;
    toAoS[i].parentId     = (fromSoA.parentId)[i]     ;
    toAoS[i].x            = (fromSoA.x)[i]            ;
    toAoS[i].y            = (fromSoA.y)[i]            ;
    toAoS[i].z            = (fromSoA.z)[i]            ;
    toAoS[i].q            = (fromSoA.q)[i]            ;
    toAoS[i].s            = (fromSoA.s)[i]            ;
    toAoS[i].px           = (fromSoA.px)[i]           ;
    toAoS[i].py           = (fromSoA.py)[i]           ;
    toAoS[i].pz           = (fromSoA.pz)[i]           ;
    toAoS[i].E            = (fromSoA.E)[i]            ;
  }
}

} // end namespace vecphys
