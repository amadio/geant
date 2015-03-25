#include "GUTrackHandler.h"
#include "GUConstants.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "mm_malloc.h"

#include "backend/Backend.h"

// add the sincos function on MAC because sincos is not part of math.h
#ifdef __APPLE__ // possibly other conditions
inline void sincos(double x, double *s, double *c){
  __sincos(x,s,c);
}
#endif

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
  // const int blockInts = blockSize / sizeof(int);
  // const int blockDbls = blockSize / sizeof(double);
  // const int blockMin = blockSizeBytes / std::max( sizeof(int), sizeof(double));

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

    rho = ecalRmim + (ecalRmax-ecalRmim)*UniformRandom<Precision>(0,-1);

    if( minP == maxP )
    {
       p= maxP;
    }
    else
    {
       do {
	 p = minP - 0.2*(maxP - minP)*log(UniformRandom<Precision>(0,-1));
       }
       while (p>maxP);
    }

    z = ecalZmax*(2*UniformRandom<Precision>(0,-1)-1.0);
    phi = 2*pi*UniformRandom<Precision>(0,-1);
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
    (fTrack_soa.E)[i]  = fTrack_aos[i].E  = p*p/(sqrt(p*p + mass*mass) + mass);

  }
}

} // end namespace vecphys
