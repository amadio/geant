double boxsize[3];
double origin[3];
#include "math.h"

bool contains( const double * point )
{
  for( int unsigned dir=0; dir < 3; ++dir )
    {
      if ( fabs ( point[dir]-origin[dir] ) > boxsize[dir] ) return false;
    }
  return true;
}


void contains_v( const double * point, bool * isin, int np )
{
  for( unsigned int k=0; k < np; ++k)
    {
      for( unsigned int dir=0; dir < 3; ++dir )
	{
	  if ( fabs ( point[3*k+dir]-origin[dir] ) > boxsize[dir] ) isin[k]=false;
	}
      isin[k]=true;
    }
}


void contains_v2( const double * point, bool * isin, int np )
{
  for( unsigned int k=0; k < np; ++k)
    {
      for( unsigned int dir=0; dir < 3; ++dir )
	{
	  if ( fabs ( point[3*k+dir]-origin[dir] ) > boxsize[dir] ) isin[k]=false;
	}
      isin[k]=true;
    }
}


void contains_v3( const double * point, bool * isin, int np )
{
  for( unsigned int k=0; k < np; ++k)
    {
      for( unsigned int dir=0; dir < 3; ++dir )
	{
	  if ( fabs ( point[3*k+dir]-origin[dir] ) > boxsize[dir] ) isin[k]=false;
	}
      isin[k]=true;
    }
}


void contains_v4( const double * point, bool * isin, int np )
{
  #pragma ivdep
  for( unsigned int k=0; k < np; ++k)
    {
      bool result=true;
      for( unsigned int dir=0; dir < 3; ++dir )
	{
	  if ( fabs ( point[3*k+dir]-origin[dir] ) > boxsize[dir] ) result=false;
	}
      isin[k]=result;
    }
}


typedef struct
{
  double *coord[3];
} P; 


void contains_v5( const P & point, bool * isin, int np )
{
  #pragma ivdep
  for( unsigned int k=0; k < np; ++k)
    {
      bool result=true;
      for( unsigned int dir=0; dir < 3; ++dir )
	{
	  if (fabs (point.coord[dir][k]-origin[dir]) > boxsize[dir]) result=false;
	}
      isin[k]=result;
    }
}


void contains_v6( const P & __restrict__ point, bool * __restrict__ isin, int np )
{
  #pragma ivdep
  for( unsigned int k=0; k < np; ++k)
    {
      bool result[3];
      for( unsigned int dir=0; dir < 3; ++dir )
	{
	  result[dir] = (fabs (point.coord[dir][k]-origin[dir]) > boxsize[dir]);
	}
      isin[k]=result[0] & result[1] & result[2];
    }
}



void contains_v7( const P & __restrict__ point, bool * __restrict__ isin, int np )
{
  for( unsigned int k=0; k < np; ++k)
    {
      bool resultx=(fabs (point.coord[0][k]-origin[0]) > boxsize[0]);
      bool resulty=(fabs (point.coord[1][k]-origin[1]) > boxsize[1]);
      bool resultz=(fabs (point.coord[2][k]-origin[2]) > boxsize[2]);
      isin[k]=resultx & resulty & resultz;
    }
}


