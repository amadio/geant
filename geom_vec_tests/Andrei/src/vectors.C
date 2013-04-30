#include <vector>
//#include <xmmintrin.h>
#include <emmintrin.h>
//#include <pmmintrin.h>
#include <cmath>
#include "TStopwatch.h"
#include "TRandom.h"

//typedef double v4sf __attribute__ ((vector_size (16))); // vector of four single doubles
//#define DO__COMPILE
//#define ALIGN

void vectors(int);

int main(int argc, char** argv)                                                                                                                                                                                  
{
   vectors(1000000);
}
   

const union __attribute__((aligned(16)))
{
       Long64_t i[2];
       __m128d d;
} absMask = {0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF};

//______________________________________________________________________________
struct point_t
{
  double x[3] __attribute__((aligned(16)));
};

//______________________________________________________________________________
struct vpoint16_t
{
// Structure of coordinate arrays for 16 points
  int    np;          // number of points
  double x[16] __attribute__((aligned(16)));
  double y[16] __attribute__((aligned(16)));
  double z[16] __attribute__((aligned(16)));
  
  void GetPoint(int i, point_t &point) {point.x[0]=x[i]; point.x[1]=y[i]; point.x[2]=z[i];}
  void LoadPoints(int n, point_t *array);
};

void vpoint16_t::LoadPoints(int n, point_t * __restrict__ array)
{
// Load n points from an array of points. Gets vectorized
  np = n;
  for (int i=0; i<n; i++) {
    x[i] = array[i].x[0];
    y[i] = array[i].x[1];
    z[i] = array[i].x[2];
  }  
}    

//______________________________________________________________________________
struct pointvec_t
{
// Structure of arrays of 4 point coordinates.
  double x[4] __attribute__((aligned(16)));
  double y[4] __attribute__((aligned(16)));
  double z[4] __attribute__((aligned(16)));
  
  pointvec_t(const point_t * __restrict__ points);
  void Fill(const point_t * __restrict__ points);
};       

pointvec_t::pointvec_t(const point_t * __restrict__ points)
{
// Copy the array of 4 points into the structure
   Fill(points);
}
   
void pointvec_t::Fill(const point_t * __restrict__ points)
{
// Copy the array of 4 points into the structure
// The points array is assumed contiguous
  const double * __restrict__ crt = points[0].x;
  for (int i=0; i<4; i++) {
    x[i] = crt[3*i];
    y[i] = crt[3*i+1];
    z[i] = crt[3*i+2];
  }
}

//______________________________________________________________________________
Bool_t ContainsEarlyReturn(point_t &a)
{
// Contains version doing early returns but having to pay with one function
// call per point.
   const double dx = 10.;
   const double dy = 20.;
   const double dz = 5.;
   if (fabs(a.x[0])>dx) return false;   
   if (fabs(a.x[1])>dy) return false;   
   if (fabs(a.x[2])>dz) return false;   
   return true;
}

//______________________________________________________________________________
void ContainsNonVect(int n, point_t *a, bool * inside)
{
// Non-vectorizable example, but taking a vector input.
   const double dx = 10.;
   const double dy = 20.;
   const double dz = 5.;
   int i;
   for (i=0; i<n; i++) {
     inside[i] =  (fabs(a[i].x[0])<=dx) & (fabs(a[i].x[1])<=dy) & (fabs(a[i].x[2])<=dz);
   }  
}

//______________________________________________________________________________
void ContainsVect(int n, point_t * __restrict__ a, bool * __restrict__ inside)
{
// Vectorized Box::contains method.
   const double dx = 10.;
   const double dy = 20.;
   const double dz = 5.;
   int i;
   for (i=0; i<n; i++) {
      double *temp = a[i].x;
      inside[i] =  (fabs(*temp)<dx) & (fabs(*(temp+1))<dy) & (fabs(*(temp+2))<dz);
  }   
}

//______________________________________________________________________________
void ContainsVect(vpoint16_t &va, bool * __restrict__ inside)
{
// Vectorized Box::contains method.
   const double dx = 10.;
   const double dy = 20.;
   const double dz = 5.;
   int i;
   for (i=0; i<va.np; i++) {
//      double *temp = a[i].x;
      inside[i] =  (fabs(va.x[i])<dx) & (fabs(va.y[i])<dy) & (fabs(va.z[i])<dz);
  }   
}

/*
//______________________________________________________________________________
void ContainsIntrinsics(int n, point_t * __restrict__ a, bool * __restrict__ inside)
{
// Contains method using intrinsics
   static const double dd[4] = {10., 20., 5., 0.};
   unsigned i;
   __m128d vboxxy = _mm_load_pd(dd);
   __m128d vboxzt = _mm_load_pd(dd+2);
   __m128d vecxy, veczt, res;
   double rr;
   double temp[2];
   double *addr;
//   for (i=0; i<=n-2; i += 2) {
   for (i=0; i<n; i++) {
      addr = a[i].x;
//      printf("a=( %f, %f, %f)\n", *addr, *(addr+1), *(addr+2));
      vecxy = _mm_load_pd(addr);        // vecxy = (x,y)
//      _mm_store_pd(temp, vecxy);
//      printf("stored xy: (%f, %f)\n", temp[0], temp[1]);
      vecxy = _mm_and_pd(vecxy, absMask.d);  // vecxy = (abs(x), abs(y))
//      _mm_store_pd(temp, vecxy);
//      printf("abs xy: (%f, %f)\n", temp[0], temp[1]);
      vecxy = _mm_cmpgt_pd(vecxy, vboxxy); // ((|x|<=dx)?1:0, /y/<=dy?1:0)
//      _mm_store_pd(temp, vecxy);
//      printf("after compare xy: (%f, %f)\n", temp[0], temp[1]);
      veczt = _mm_load_pd(addr+2);        // Same for z
//      _mm_store_pd(temp, veczt);
//      printf("stored zt: (%f, %f)\n", temp[0], temp[1]);
      veczt = _mm_and_pd(veczt, absMask.d);
//      _mm_store_pd(temp, veczt);
//      printf("abs zt: (%f, %f)\n", temp[0], temp[1]);
      veczt = _mm_cmpgt_pd(veczt, vboxzt); // ((|z|<=dz)?1:0, 1)
//      _mm_store_pd(temp, veczt);
//      printf("after compare zt: (%f, %f)\n", temp[0], temp[1]);
      res   = _mm_or_pd(vecxy, veczt);    // Make bitwise OR
//      _mm_store_pd(temp, res);
//      printf("bitwise OR xy zt: (%f, %f)\n", temp[0], temp[1]);
      res   = _mm_hadd_pd(res, res);
//      _mm_store_pd(temp, res);
//      printf("xy + zt: (%f, %f)\n", temp[0], temp[1]);
      _mm_storeh_pd(&rr, res);
      inside[i] = rr;
      inside[i] = !inside[i];
//      if (inside[i]) printf("INSIDE\n");
  }   
}
*/   

//______________________________________________________________________________
//#ifdef DO__COMPILE
void vectors(int nloop)
{
   TStopwatch timer;
   int i,j, checksum;
   point_t points[1024];
   vpoint16_t vpoints[64] __attribute__((aligned(16)));
   bool inside[1024];
   double *crt;

   // Fill the vector of points
   double rnd[4] __attribute__((aligned(16)));
   for (i=0; i<1024; i++) {
     inside[i] = false;
     gRandom->RndmArray(3, rnd);
     crt = points[i].x;
     for (j=0; j<3; j++) crt[j] = 15.*(1.-2.*rnd[j]);
   }  
   for (i=0; i<64; i++) {
     j = 16*i;
     vpoints[i].LoadPoints(16, &points[j]);
   }  
   printf("Testing %d points...\n", 1024*nloop);
   timer.Start(kTRUE);
   for (i=0; i<nloop; i++) {
     for (j=0; j<1024; j++) inside[j] = ContainsEarlyReturn(points[j]);
   }  
   printf("ContainsEarlyReturn input= one point (x,y,z) per call:\n");
   timer.Stop();
   timer.Print();  
   checksum = 0;
   for (j=0; j<1024; j++) {checksum += inside[j];}
   memset(inside, 0, 1024*sizeof(bool));
   printf("  checksum=%d\n", checksum);

   timer.Start(kTRUE);
   for (i=0; i<nloop; i++) {
     for (j=0; j<64; j++) ContainsNonVect(16, points+16*j, inside+16*j);
   }  
   printf("ContainsNonVect input=array of 16 points per call:\n");
   timer.Stop();
   timer.Print();  
   checksum = 0;
   for (j=0; j<1024; j++) {checksum += inside[j];}
   memset(inside, 0, 1024*sizeof(bool));
   printf("  checksum=%d\n", checksum);

   timer.Start(kTRUE);
   for (i=0; i<nloop; i++) {
     for (j=0; j<64; j++) ContainsVect(16, points+16*j, inside+16*j);
   }  
   printf("ContainsVect input=array of 16 points per call, autovectorized:\n");
   timer.Stop();
   timer.Print();  
   checksum = 0;
   for (j=0; j<1024; j++) {checksum += inside[j];}
   memset(inside, 0, 1024*sizeof(bool));
   printf("  checksum=%d\n", checksum);
  
   timer.Start(kTRUE);
   for (i=0; i<nloop; i++) {
     for (j=0; j<64; j++) ContainsVect(vpoints[j], inside+16*j);
   }  
   printf("ContainsVect input=struct of arrays (x[16,y[16,z[16]) autovectorized:\n");
   timer.Stop();
   timer.Print();  
   checksum = 0;
   for (j=0; j<1024; j++) {checksum += inside[j];}
   memset(inside, 0, 1024*sizeof(bool));
   printf("  checksum=%d\n", checksum);

/*
   timer.Start(kTRUE);
   for (i=0; i<nloop*64; i++) ContainsIntrinsics(16, points, inside);
   printf("ContainsIntrinsics:\n");
   timer.Stop();
   timer.Print();  
   checksum = 0;
   for (j=0; j<1024; j++) {checksum += inside[j]; inside[j] = false;}
   printf("  checksum=%d\n", checksum);
*/   
   printf("ContainsIntrinsics: coming soon\n");
}
//#endif
