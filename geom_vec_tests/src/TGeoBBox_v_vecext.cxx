// // This is an experimental playground for using vecextensions for SIMD vectorization
// // typedef for vector extensions

typedef double vd4 __attribute__((vector_size(4*sizeof(double)))) __attribute__((aligned(32)));
typedef double vd2 __attribute__((vector_size(2*sizeof(double)))) __attribute__((aligned(32)));
typedef unsigned int vui2 __attribute__((vector_size(2*sizeof(unsigned int)))) __attribute__((aligned(32)));
typedef unsigned int vui4 __attribute__((vector_size(4*sizeof(unsigned int)))) __attribute__((aligned(32)));
typedef unsigned long int vli4 __attribute__((vector_size(4*sizeof(long int)))) __attribute__((aligned(32)));
static const vd4 zerod4={0.,0.,0.,0.};
static const vd4 oned4={1.,1.,1.,1.};
static const vli4 oneli4={1,1,1,1};


// // this is the actual kernel doing the computation with possibility of early return
// void TGeoBBox_v::DistFromOutside_VECEXT(StructOfCoord const & point, StructOfCoord const & dir, 
//                          Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t * stepmax, Double_t * distance, Int_t np)
// {
//   // note: we take stepmax as a maxstep PER particle
//   // ALIGNMENTSTUFF HERE (TODO: should be generated automatically, MACRO? )
// #ifndef __INTEL_COMPILER
//   const double * p[3];
//   p[0] = (const double *) __builtin_assume_aligned (point.x, 32); // 32 bit for AVX 
//   p[1] = (const double *) __builtin_assume_aligned (point.y, 32); 
//   p[2] = (const double *) __builtin_assume_aligned (point.z, 32); 
//   const double * d[3];
//   d[0] = (const double *) __builtin_assume_aligned (dir.x, 32); 
//   d[1] = (const double *) __builtin_assume_aligned (dir.y, 32); 
//   d[2] = (const double *) __builtin_assume_aligned (dir.z, 32); 
//   double * dist = (double *) __builtin_assume_aligned (distance, 32); 
//   const double * stepm = (const double *) __builtin_assume_aligned(stepmax,32); 
// #else
//   const double * p[3];
//   p[0] = (const double *) point.x; 
//   p[1] = (const double *) point.y; 
//   p[2] = (const double *) point.z; 
//   const double * d[3];
//   d[0] = (const double *) dir.x; 
//   d[1] = (const double *) dir.y; 
//   d[2] = (const double *) dir.z; 
//   double * dist = (double *) distance; 
// #endif

// }

// a fake and probably inefficient implementation
vd4 Absd4( vd4 const &x )
{
  vd4 r=(x<0);
  return x*r;
}

extern void foo();
 // this is the actual kernel doing the computation with possibility of early return
vd4 DistFromOutside_VECEXT_P4( vd4 const & x , vd4 const & y, vd4 const & z, vd4 const & dirx, vd4 const & diry, vd4 const & dirz,  
                          double dx, double dy, double dz, const vd4 origin[3], vd4 const & stepmax, vd4 const & distance )
 {
  // what do we have as vectors ?
   vli4 in; // so maybe this is going to be vd4 as well
   vd4 saf[3];
   vd4 newpt[3];
  
   vli4 factor=oneli4; // initializing all components to zero
   vd4 infactor=oned4; // will be zero or one depending if we are inside/outside

   vd4 parx = {dx,dx,dx,dx};
   vd4 pary = {dy,dy,dy,dy};
   vd4 parz = {dz,dz,dz,dz};
   foo();
   // unrolled above block manually: ( it would be nice to have a template unrool loop and a lambda function ?? )
   newpt[0] = x - origin[0];
   //saf[0] = Absd4(newpt[0])-parx; // can we profit here from abs function in array form?
   // factor = (saf[0]>=stepmax); 
   // in =  (saf[0] < zerod4);
   foo();
   newpt[1] = y - origin[1];
   //saf[1] = Absd4(newpt[1])-pary;
   // factor += (saf[1]>=stepmax);
   // in = in & (saf[1] < zerod4);
   foo();
   newpt[2] = z - origin[2];
   //saf[2] = Absd4(newpt[2])-parz;
   //factor += (saf[2]>=stepmax);
   // in = in & (saf[2] < zerod4);
   //in = (saf[0] < zerod4) + (saf[1] < zerod4) + (saf[2] < zerod4);
        
   return  newpt[0]+newpt[1]+newpt[2];
 }


// // SOA version of static method DistFromOutside treating 4 particles ( if HAVE AVX )
// void TGeoBBox_v::DistFromOutside_v(const StructOfCoord & __restrict__  point,const StructOfCoord & __restrict__ dir,
// 				     Double_t dx, Double_t dy, Double_t dz, const Double_t * __restrict__ origin, const Double_t * __restrict__ stepmax, Double_t * __restrict__ distance, Int_t np)
// {

//   Double_t par[3];
//   par[0] = dx;
//   par[1] = dy;
//   par[2] = dz;
// #pragma ivdep
//   for(unsigned int k=0;k<np;++k) // @EXPECTVEC
//      {
//        Bool_t in;
//        Double_t saf[3];
//        Double_t newpt[3];
  
//        Double_t factor=1.;
//        Double_t infactor; // will be zero or one depending if we are inside/outside

//        // unrolled above block manually: ( it would be nice to have a template unrool loop and a lambda function ?? )
//        newpt[0] = p[0][k] - origin[0];
//        saf[0] = fabs(newpt[0])-par[0];
//        factor = (saf[0]>=stepmax[k]) ? TGeoShape::Big() : 1.; // this might be done at the end
//        in = (saf[0]<0);
       
//        newpt[1] = p[1][k] - origin[1];
//        saf[1] = fabs(newpt[1])-par[1];
//        factor *= (saf[1]>=stepmax[k]) ? TGeoShape::Big() : 1.; // this might be done at the end
//        in = in & (saf[1]<0);
	 
//        newpt[2] = p[2][k] - origin[2];
//        saf[2] = fabs(newpt[2])-par[2];
//        factor *= (saf[2]>=stepmax[k]) ? TGeoShape::Big() : 1.; // this might be done at the end
//        in = in & (saf[2]<0);
 
//        infactor = (double) !in;
      
//        // NOW WE HAVE THE SAFETYS AND IF IN OUT
     
//        Double_t coord;
//        Double_t snxt[3]={TGeoShape::Big(),TGeoShape::Big(),TGeoShape::Big()};
     
//        /*
//        Int_t ibreak=0;
//        for (i=0; i<3; i++) {
// 	 if (saf[i]<0) continue;
// 	 if (newpt[i]*dir[i] >= 0) continue;
// 	 snxt = saf[i]/fabs(dir[i]);
// 	 ibreak = 0;
// 	 for (j=0; j<3; j++) {
// 	   if (j==i) continue;
// 	   coord=newpt[j]+snxt*dir[j];
// 	   if (fabs(coord)>par[j]) {
// 	     ibreak=1;
// 	     break;
// 	   }
// 	 }
// 	 if (!ibreak) return snxt;
// 	 }
//        */
       
//        // THE FOLLOWING IS AN UNROLLED VERSION OF ABOVE CONSTRUCT WITHOUT EARLY RETURNS

//        // i=0
//        Int_t hit0=0;
//        if ( saf[0] > 0 & newpt[0]*d[0][k] < 0 ) // if out and right direction
// 	 {
// 	   snxt[0] = saf[0]/fabs(d[0][k]); // distance to y-z face
	   
// 	   double coord1=newpt[1]+snxt[0]*d[1][k]; // calculate new y and z coordinate
// 	   double coord2=newpt[2]+snxt[0]*d[2][k];
// 	   hit0 = (fabs(coord1)>par[1] | fabs(coord2)>par[2])? 0 : 1; // 0 means miss, 1 means hit
// 	 }
            
//        Int_t hit1=0;
//        if ( saf[1] > 0 & newpt[1]*d[1][k] < 0 )
// 	 {
// 	   snxt[1] = saf[1]/fabs(d[1][k]);
	  
// 	   double coord0=newpt[0]+snxt[1]*d[0][k];
// 	   double coord2=newpt[2]+snxt[1]*d[2][k];
// 	   hit1 = (fabs(coord0)>par[0] | fabs(coord2)>par[2])? 0 : 1;
// 	 }
     
//        Int_t hit2=0;
//        if ( saf[2] > 0 & newpt[2]*d[2][k] < 0 )
// 	 {
// 	   snxt[2] = saf[2]/fabs(d[2][k]);
// 	   double coord0 = newpt[0]+snxt[2]*d[0][k];
// 	   double coord1 = newpt[1]+snxt[2]*d[1][k];
// 	   hit2 = (fabs(coord0)>par[0] | fabs(coord1)>par[1])? 0 : 1;
// 	 }

//        distance[k]= ( hit0 | hit1 | hit2  )? factor*infactor*(hit0*snxt[0] + hit1*snxt[1] + hit2*snxt[2]) : infactor*TGeoShape::Big();
//      }
// }
