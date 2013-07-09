#include "Vc/vector.h"
#include <iostream>
// #include "Vc/IO"

/*
template<typename T>
static std::ostream &operator<<(std::ostream &out, const Vc::SSE::Vector<T> &v)
{
  // out << AnsiColor::green << "[";
  out << "[" << v[0];
    for (int i = 1; i < v.Size; ++i) {
        out << ", " << v[i];
    }
    out << "]";
    return out;
}
*/
typedef Vc::double_v vd; // short for vector double

// // this is the actual kernel doing the computation with possibility of early return
void TGeoBBox_v::DistFromOutsideS_v(StructOfCoord const & point, StructOfCoord const & dir, 
                          double dx, double dy, double dz, const double *origin, const double * stepmax, double * distance, int np)
 {
   vd zero=0.;
   vd vcorigin[3]={origin[0],origin[1],origin[2]};
   int tailsize =np % Vc::double_v::Size; 
   for(volatile unsigned int i = 0; i < np-tailsize; i+= Vc::double_v::Size) 
    {
      // this is a lot of memory movement ( maybe a cast would be better )?
      vd x(&point.x[i]);
      vd y(&point.y[i]);
      vd z(&point.z[i]);
      vd dirx(&dir.x[i]);
      vd diry(&dir.y[i]);
      vd dirz(&dir.z[i]);
      vd s(&stepmax[i]);
      vd dist(0.);
      DistFromOutsideS_v4( x, y, z, dirx, diry, dirz, dx, dy, dz, vcorigin, s, dist );
      dist.store(&distance[i]);
   }
   // do the tail part for the moment, we just call the old static version
   for(unsigned int i = 0; i < tailsize; ++i)
     {
       double p[3]={point.x[np-tailsize+i], point.y[np-tailsize+i], point.z[np-tailsize+i]};
       double d[3]={dir.x[np-tailsize+i], dir.y[np-tailsize+i], dir.z[np-tailsize+i]};
       distance[np-tailsize+i]=TGeoBBox_v::DistFromOutsideS(p, d, dx, dy, dz, origin, stepmax[np-tailsize+i] );
     }
 }


 // this is the actual kernel doing the computation with possibility of early return
void TGeoBBox_v::DistFromOutsideS_v4( Vc::double_v const & x , Vc::double_v const & y, Vc::double_v const & z, Vc::double_v const & dirx, Vc::double_v const & diry, Vc::double_v const & dirz, double dx, double dy, double dz, Vc::double_v const origin[3], Vc::double_v const & stepmax, Vc::double_v & distance )
 {
   Vc::double_v in=Vc::double_v(1.);
   Vc::double_v saf[3];
   Vc::double_v newpt[3];
   Vc::double_v tiny(1e-20);
   Vc::double_v big(1e30);
   Vc::double_v faraway=Vc::double_v(0.); // initializing all components to zero
   Vc::double_v par[3]={dx,dy,dz}; // very convenient

   newpt[0] = x-origin[0];
   saf[0] = Vc::abs(newpt[0])-par[0]; // can we profit here from abs function in array form?
   newpt[1] = y-origin[1];
   saf[1] = Vc::abs(newpt[1])-par[1];
   newpt[2] = z-origin[2];
   saf[2] = Vc::abs(newpt[2])-par[2];
   faraway(saf[0]>=stepmax | saf[1]>=stepmax | saf[2]>=stepmax)=1; 
   in(saf[0]<0. & saf[1]<0. & saf[2]<0.)=0;
   distance=big;

   if( faraway > Vc::double_v(0.) )
     {
       return; // return big
     }
  
   // proceed to analysis of hits
   Vc::double_v snxt[3];
   Vc::double_v hit0=Vc::double_v(0.);
   snxt[0] = saf[0]/Vc::abs(dirx+tiny); // distance to y-z face
   Vc::double_v coord1=newpt[1]+snxt[0]*diry; // calculate new y and z coordinate
   Vc::double_v coord2=newpt[2]+snxt[0]*dirz;
   hit0( saf[0] > 0 & newpt[0]*dirx < 0 & ( Vc::abs(coord1)<= par[1] & Vc::abs(coord2)<= par[2] ) ) = 1; // if out and right direction
   
   Vc::double_v hit1=Vc::double_v(0.);
   snxt[1] = saf[1]/Vc::abs(diry+tiny); // distance to x-z face
   coord1=newpt[0]+snxt[1]*dirx; // calculate new x and z coordinate
   coord2=newpt[2]+snxt[1]*dirz;
   hit1( saf[1] > 0 & newpt[1]*diry < 0 & ( Vc::abs(coord1)<= par[0] & Vc::abs(coord2)<= par[2] ) ) = 1; // if out and right direction

   Vc::double_v hit2=Vc::double_v(0.);
   snxt[2] = saf[2]/Vc::abs(dirz+tiny); // distance to x-y face
   coord1=newpt[0]+snxt[2]*dirx; // calculate new x and y coordinate
   coord2=newpt[1]+snxt[2]*diry;
   hit2( saf[2] > 0 & newpt[2]*dirz < 0 & ( Vc::abs(coord1)<= par[0] & Vc::abs(coord2)<= par[1] ) ) = 1; // if out and right direction

   distance( hit0>0 | hit1>0 | hit2>0 )=(hit0*snxt[0] + hit1*snxt[1] + hit2*snxt[2]);
   distance=in*distance;
   return;
}

// VC implementation of Contains function
void TGeoBBox_v::Contains_v(const StructOfCoord & pointi, Bool_t * isin, Int_t np) const
{
   vd vcorigin[3]={fOrigin[0],fOrigin[1],fOrigin[2]};
   vd vfDX(fDX);
   vd vfDY(fDY);
   vd vfDZ(fDZ);

   int tailsize =np % Vc::double_v::Size; 
   for(unsigned int i = 0; i < np-tailsize; i+= Vc::double_v::Size) 
    {
      vd x(&pointi.x[i]); // will copy a certain number of x's into vc vector x;
      vd y(&pointi.y[i]); 
      vd z(&pointi.z[i]);
      vd particleinside(0.); // we initialize it to zero meaning outside 
      x=Vc::abs(x-vcorigin[0]);
      y=Vc::abs(y-vcorigin[1]);
      z=Vc::abs(z-vcorigin[2]);
      particleinside( (x < vfDZ) & (y < vfDY) & (z<vfDZ) ) = 1.; // we set it to one if inside
      
      for(unsigned int j=0; j<Vc::double_v::Size; ++j)
	{
	  isin[i+j] = (bool) particleinside[j];
	}
   }
   // do the tail part for the moment, we just call the old static version
   for(unsigned int i = 0; i < tailsize; ++i)
     {
       Double_t xx, yy, zz;
       xx= pointi.x[np-tailsize+i] - fOrigin[0];
       yy= pointi.y[np-tailsize+i] - fOrigin[1];
       zz= pointi.z[np-tailsize+i] - fOrigin[2];
       isin[i]=(TMath::Abs(xx)<fDX) & (TMath::Abs(yy)<fDY) & (TMath::Abs(zz)<fDZ); 
     }
}



/*
void testKernel()
{
  vd x,y,z,dirx,diry,dirz;
  vd zero=0.;
  vd origin[3]={zero,zero,zero};
  x[0]=1000;x[1]=2000;
  y[0]=300;y[1]=49;
  z[0]=4;z[1]=49;

  dirx[0]=1;dirx[1]=-1;
  diry[0]=0;diry[1]=0;
  dirz[0]=0;dirz[1]=0;

  vd s;
  s[0]=10000;s[1]=10000;
  vd distance;

  DistFromOutside_VECEXT_P4( x, y, z, dirx, diry, dirz, 100, 50, 50, origin, s, distance );
  std::cerr << distance[0] << " " << distance[1] << std::endl;
}
*/

//int main()
//{
//  testKernel();
//}
