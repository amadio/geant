#include <iostream>
#include "TGeoCone_v.h"
#include "TMath.h"
#ifdef VEC_EXTENSIONS
#include "Vc/vector.h"
//#include <Vc/double_v>
typedef Vc::double_v vd;
#endif

Vc::int_v LocMin(Vc::int_v n, const Vc::double_v *a)
{
    // Return index of array with the minimum element.
    // If more than one element is minimum returns first found.
    
    if  (n <= 0 || !a) return -1;
    Vc::double_v xmin = a[0];
    Vc::int_v loc = 0;
    for  (int i = 1; i < n; i++) {
        if (xmin > a[i])  {
            xmin = a[i];
            loc = i;
        }
    }
    return loc;
}

/*
 TGEOSHAPE::SAFETYSEG
 
 Double_t TGeoShape::SafetySeg(Double_t r, Double_t z, Double_t r1, Double_t z1, Double_t r2, Double_t z2, Bool_t outer)
 {
 
 Double_t crossp = (z2-z1)*(r-r1)-(z-z1)*(r2-r1);
 crossp *= (outer) ? 1. : -1.;
 if (crossp < 0) {
 if (((z-z1)*(z2-z)) > 0) return 0;
 return TGeoShape::Big();
 }   
 // Compute (1,P) dot (1,2)
 Double_t c1 = (z-z1)*(z2-z1)+(r-r1)*(r2-r1);
 // Negative c1 means point (1) is closest
 if (c1<1.E-10) return TMath::Sqrt((r-r1)*(r-r1)+(z-z1)*(z-z1));
 // Compute (2,P) dot (1,2)
 Double_t c2 = (z-z2)*(z2-z1)+(r-r2)*(r2-r1);
 // Positive c2 means point (2) is closest
 if (c2>-1.E-10) return TMath::Sqrt((r-r2)*(r-r2)+(z-z2)*(z-z2));
 // The closest point is between (1) and (2)
 c2 = (z2-z1)*(z2-z1)+(r2-r1)*(r2-r1);
 // projected length factor with respect to (1,2) length
 Double_t alpha = c1/c2;
 Double_t rp = r1 + alpha*(r2-r1);
 Double_t zp = z1 + alpha*(z2-z1);
 return TMath::Sqrt((r-rp)*(r-rp)+(z-zp)*(z-zp));
 }
 
 */

void geoShapeSafetySeg(Vc::double_v r,Vc::double_v z, Vc::double_v r1, Vc::double_v z1, Vc::double_v r2, Vc::double_v z2, Bool_t outer, Vc::double_v & return_val)
{
    // Compute distance from point of coordinates (r,z) to segment (r1,z1):(r2,z2)
    Vc::double_v crossp = (z2-z1)*(r-r1)-(z-z1)*(r2-r1);
    
    /*for(int k=0;k<Vc::double_v::Size; k++)
    {
        std::cout<<"crossp ["<<k<<"] "<<crossp[k]<<std::endl;
        
    }
    */

    crossp *= (outer) ? 1. : -1.;
    
    
    /*for(int k=0;k<Vc::double_v::Size; k++)
    {
        std::cout<<"crossp ["<<k<<"] "<<crossp[k]<<std::endl;
        std::cout<<"((z-z1)*(z2-z)) ["<<k<<"] "<<(z[k]-z1[k])*(z2[k]-z[k])<<std::endl;
        std::cout<<"exit0:0 "<<std::endl;
        std::cout<<"exit1: "<<TGeoShape::Big()<<std::endl;
        std::cout<<"exit2: "<< TMath::Sqrt((r[k]-r1[k])*(r[k]-r1[k])+(z[k]-z1[k])*(z[k]-z1[k]))<<std::endl;
        std::cout<<"exit3: "<< TMath::Sqrt((r[k]-r2[k])*(r[k]-r2[k])+(z[k]-z2[k])*(z[k]-z2[k]))<<std::endl;
        //std::cout<<"exit4: "<< TMath::Sqrt((r[k]-rp[k])*(r[k]-rp[k])+(z[k]-zp[k])*(z[k]-zp[k]))<<"\n"<<std::endl;
    }*/
    // Positive crossp means point on the requested side of the (1,2) segment
    if((crossp<0))
    {
        if(((z-z1)*(z2-z)) > 0) {return_val=0.;
            return;
        }
        else
        {
            return_val = TGeoShape::Big();
            return;
        }
        
    }
    else{
        // Compute (1,P) dot (1,2)
        Vc::double_v c1 = (z-z1)*(z2-z1)+(r-r1)*(r2-r1);
        // Negative c1 means point (1) is closest
        if(c1<1.E-10){
            return_val=Vc::sqrt((r-r1)*(r-r1)+(z-z1)*(z-z1));
            return;
            
        }else{
            // Compute (2,P) dot (1,2)
            Vc::double_v c2 = (z-z2)*(z2-z1)+(r-r2)*(r2-r1);
            // Positive c2 means point (2) is closest
            if(c2>-1.E-10){
                return_val=Vc::sqrt((r-r2)*(r-r2)+(z-z2)*(z-z2));
                return;
            }
            else{
                // The closest point is between (1) and (2)
                c2 = (z2-z1)*(z2-z1)+(r2-r1)*(r2-r1);
                // projected length factor with respect to (1,2) length
                Vc::double_v alpha = c1/c2;
                Vc::double_v rp = r1 + alpha*(r2-r1);
                Vc::double_v zp = z1 + alpha*(z2-z1);
                return_val=Vc::sqrt((r-rp)*(r-rp)+(z-zp)*(z-zp));
                return;
            }
        }
    }
}

#ifndef VEC_EXTENSIONS 
//_____________________________________________________________________________                                                                                                     
void TGeoCone_v::Contains_v(const StructOfCoord  & pointi, Bool_t * isin, Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      isin[i] = TGeoCone::Contains(point);
    }
}
#else
//_____________________________________________________________________________                                                                                                     
void TGeoCone_v::Contains_v(const StructOfCoord  & pointi, Bool_t * isin, Int_t np) const
{
    static vd vfDz(fDz);
    static vd vfRmin1(fRmin1);
    static vd vfRmin2(fRmin2);
    static vd vfRmax1(fRmax1);
    static vd vfRmax2(fRmax2);
    
    Vc::double_m particleinside;
    
    int tailsize =np % Vc::double_v::Size;
    for(unsigned int i = 0; i < np-tailsize; i+= Vc::double_v::Size)
    {
        vd x(&pointi.x[i]); // will copy a certain number of x's into vc vector x;
        vd y(&pointi.y[i]);
        vd z(&pointi.z[i]);
        
        vd me=Vc::abs(z);
        Vc::double_m c1 = (me > vfDz);
        //  if( c1 ) continue;
        
        vd r2=x*x+y*y;
        vd rl=0.5*(vfRmin2*(z+vfDz)+vfRmin1*(vfDz-z))/vfDz;
        vd rh= 0.5*(vfRmax2*(z+vfDz)+vfRmax1*(vfDz-z))/vfDz;
        Vc::double_m c2 = (r2<rl*rl);
        // if( c2 ) continue;
        
        Vc::double_m c3 = (r2>rh*rh);
        // if( c3 ) continue;
        
        particleinside = !c1 && !c2 && !c3;
        for(unsigned int j=0;j<Vc::double_v::Size;++j) isin[i+j]= particleinside[j];
    }
    
    // do the tail part for the moment, we just call the old static version
    for(unsigned int i = 0; i < tailsize; ++i)
    {
        Double_t xx, yy, zz;
        xx= pointi.x[np-tailsize+i];
        yy= pointi.y[np-tailsize+i];
        zz= pointi.z[np-tailsize+i];
        
        Double_t rr2 = xx*xx+yy*yy;
        Double_t rrl = 0.5*(fRmin2*(zz+fDz)+fRmin1*(fDz-zz))/fDz;
        Double_t rrh = 0.5*(fRmax2*(zz+fDz)+fRmax1*(fDz-zz))/fDz;
        isin[np-tailsize+i]=( (TMath::Abs(zz) <= fDz) & (rr2>=rrl*rrl) & (rr2<=rrh*rrh) );
    }
}
#endif

#ifndef VEC_EXTENSIONS 
//_____________________________________________________________________________
void TGeoCone_v::Safety_v(const StructOfCoord & pointi, Bool_t in, Double_t * safety, Int_t np ) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      safety[i] = TGeoCone::Safety(point, in);
    }
}

#else
/*
 ---------------
 SAFETY
 ---------------
 
 Double_t TGeoCone::Safety(const Double_t *point, Bool_t in) const
 {
 // computes the closest distance from given point to this shape, according
 // to option. The matching point on the shape is stored in spoint.
 Double_t saf[4];
 Double_t r=TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
 saf[0] = TGeoShape::SafetySeg(r,point[2], fRmin1, -fDz, fRmax1, -fDz, !in);
 saf[1] = TGeoShape::SafetySeg(r,point[2], fRmax2, fDz, fRmin2, fDz, !in);
 saf[2] = TGeoShape::SafetySeg(r,point[2], fRmin2, fDz, fRmin1, -fDz, !in);
 saf[3] = TGeoShape::SafetySeg(r,point[2], fRmax1, -fDz, fRmax2, fDz, !in);
 return saf[TMath::LocMin(4,saf)];
 } 
 */

//_____________________________________________________________________________
void TGeoCone_v::Safety_v(const StructOfCoord & pointi, Bool_t in, Double_t * safety, Int_t np ) const
{
    
    //Vc implementation
    vd const vfRmin1(fRmin1);
    static vd vfRmin2(fRmin2);
    static vd vfRmax1(fRmax1);
    static vd vfRmax2(fRmax2);
    static vd vfDz(fDz);
    
    Vc::double_m particleinside;
    int vectorsize=Vc::double_v::Size;
    int tailsize =np % vectorsize;
    
    for(unsigned int i = 0; i < np-tailsize; i+= vectorsize)
    {
        vd x(&pointi.x[i]); // will copy a certain number of x's into vc vector x;
        vd y(&pointi.y[i]);
        vd z(&pointi.z[i]);
        
        /*Trivial and not really vectorized version*/
        Double_t cr[vectorsize];
        Double_t csaf[4];
    
        for(int k=0;k<Vc::double_v::Size; k++)
        {
            cr[k]=TMath::Sqrt(pointi.x[i+k]*pointi.x[i+k]+pointi.y[i+k]*pointi.y[i+k]);
            csaf[0] = TGeoShape::SafetySeg(cr[k],pointi.z[i+k], fRmin1, -fDz, fRmax1, -fDz, !in);
            csaf[1] = TGeoShape::SafetySeg(cr[k],pointi.z[i+k], fRmax2, fDz, fRmin2, fDz, !in);
            csaf[2] = TGeoShape::SafetySeg(cr[k],pointi.z[i+k], fRmin2, fDz, fRmin1, -fDz, !in);
            csaf[3] = TGeoShape::SafetySeg(cr[k],pointi.z[i+k], fRmax1, -fDz, fRmax2, fDz, !in);
            
            safety[i+k]=csaf[TMath::LocMin(4,csaf)]; 
            
        }
    }
    
    // do the tail part for the moment, we just call the old static version
   for(unsigned int i = 0; i < tailsize; ++i)
    {
    
        Double_t xx, yy, zz;
        xx= pointi.x[np-tailsize+i];
        yy= pointi.y[np-tailsize+i];
        zz= pointi.z[np-tailsize+i];
        Double_t xsaf[4];
        
        Double_t xr=TMath::Sqrt(xx*xx+yy*yy);
        xsaf[0] = TGeoShape::SafetySeg(xr,zz, fRmin1, -fDz, fRmax1, -fDz, !in);
        xsaf[1] = TGeoShape::SafetySeg(xr,zz, fRmax2, fDz, fRmin2, fDz, !in);
        xsaf[2] = TGeoShape::SafetySeg(xr,zz, fRmin2, fDz, fRmin1, -fDz, !in);
        xsaf[3] = TGeoShape::SafetySeg(xr,zz, fRmax1, -fDz, fRmax2, fDz, !in);
        safety[np-tailsize+i]=xsaf[TMath::LocMin(4,xsaf)];
    }
}

#endif


//_____________________________________________________________________________                                                                                                     
void TGeoCone_v::DistFromInside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/,const Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoCone::DistFromInside(point, dir, 3, step[i], 0);
    }
}

//_____________________________________________________________________________                                                                                                    

void TGeoCone_v::DistFromOutside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/, const  Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
    for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoCone::DistFromOutside(point, dir, 3, step[i], 0);
    }

}

#ifndef VEC_EXTENSIONS
//_____________________________________________________________________________                                                                                                     
void TGeoConeSeg_v::Contains_v(const StructOfCoord  & pointi, Bool_t * isin, Int_t np) const
{
    for(unsigned int i = 0; i < np; i++)
    {
        
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      isin[i] = TGeoConeSeg::Contains(point);
    }
}
#else

//_____________________________________________________________________________
void TGeoConeSeg_v::Contains_v4( Vc::double_v const & x , Vc::double_v const & y, Vc::double_v const & z,  Vc::double_m & c1)
const{
    vd vfDz(fDz);
    vd vfRmin1(fRmin1);
    vd vfRmin2(fRmin2);
    vd vfRmax1(fRmax1);
    vd vfRmax2(fRmax2);
    
    vd absZ=Vc::abs(z);
    vd r2(x*x+y*y);
    vd rl(0.5*(vfRmin2*(z+vfDz)+vfRmin1*(vfDz-z))/vfDz);
    vd rh(0.5*(vfRmax2*(z+vfDz)+vfRmax1*(vfDz-z))/vfDz);
    c1=((absZ<=vfDz) && (r2>=rl*rl) && (r2<=rh*rh));
}

//_____________________________________________________________________________
void TGeoConeSeg_v::Contains_v(const StructOfCoord  & pointi, Bool_t * isin, Int_t np) const
{
    static vd vfPhi1(fPhi1);
    static vd vfPhi2(fPhi2);
    int vectorsize=Vc::double_v::Size;
    int tailsize =np % vectorsize;
    Vc::double_m particleinside;
    
    for(unsigned int i = 0; i < np-tailsize; i+= vectorsize)
    {
        vd x(&pointi.x[i]); // will copy a certain number of x's into vc vector x
        vd y(&pointi.y[i]);
        vd z(&pointi.z[i]);
        vd dphi(vfPhi2 - vfPhi1);
        vd phi = Vc::atan2(y, x) * 57.295780181884765625f; // 180/phi
        //vd phi = Vc::atan2(y, x) * TMath::RadToDeg();
        
        phi(phi<0.f)+=360.f;
        //if (phi < 0 ) phi+=360.;
        vd ddp = phi-vfPhi1;
        
        ddp(ddp<0.f)+=360.f;
        //if (ddp < 0) ddp+=360.;
        Vc::double_m c1;
        TGeoConeSeg_v::Contains_v4(x,y,z,c1);
        
        Vc::double_m c2 = (dphi>=360.);
        Vc::double_m c3 = (ddp>dphi) ;
        particleinside = c1 && (c2||(!c2 && !c3));
        //particleinside = !c1;
        for(unsigned int j=0; j<vectorsize; ++j)
            isin[i+j] = (bool) particleinside[j];
    }
    // do the tail part for the moment, we just call the old static version
    for(unsigned int i = 0; i < tailsize; ++i)
    {
        Double_t point[3];
        
        point[0]= pointi.x[np-tailsize+i];
        point[1]= pointi.y[np-tailsize+i];
        point[2]= pointi.z[np-tailsize+i];
        
        Double_t t_dphi = fPhi2 - fPhi1;
        Double_t t_phi = TMath::ATan2(point[1], point[0]) * TMath::RadToDeg();
        if (t_phi < 0 ) t_phi+=360.;
        Double_t t_ddp = t_phi-fPhi1;
        if (t_ddp < 0) t_ddp+=360.;
        //isin[i]=( (TGeoCone::Contains(point) & (t_ddp<=t_dphi)) | (t_dphi>=360.));
        isin[np-tailsize+i]=( (TGeoCone::Contains(point) & (t_ddp<=t_dphi)) | (t_dphi>=360.));
        
    }
}
#endif

#ifndef VEC_EXTENSIONS

//_____________________________________________________________________________
void TGeoConeSeg_v::Safety_v(const StructOfCoord & pointi, Bool_t in, Double_t * safety, Int_t np ) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      safety[i] = TGeoConeSeg::Safety(point, in);
    }

}

#else
//_____________________________________________________________________________   
void TGeoConeSeg_v::Safety_v4(Vc::double_v const & x , Vc::double_v const & y, Vc::double_v const & z , Bool_t in, Vc::double_v &safety) const
{
    
    //Vc implementation
    static vd vfRmin1(fRmin1);
    static vd vfRmin2(fRmin2);
    static vd vfRmax1(fRmax1);
    static vd vfRmax2(fRmax2);
    static vd vfDz(fDz);
    
    int vectorsize=Vc::double_v::Size;
    Double_t cr[vectorsize];
    Double_t csaf[4];
        
    for(int k=0;k<vectorsize; k++)
    {
        cr[k]=TMath::Sqrt(x[k]*x[k]+y[k]*y[k]);
        csaf[0] = TGeoShape::SafetySeg(cr[k],z[k], fRmin1, -fDz, fRmax1, -fDz, !in);
        csaf[1] = TGeoShape::SafetySeg(cr[k],z[k], fRmax2, fDz, fRmin2, fDz, !in);
        csaf[2] = TGeoShape::SafetySeg(cr[k],z[k], fRmin2, fDz, fRmin1, -fDz, !in);
        csaf[3] = TGeoShape::SafetySeg(cr[k],z[k], fRmax1, -fDz, fRmax2, fDz, !in);
        safety[k]=csaf[TMath::LocMin(4,csaf)];
    }

}

/*
 
 ---------------
 SAFETY
 ---------------
 
 Double_t TGeoConeSeg::Safety(const Double_t *point, Bool_t in) const
 {
 // computes the closest distance from given point to this shape, according
 // to option. The matching point on the shape is stored in spoint.
 
 Double_t safe = TGeoCone::Safety(point,in);
 if ((fPhi2-fPhi1)>=360.) return safe; //c1
 Double_t safphi = TGeoShape::SafetyPhi(point, in, fPhi1, fPhi2);
 if (in) return TMath::Min(safe, safphi); //c2&& safeIsmin(c3) 
 if (safe>1.E10) return safphi; //!c2&& safe<1.E10 (c4)&& !safeIsMin (!c3) 
 return TMath::Max(safe, safphi);
 
 }
 
 */

//_____________________________________________________________________________   
void TGeoConeSeg_v::Safety_v(const StructOfCoord & pointi, Bool_t in, Double_t * safety, Int_t np ) const
{
    //Vc implementation
    static vd vfPhi1(fPhi1);
    static vd vfPhi2(fPhi2);
    int vectorsize=Vc::double_v::Size;
    int tailsize =np % vectorsize;
    
    for(unsigned int i = 0; i < np-tailsize; i+= vectorsize)
    {
        vd x(&pointi.x[i]); // will copy a certain number of x's into vc vector x;
        vd y(&pointi.y[i]);
        vd z(&pointi.z[i]);
        
        /*Trivial and not really vectorized version*/
        Vc::double_v vsafe;
        Vc::double_v vresult;
        
        TGeoConeSeg_v::Safety_v4(x,y,z,in,vsafe);
        Vc::double_m c1((vfPhi2-vfPhi1)>=360.); //if c1 return vsafe
        if (c1) vresult=vsafe; //if c1 is valid is better to exit here
        else{
            Vc::double_v vsafphi;
            for(int k=0;k<Vc::double_v::Size; k++)
            {
                Double_t point[3]={pointi.x[i+k], pointi.y[i+k], pointi.z[i+k]}; 
                vsafphi[k] = TGeoShape::SafetyPhi(point, in, fPhi1, fPhi2); //just to try if it works
            }
            Vc::double_m c2(in); //if c2 return min(vsafe, vsafphi)
            Vc::double_m c3(vsafe<=vsafphi);
            Vc::double_m c4(vsafe>1.E10); //if c4 return vsafphi;
            //else return max(vsafe, vsafphi);
            vresult = ((c2&c3)|((!c2)&(!c3)&(!c4))) ? vsafe  : vsafphi;
        }
        for(unsigned int k=0; k<vectorsize; ++k)
            safety[i+k] = vresult[k];
     
    }
    
    // do the tail part for the moment, we just call the old static version
    for(unsigned int i = 0; i < tailsize; ++i)
    {
                
        Double_t xx, yy, zz;
        xx= pointi.x[np-tailsize+i];
        yy= pointi.y[np-tailsize+i];
        zz= pointi.z[np-tailsize+i];
        Double_t xpoint[3]={xx, yy,zz};
        Double_t xsafe = TGeoCone::Safety(xpoint,in);
        if ((fPhi2-fPhi1)>=360.) safety[np-tailsize+i]=xsafe; 
        else
        {
            Double_t xsafphi = TGeoShape::SafetyPhi(xpoint, in, fPhi1, fPhi2);
            
            /*if(c1 || (!c1&& ((c2&&c3) || (!c2&&c4&&!c3) )))
            where 
            c2 in
            c3 (safeIsmin)
            c4 (safe>1.E10)*/
            //if ((in && xsafe<xsafphi) || (!in && (xsafe<1.E10) && (xsafe>xsafphi)) )
            
            if ( (in && xsafe<=xsafphi) || ((!in) && (xsafe<1.E10) && (xsafe>xsafphi)) ) safety[np-tailsize+i]=xsafe;
                else safety[np-tailsize+i]=xsafphi;
        }
    }
}
#endif
//_____________________________________________________________________________                                                                                                     
void TGeoConeSeg_v::DistFromInside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/,const Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoConeSeg::DistFromInside(point, dir, 3, step[i], 0);
    }
}

//_____________________________________________________________________________                                                                                                    

void TGeoConeSeg_v::DistFromOutside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/, const  Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
    for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoConeSeg::DistFromOutside(point, dir, 3, step[i], 0);
    }

}
