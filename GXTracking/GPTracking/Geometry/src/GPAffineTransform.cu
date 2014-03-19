#include "GPAffineTransform.h"

FQUALIFIER 
void GPAffineTransform_Constructor( GPAffineTransform *This )
{
  This->rxx = 1;
  This->rxy = 0;
  This->rxz = 0;

  This->ryx = 0;
  This->ryy = 1;
  This->ryz = 0;

  This->rzx = 0;
  This->rzy = 0;
  This->rzz = 1;

  This->tx = 0;
  This->ty = 0;
  This->tz = 0;
}

FQUALIFIER 
GPAffineTransform GPAffineTransform_create_id(void)
{
  GPAffineTransform t;
  GPAffineTransform_Constructor(&t);
  return t;
}

FQUALIFIER 
void GPAffineTransform_Set(GPAffineTransform *This,
			   GPAffineTransform rhs)
{
  //  GPAffineTransform_Constructor(This);
  This->rxx = rhs.rxx;
  This->ryy = rhs.ryy;
  This->rzz = rhs.rzz;
  This->rxy = rhs.rxy;
  This->rxz = rhs.rxz;
  This->ryx = rhs.ryx;
  This->ryz = rhs.ryz;
  This->rzx = rhs.rzx;
  This->rzy = rhs.rzy;
}

FQUALIFIER 
void GPAffineTransform_Vector(GPAffineTransform *This,
			      const GPThreeVector tlate)
{
  GPAffineTransform_Constructor(This);
  This->tx = tlate.x;
  This->ty = tlate.y;
  This->tz = tlate.z;
}

FQUALIFIER 
void GPAffineTransform_Matrix(GPAffineTransform *This,
			      const GPRotationMatrix rot)
{
  GPAffineTransform_Constructor(This);
  This->rxx = rot.rxx;
  This->ryy = rot.ryy;
  This->rzz = rot.rzz;
  This->rxy = rot.rxy;
  This->rxz = rot.rxz;
  This->ryx = rot.ryx;
  This->ryz = rot.ryz;
  This->rzx = rot.rzx;
  This->rzy = rot.rzy;
}

FQUALIFIER 
void GPAffineTransform_Constructor2( GPAffineTransform *This,
				     const GPRotationMatrix rot,
				     const GPThreeVector tlate )
{
  This->rxx = rot.rxx;
  This->rxy = rot.rxy;
  This->rxz = rot.rxz;

  This->ryx = rot.ryx;
  This->ryy = rot.ryy;
  This->ryz = rot.ryz;

  This->rzx = rot.rzx;
  This->rzy = rot.rzy;
  This->rzz = rot.rzz;

  This->tx = tlate.x;
  This->ty = tlate.y;
  This->tz = tlate.z;

}

FQUALIFIER 
void GPAffineTransform_Constructor3( GPAffineTransform *This, 
				     const GPRotationMatrix *rot,
				     const GPThreeVector tlate)
{
  if (rot) GPAffineTransform_Constructor2( This, *rot, tlate );
  else GPAffineTransform_Vector( This, tlate );

}

FQUALIFIER 
void GPAffineTransform_Elements( GPAffineTransform *This, 
		  const G4double prxx,const G4double prxy,const G4double prxz,
	          const G4double pryx,const G4double pryy,const G4double pryz,
                  const G4double przx,const G4double przy,const G4double przz,
                  const G4double ptx,const G4double pty,const G4double ptz)
{
  This->rxx = prxx;
  This->rxy = prxy;
  This->rxz = prxz;

  This->ryx = pryx;
  This->ryy = pryy;
  This->ryz = pryz;

  This->rzx = przx;
  This->rzy = przy;
  This->rzz = przz;

  This->tx = ptx;
  This->ty = pty;
  This->tz = ptz;
}

/*
GPAffineTransform_GPAffineTransform_operator * (const GPAffineTransform& tf) const
{
        return GPAffineTransform(
        rxx*tf.rxx+rxy*tf.ryx+rxz*tf.rzx,
        rxx*tf.rxy+rxy*tf.ryy+rxz*tf.rzy,
        rxx*tf.rxz+rxy*tf.ryz+rxz*tf.rzz,

        ryx*tf.rxx+ryy*tf.ryx+ryz*tf.rzx,
        ryx*tf.rxy+ryy*tf.ryy+ryz*tf.rzy,
        ryx*tf.rxz+ryy*tf.ryz+ryz*tf.rzz,

        rzx*tf.rxx+rzy*tf.ryx+rzz*tf.rzx,
        rzx*tf.rxy+rzy*tf.ryy+rzz*tf.rzy,
        rzx*tf.rxz+rzy*tf.ryz+rzz*tf.rzz,
        
        tx*tf.rxx+ty*tf.ryx+tz*tf.rzx+tf.tx,
        tx*tf.rxy+ty*tf.ryy+tz*tf.rzy+tf.ty,
        tx*tf.rxz+ty*tf.ryz+tz*tf.rzz+tf.tz);
}

GPAffineTransform& GPAffineTransform_operator *= (const GPAffineTransform& tf)
{
  // Use temporaries for `in place' compound transform computation

  G4double nrxx=rxx*tf.rxx+rxy*tf.ryx+rxz*tf.rzx;
  G4double nrxy=rxx*tf.rxy+rxy*tf.ryy+rxz*tf.rzy;
  G4double nrxz=rxx*tf.rxz+rxy*tf.ryz+rxz*tf.rzz;

  G4double nryx=ryx*tf.rxx+ryy*tf.ryx+ryz*tf.rzx;
  G4double nryy=ryx*tf.rxy+ryy*tf.ryy+ryz*tf.rzy;
  G4double nryz=ryx*tf.rxz+ryy*tf.ryz+ryz*tf.rzz;
  
  G4double nrzx=rzx*tf.rxx+rzy*tf.ryx+rzz*tf.rzx;
  G4double nrzy=rzx*tf.rxy+rzy*tf.ryy+rzz*tf.rzy;
  G4double nrzz=rzx*tf.rxz+rzy*tf.ryz+rzz*tf.rzz;
  
  G4double ntx=tx*tf.rxx+ty*tf.ryx+tz*tf.rzx+tf.tx;
  G4double nty=tx*tf.rxy+ty*tf.ryy+tz*tf.rzy+tf.ty;
  G4double ntz=tx*tf.rxz+ty*tf.ryz+tz*tf.rzz+tf.tz;
  
  tx=ntx; ty=nty; tz=ntz;
  rxx=nrxx; rxy=nrxy; rxz=nrxz;
  ryx=nryx; ryy=nryy; ryz=nryz;
  rzx=nrzx; rzy=nrzy; rzz=nrzz;
  
  return *this;
}
*/

FQUALIFIER 
GPAffineTransform GPAffineTransform_Product(GPAffineTransform *This,
					    const GPAffineTransform *tf1,
					    const GPAffineTransform *tf2)
{
  This->rxx=tf1->rxx*tf2->rxx + tf1->rxy*tf2->ryx + tf1->rxz*tf2->rzx;
  This->rxy=tf1->rxx*tf2->rxy + tf1->rxy*tf2->ryy + tf1->rxz*tf2->rzy;
  This->rxz=tf1->rxx*tf2->rxz + tf1->rxy*tf2->ryz + tf1->rxz*tf2->rzz;
  
  This->ryx=tf1->ryx*tf2->rxx + tf1->ryy*tf2->ryx + tf1->ryz*tf2->rzx;
  This->ryy=tf1->ryx*tf2->rxy + tf1->ryy*tf2->ryy + tf1->ryz*tf2->rzy;
  This->ryz=tf1->ryx*tf2->rxz + tf1->ryy*tf2->ryz + tf1->ryz*tf2->rzz;
  
  This->rzx=tf1->rzx*tf2->rxx + tf1->rzy*tf2->ryx + tf1->rzz*tf2->rzx;
  This->rzy=tf1->rzx*tf2->rxy + tf1->rzy*tf2->ryy + tf1->rzz*tf2->rzy;
  This->rzz=tf1->rzx*tf2->rxz + tf1->rzy*tf2->ryz + tf1->rzz*tf2->rzz;
  
  This->tx=tf1->tx*tf2->rxx + tf1->ty*tf2->ryx + tf1->tz*tf2->rzx   + tf2->tx;
  This->ty=tf1->tx*tf2->rxy + tf1->ty*tf2->ryy + tf1->tz*tf2->rzy   + tf2->ty;
  This->tz=tf1->tx*tf2->rxz + tf1->ty*tf2->ryz + tf1->tz*tf2->rzz   + tf2->tz; 
  
  return *This;
}

FQUALIFIER 
GPAffineTransform GPAffineTransform_InverseProduct( GPAffineTransform *This,
						    GPAffineTransform *tf1,
						    GPAffineTransform *tf2)
{
  G4double itf2tx = - tf2->tx*tf2->rxx - tf2->ty*tf2->rxy - tf2->tz*tf2->rxz;
  G4double itf2ty = - tf2->tx*tf2->ryx - tf2->ty*tf2->ryy - tf2->tz*tf2->ryz;
  G4double itf2tz = - tf2->tx*tf2->rzx - tf2->ty*tf2->rzy - tf2->tz*tf2->rzz;

  This->rxx=tf1->rxx*tf2->rxx+tf1->rxy*tf2->rxy+tf1->rxz*tf2->rxz;
  This->rxy=tf1->rxx*tf2->ryx+tf1->rxy*tf2->ryy+tf1->rxz*tf2->ryz;
  This->rxz=tf1->rxx*tf2->rzx+tf1->rxy*tf2->rzy+tf1->rxz*tf2->rzz;

  This->ryx=tf1->ryx*tf2->rxx+tf1->ryy*tf2->rxy+tf1->ryz*tf2->rxz;
  This->ryy=tf1->ryx*tf2->ryx+tf1->ryy*tf2->ryy+tf1->ryz*tf2->ryz;
  This->ryz=tf1->ryx*tf2->rzx+tf1->ryy*tf2->rzy+tf1->ryz*tf2->rzz;
  
  This->rzx=tf1->rzx*tf2->rxx+tf1->rzy*tf2->rxy+tf1->rzz*tf2->rxz;
  This->rzy=tf1->rzx*tf2->ryx+tf1->rzy*tf2->ryy+tf1->rzz*tf2->ryz;
  This->rzz=tf1->rzx*tf2->rzx+tf1->rzy*tf2->rzy+tf1->rzz*tf2->rzz;
        
  This->tx=tf1->tx*tf2->rxx+tf1->ty*tf2->rxy+tf1->tz*tf2->rxz+itf2tx;
  This->ty=tf1->tx*tf2->ryx+tf1->ty*tf2->ryy+tf1->tz*tf2->ryz+itf2ty;
  This->tz=tf1->tx*tf2->rzx+tf1->ty*tf2->rzy+tf1->tz*tf2->rzz+itf2tz;
  
  return *This;
}

FQUALIFIER
GPThreeVector GPAffineTransform_TransformPoint(const GPAffineTransform *This, 
					       const GPThreeVector vec)
{
  return GPThreeVector_create( vec.x*This->rxx + vec.y*This->ryx + vec.z*This->rzx   + This->tx,
			       vec.x*This->rxy + vec.y*This->ryy + vec.z*This->rzy   + This->ty,
			       vec.x*This->rxz + vec.y*This->ryz + vec.z*This->rzz   + This->tz  );
}

FQUALIFIER
GPThreeVector GPAffineTransform_TransformAxis(const GPAffineTransform *This,
					      const GPThreeVector axis)
{
  return GPThreeVector_create( axis.x*This->rxx + axis.y*This->ryx + axis.z*This->rzx,
			       axis.x*This->rxy + axis.y*This->ryy + axis.z*This->rzy,
			       axis.x*This->rxz + axis.y*This->ryz + axis.z*This->rzz  );
}

FQUALIFIER
GPThreeVector GPAffineTransform_ApplyPointTransform(const GPAffineTransform *This,
						    GPThreeVector vec)
{
  G4double x = vec.x*This->rxx + vec.y*This->ryx + vec.z*This->rzx    + This->tx;
  G4double y = vec.x*This->rxy + vec.y*This->ryy + vec.z*This->rzy    + This->ty;
  G4double z = vec.x*This->rxz + vec.y*This->ryz + vec.z*This->rzz    + This->tz;

  return GPThreeVector_create(x,y,z);
}

FQUALIFIER
GPThreeVector GPAffineTransform_ApplyAxisTransform(const GPAffineTransform *This,
						   GPThreeVector axis)
{
  G4double x = axis.x*This->rxx + axis.y*This->ryx + axis.z*This->rzx;
  G4double y = axis.x*This->rxy + axis.y*This->ryy + axis.z*This->rzy;
  G4double z = axis.x*This->rxz + axis.y*This->ryz + axis.z*This->rzz;

  return GPThreeVector_create(x,y,z);
}

FQUALIFIER
GPAffineTransform GPAffineTransform_Inverse(GPAffineTransform *This)
{
  GPAffineTransform aT;
  GPAffineTransform_Elements( &aT , 
			      This->rxx, This->ryx, This->rzx,
			      This->rxy, This->ryy, This->rzy,
			      This->rxz, This->ryz, This->rzz,
			      
			      -This->tx*This->rxx - This->ty*This->rxy - This->tz*This->rxz,
			      -This->tx*This->ryx - This->ty*This->ryy - This->tz*This->ryz,
			      -This->tx*This->rzx - This->ty*This->rzy - This->tz*This->rzz  );

  return aT;
}

FQUALIFIER
GPAffineTransform GPAffineTransform_Invert(GPAffineTransform *This)
{
  G4double v1 = -This->tx*This->rxx - This->ty*This->rxy - This->tz*This->rxz;
  G4double v2 = -This->tx*This->ryx - This->ty*This->ryy - This->tz*This->ryz;
  G4double v3 = -This->tx*This->rzx - This->ty*This->rzy - This->tz*This->rzz;

  This->tx=v1; This->ty=v2; This->tz=v3;

  G4double tmp1=This->ryx; This->ryx=This->rxy; This->rxy=tmp1;
  G4double tmp2=This->rzx; This->rzx=This->rxz; This->rxz=tmp2;
  G4double tmp3=This->rzy; This->rzy=This->ryz; This->ryz=tmp3;

  return *This;

}

/*
GPAffineTransform& GPAffineTransform_operator +=(const GPThreeVector& tlate)
{
        tx += tlate.x;
        ty += tlate.y;
        tz += tlate.z;

        return *this;
}

GPAffineTransform& GPAffineTransform_operator -=(const GPThreeVector& tlate)
{
        tx -= tlate.x;
        ty -= tlate.y;
        tz -= tlate.z;

        return *this;
}

G4bool GPAffineTransform_operator == (const GPAffineTransform& tf) const
{
        return (tx==tf.tx&&ty==tf.ty&&tz==tf.tz&&
                rxx==tf.rxx&&rxy==tf.rxy&&rxz==tf.rxz&&
                ryx==tf.ryx&&ryy==tf.ryy&&ryz==tf.ryz&&
                rzx==tf.rzx&&rzy==tf.rzy&&rzz==tf.rzz) ? true : false;
}
G4bool GPAffineTransform_operator != (const GPAffineTransform& tf) const
{
        return (tx!=tf.tx||ty!=tf.ty||tz!=tf.tz||
                rxx!=tf.rxx||rxy!=tf.rxy||rxz!=tf.rxz||
                ryx!=tf.ryx||ryy!=tf.ryy||ryz!=tf.ryz||
                rzx!=tf.rzx||rzy!=tf.rzy||rzz!=tf.rzz) ? true : false;
}

G4double GPAffineTransform_operator [] (const G4int n) const
{
        G4double v = 0.0;
        switch(n)
                {
                case 0:
                        v=rxx;
                        break;
                case 1:
                        v=rxy;
                        break;
                case 2:
                        v=rxz;
                        break;
                case 4:
                        v=ryx;
                        break;
                case 5:
                        v=ryy;
                        break;
                case 6:
                        v=ryz;
                        break;
                case 8:
                        v=rzx;
                        break;
                case 9:
                        v=rzy;
                        break;
                case 10:
                        v=rzz;
                        break;
                case 12:
                        v=tx;
                        break;
                case 13:
                        v=ty;
                        break;
                case 14:
                        v=tz;
                        break;
                case 3:
                case 7:
                case 11:
                        break;
                case 15:
                        v=1.0;
                        break;
                }
        return v;
}
*/

FQUALIFIER
G4bool GPAffineTransform_IsRotated( GPAffineTransform *This )
{
  return (This->rxx==1.0 && This->ryy==1.0 && This->rzz==1.0) ? false : true;
}

FQUALIFIER 
G4bool GPAffineTransform_IsTranslated(GPAffineTransform *This)
{
  return (This->tx || This->ty || This->tz) ? true:false;
}

FQUALIFIER 
GPRotationMatrix GPAffineTransform_NetRotation( GPAffineTransform *This ) 
{
  GPRotationMatrix m = GPRotationMatrix_create(1,0,0,
					       0,1,0,
					       0,0,1);
  return GPRotationMatrix_rotateAxes(&m,
				     GPThreeVector_create(This->rxx,This->ryx,This->rzx),
				     GPThreeVector_create(This->rxy,This->ryy,This->rzy),
				     GPThreeVector_create(This->rxz,This->ryz,This->rzz));
}

FQUALIFIER
GPThreeVector GPAffineTransform_NetTranslation(GPAffineTransform *This)
{
  return GPThreeVector_create(This->tx,This->ty,This->tz);
}

FQUALIFIER 
void GPAffineTransform_SetNetRotation(GPAffineTransform *This,
				      const GPRotationMatrix rot)
{
  This->rxx=rot.rxx;
  This->rxy=rot.rxy;
  This->rxz=rot.rxz;
  This->ryx=rot.ryx;
  This->ryy=rot.ryy;
  This->ryz=rot.ryz;
  This->rzx=rot.rzx;
  This->rzy=rot.rzy;
  This->rzz=rot.rzz;
}

FQUALIFIER
void GPAffineTransform_SetNetTranslation(GPAffineTransform *This,
					 const GPThreeVector tlate)
{
  This->tx=GPThreeVector_x(tlate);
  This->ty=GPThreeVector_y(tlate);
  This->tz=GPThreeVector_z(tlate);
}
