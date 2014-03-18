//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4MagneticField.cc,v 1.3 2006-06-29 18:24:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------

#include "GPMagneticField.h"

FQUALIFIER 
void GPMagneticField_Constructor(GPMagneticField *This,
				 GPFieldMap *fieldMapArray)
{
  This->fGravityActive = false;
  This->fieldMap = fieldMapArray;

  //parametric b-field
  const G4double CMS_Bp_3_8T[9] = {4.24326,15.0201,3.81492,0.0178712,0.000656527,
                               2.45818,0.00778695,2.12500,1.77436};

  for (G4int i=0; i<9; i++) This->prm[i] = CMS_Bp_3_8T[i];

  G4double ap2=4*This->prm[0]*This->prm[0]/(This->prm[1]*This->prm[1]);  

  This->hb0=0.5*This->prm[2]*sqrt(1.0+ap2);
  This->hlova=1/sqrt(ap2);
  This->ainv=2*This->hlova/This->prm[1];
  This->coeff=1/(This->prm[3]*This->prm[3]);

}

FQUALIFIER 
void GPMagneticField_Constructor_Parametric(GPMagneticField *This)
{
  
  This->fGravityActive = false;
  This->fieldMap = NULL;

  //parametric b-field
  const G4double CMS_Bp_3_8T[9] = {4.24326,15.0201,3.81492,0.0178712,0.000656527,
                               2.45818,0.00778695,2.12500,1.77436};

  for (G4int i=0; i<9; i++) This->prm[i] = CMS_Bp_3_8T[i];

  G4double ap2=4*This->prm[0]*This->prm[0]/(This->prm[1]*This->prm[1]);  

  This->hb0=0.5*This->prm[2]*sqrt(1.0+ap2);
  This->hlova=1/sqrt(ap2);
  This->ainv=2*This->hlova/This->prm[1];
  This->coeff=1/(This->prm[3]*This->prm[3]);
}

FQUALIFIER 
void  GPMagneticField_GetFieldValue( GPMagneticField *This,
				     const  G4double point[4],
				     G4double *bField ) 
{
  //default values
  bField[0] = 0.; 
  bField[1] = 0.;
  bField[2] = 0.;

  //  double rho = std::sqrt(point[0]*point[0]+point[1]*point[1]); 

  //magnetic field map array 
  //  const int nbinR = 900+1;
  const int nbinZ = 2*1600+1;
  const int noffZ = 1600;

  //magnetic field unit
  const G4double tesla = 0.001;

  G4double x        = point[0];        
  G4double y        = point[1];    

  //conversion: [mm] to [cm]
  G4double rho      = 0.1*sqrt(x*x+y*y);
  G4double z        = 0.1*point[2]; 

  if(abs(z)>1600.0 || rho > 900.0) {
    bField[0] = bField[1] = bField[2] = 0.0;
  }
  else {
    int ir = int(rho);
    G4double dr = rho - int(rho);
  
    int iz = int(z + noffZ);
    G4double dz = z + noffZ - iz;
    
    int index = iz+ir*nbinZ;
    
    G4double Bz_lb = This->fieldMap[index].Bz;
    G4double Bz_ub = This->fieldMap[index+nbinZ].Bz;
    
    bField[2] = tesla*(Bz_lb + (Bz_ub-Bz_lb)*dz);  
    
    if(rho>0) {
      G4double Br_lb = This->fieldMap[index].Br;
      G4double Br_ub = This->fieldMap[index+1].Br;
      bField[0] = tesla*(Br_lb + (Br_ub-Br_lb)*dr)*x/rho;
      bField[1] = tesla*(Br_lb + (Br_ub-Br_lb)*dr)*y/rho;
    }
  }
  //  printf("show VolumeBased Bxyz = %f %f %f \n",bField[0],bField[1],bField[2]);
  //use cuda texture memory for fieldMap

}

FQUALIFIER 
void  GPMagneticField_GetFieldValue_Parametric( GPMagneticField *This,
						const  G4double point[4],
						G4double *bField ) 
{

  //default values
  bField[0] = 0.; 
  bField[1] = 0.;
  bField[2] = 0.;

  //convert to m (cms magnetic field parameterization)
  G4double r = 0.001*std::sqrt(point[0]*point[0]+point[1]*point[1]); 
  G4double z = 0.001*point[2]; 
  
  //inside solenoid volume
  if(abs(z)>5.0 || r > 3.0) {
    bField[0] = bField[1] = bField[2] = 0.0;
  }
  else {

    //parameterized magnetic field in r,z (phi symmetric)
    G4double Bw[3] = {0,0,0};

    z -= This->prm[3];                    // max Bz point is shifted in z

    G4double az=abs(z);
    G4double zainv=z*This->ainv;
    G4double u=This->hlova-zainv;
    G4double v=This->hlova+zainv;
    G4double fu[4],gv[4];

    GPMagneticField_ffunkti(u,fu);
    GPMagneticField_ffunkti(v,gv);
    
    G4double rat=0.5*r*This->ainv;
    G4double rat2=rat*rat;

    Bw[0]=This->hb0*rat*(fu[1]-gv[1]-(fu[3]-gv[3])*rat2*0.5);
    Bw[2]=This->hb0*(fu[0]+gv[0]-(fu[2]+gv[2])*rat2);

    G4double corBr= This->prm[4]*r*z*(az-This->prm[5])*(az-This->prm[5]);
    G4double corBz=-This->prm[6]*(exp(-(z-This->prm[7])*(z-This->prm[7])*This->coeff) +
				  exp(-(z+This->prm[7])*(z+This->prm[7])*This->coeff) );
    Bw[0] += corBr;
    Bw[2] += corBz;

    //convert to xyz in tesla unit
    const G4double tesla = 0.001;

    G4double rinv=(r>0) ? 1/r : 0;
    bField[0]=tesla*Bw[0]*point[0]*rinv;
    bField[1]=tesla*Bw[0]*point[1]*rinv;
    bField[2]=tesla*Bw[2];

    //    printf("show Bxyz = %f %f %f \n",bField[0],bField[1],bField[2]);

  }
}

FQUALIFIER
void GPMagneticField_ffunkti(G4double u, G4double *ff)
{
  // Function and its 3 derivatives
  G4double a = 1.0/(1.0+u*u);
  G4double b = sqrt(a);

  ff[0]=u*b;
  ff[1]=a*b;
  ff[2]=-3.0*u*a*ff[1];
  ff[3]=a*ff[2]*((1.0/u)-4.0*u);
}

FQUALIFIER
void  GPMagneticField_SetGravityActive(GPMagneticField *This, 
				       G4bool OnOffFlag )
{ 
  This->fGravityActive = OnOffFlag; 
} 

FQUALIFIER
G4bool GPMagneticField_IsGravityActive(GPMagneticField *This) 
{ 
  return This->fGravityActive;
}

FQUALIFIER
G4bool GPMagneticField_DoesFieldChangeEnergy(GPMagneticField *This) 
{ 
  return false; 
}
