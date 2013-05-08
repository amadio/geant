#include "TGeoManager.h"
   
void a() 
{
   TGeoManager *testvec = new TGeoManager("testvec","This is a naive test");
   TGeoMaterial *vacmat=new TGeoMaterial("vacuum",0,0,0);  
   TGeoMedium *vacmed=new TGeoMedium("Vacuum",0,vacmat);

   TGeoMaterial *Fe=new TGeoMaterial("Fe",55.845,26,7.87); 
   TGeoMedium *Iron=new TGeoMedium("Iron",1,Fe);

   // create volume

   TGeoVolume *world=testvec->MakeBox("world",vacmed,1000,1000,1000);   
   testvec->SetTopVolume(world);     
   testvec->SetTopVisible(1); 


   TGeoVolume *Band=testvec->MakeBox("tbox",Iron,20,20,2.5);
     Band->SetLineColor(12);
     Band->SetFillColor(12);

   // drawing head
   world->AddNode(Band,1,new TGeoTranslation(0,0,90));


   //close geometry
   world->SetVisibility(0);
   testvec->CloseGeometry();

   // in GL viewer
   world->Draw("ogl");
}
