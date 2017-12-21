void geomPoly()
{
   TGeoManager::Import("./polycone.root");
   gGeoManager->DefaultColors();

   //gGeoManager->GetVolume("ALIC")->Draw("ogl");
   new TBrowser;
}
