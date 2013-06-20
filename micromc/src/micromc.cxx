#include <TGeoManager.h>

Int_t main () {
   Char_t *geofile = "geom.root";
   TGeoManager *geom = TGeoManager::Import(geofile);
   
   // loop materials

   TObjArray *matlist = geom->GetListOfMaterials();
   TIter next(matlist);
   TGeoMaterial *mat=0;
   TGeoMixture *mix=0;
   while((mat = (TGeoMaterial*) next())) {
      if(!mat->IsUsed()) continue;
      Int_t nelem = mat->GetNelements();
      Int_t *z = new Int_t[nelem];
      Int_t *a = new Int_t[nelem];
      Float_t *w = new Float_t[nelem];
      for(Int_t iel=0; iel<nelem; ++iel) {
	 mat->GetElementProp(a[iel],z[iel],w[iel],iel);
	 printf("Mixture %s element %s z %d a %d\n",
		mat->GetName(), mat->GetElement(iel)->GetName(),
		z[iel],a[iel]);
      }
   }
   return 0;
}
