void graphxsec(const char* mat, const char* pnam, const char* reac,
	       Float_t emin=1e-3,Float_t emax=1e6,Int_t nbin=100)
{
   gSystem->Load("libXSec");
   TFile *fx = new TFile("xsec.root");
   TPartIndex* tp = fx->Get("PartIndex");
   TMXsec *mate = fx->Get(mat);
   //TGraph *gr = mate->XSGraph(pnam,reac,emin,emax,nbin); 
   //   gr->Draw("ACP");
   mate->XSDraw(pnam,reac,emin,emax,nbin);
}
