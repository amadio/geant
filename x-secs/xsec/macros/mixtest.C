void mixtest() {
   gSystem->Load("libXsec");
   TFile *f = new TFile("xsec.root");
   Int_t z[2]={1,16};
   Float_t w[2]={2,1};
   Int_t a[2];
   tm=new TMXsec(z,a,w,2,1,kFALSE);
   w[0] = 2*TEXsec::WEle(1)/(2*TEXsec::WEle(1)+TEXsec::WEle(16));
   w[1] = TEXsec::WEle(16)/(2*TEXsec::WEle(1)+TEXsec::WEle(16));
   tm=new TMXsec(z,a,w,2,1,kTRUE);
}
