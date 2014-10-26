void lvid(Int_t nmax=1000) {
   if (gFile == 0) new TFile ("event03.root");
   if (gFile == 0) return;
   TTree *T = (TTree*)gFile->Get("Event0000");
   TH1D *hlvid = new TH1D("hlvid","number of steps per lvid;lvid;number of steps",3500,0.01,3500.01);
   T->Draw("lvid>>hlvid","","nodraw");
   Int_t *index = new Int_t[3500];
   Double_t *array = (Double_t*)hlvid->GetArray();
   TMath::Sort(3500,array,&index[1]);
   TH1D *hlvidtop = new TH1D("hlvidtop","volume ids sorted by number of steps;lvidord;number of steps",nmax,1,nmax+1);
   TObjArray *volumes = (TObjArray*)gFile->Get("LogicalVolumes");
   for (Int_t i=0;i<nmax;i++) {
      Int_t ind = index[i];
      hlvidtop->SetBinContent(i,array[ind]);
      hlvidtop->GetXaxis()->SetBinLabel(i+1,volumes->At(ind)->GetName());
   }
   hlvidtop->GetXaxis()->SetRange(1,20);
   TCanvas *c1 = new TCanvas("c1","c1",1000,800);
   c1->Divide(1,2);
   TPad *pad1 = (TPad*)c1->cd(1);
   pad1->SetLogy();
   hlvid->Draw();
   TPad *pad2 =(TPad*)c1->cd(2);
   pad2->SetLogy();
   hlvidtop->SetStats(kFALSE);
   hlvidtop->Draw();
   Double_t percent = hlvidtop->Integral()/hlvid->Integral();
   printf("percent = %g\n",percent);
   TPaveText *pave = new TPaveText(0.5,0.4,0.85,0.85,"brNDC");
   Double_t nsteps = hlvid->GetEntries();
   pave->AddText(Form("per cent of steps in    1 lvids = %g",100*hlvidtop->Integral(1,   1)/nsteps));
   pave->AddText(Form("per cent of steps in    5 lvids = %g",100*hlvidtop->Integral(1,   5)/nsteps));
   pave->AddText(Form("per cent of steps in   10 lvids = %g",100*hlvidtop->Integral(1,  10)/nsteps));
   pave->AddText(Form("per cent of steps in   50 lvids = %g",100*hlvidtop->Integral(1,  50)/nsteps));
   pave->AddText(Form("per cent of steps in  100 lvids = %g",100*hlvidtop->Integral(1, 100)/nsteps));
   pave->AddText(Form("per cent of steps in  200 lvids = %g",100*hlvidtop->Integral(1, 200)/nsteps));
   pave->AddText(Form("per cent of steps in  500 lvids = %g",100*hlvidtop->Integral(1, 500)/nsteps));
   pave->AddText(Form("per cent of steps in 1000 lvids = %g",100*hlvidtop->Integral(1,1000)/nsteps));
   pave->Draw();
   c1->Print("c1.pdf");
   TCanvas *c2 = new TCanvas("c2","c2");
   T->Draw("sign(y)*sqrt(sqrt(x*x+y*y)):sign(z)*sqrt(abs(z))","sqrt(x*x+y*y)<2500"); //fish eye view
   //T->Draw("sign(y)*sqrt(x*x+y*y):z","sqrt(x*x+y*y)<2500");  //normal view
   c2->Print("c2.gif");
}
