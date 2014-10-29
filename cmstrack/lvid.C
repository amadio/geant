void lvid(Int_t event=0, Int_t nmax=1000) {
   if (gFile == 0) new TFile ("event03.root");
   if (gFile == 0) return;
   TTree *T = (TTree*)gFile->Get(Form("Event%4.4d",event));
   TH1D *hlvid = new TH1D("hlvid","number of steps per lvid;lvid;number of steps",3500,0.01,3500.01);
   T->Draw("lvid>>hlvid","","nodraw");
   Int_t *index = new Int_t[3500];
   Double_t *array = (Double_t*)hlvid->GetArray();
   TMath::Sort(3500,array,index);
   TH1D *hlvidtop = new TH1D("hlvidtop","volume ids sorted by number of steps;lvidord;number of steps",nmax,1,nmax+1);
   TObjArray *volumes = (TObjArray*)gFile->Get("LogicalVolumes");
   for (Int_t i=0;i<nmax;i++) {
      Int_t ind = index[i];
      hlvidtop->SetBinContent(i+1,array[ind]);
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
   T->SetMarkerColor(kBlack);
   T->Draw("sign(y)*sqrt(sqrt(x*x+y*y)):sign(z)*sqrt(abs(z))","sqrt(x*x+y*y)<2500"); //fish eye view
   T->SetMarkerColor(kRed);
   T->Draw("sign(y)*sqrt(sqrt(x*x+y*y)):sign(z)*sqrt(abs(z))",Form("sqrt(x*x+y*y)<2500&&lvid==%-d",index[0]),"same");
   T->SetMarkerColor(kBlue);
   T->Draw("sign(y)*sqrt(sqrt(x*x+y*y)):sign(z)*sqrt(abs(z))",Form("sqrt(x*x+y*y)<2500&&lvid==%-d",index[1]),"same");
   T->SetMarkerColor(kGreen);
   T->Draw("sign(y)*sqrt(sqrt(x*x+y*y)):sign(z)*sqrt(abs(z))",Form("sqrt(x*x+y*y)<2500&&lvid==%-d",index[2]),"same");
   T->SetMarkerColor(kCyan);
   T->Draw("sign(y)*sqrt(sqrt(x*x+y*y)):sign(z)*sqrt(abs(z))",Form("sqrt(x*x+y*y)<2500&&lvid==%-d",index[3]),"same");

   pave1 = new TPave(0.75,0.76,0.98,0.94,0,"brNDC");
   pave1->SetShadowColor(0);
   pave1->SetLineColor(0);
   pave1->Draw();
   TText *t1 = new TText();
   t1->SetTextColor(kRed);
   t1->SetTextSize(0.03);
   t1->DrawTextNDC(0.77,0.90,Form("#1 volume %s",((TObjString*)volumes->At(index[0]))->GetString().Data()));
   TText *t2 = new TText();
   t2->SetTextColor(kBlue);
   t2->SetTextSize(0.03);
   t2->DrawTextNDC(0.77,0.86,Form("#2 volume %s",((TObjString*)volumes->At(index[1]))->GetString().Data()));
   TText *t3 = new TText();
   t3->SetTextColor(kGreen);
   t3->SetTextSize(0.03);
   t3->DrawTextNDC(0.77,0.82,Form("#3 volume %s",((TObjString*)volumes->At(index[2]))->GetString().Data()));
   TText *t4 = new TText();
   t4->SetTextColor(kCyan);
   t4->SetTextSize(0.03);
   t4->DrawTextNDC(0.77,0.78,Form("#4 volume %s",((TObjString*)volumes->At(index[3]))->GetString().Data()));
   //T->Draw("sign(y)*sqrt(x*x+y*y):z","sqrt(x*x+y*y)<2500");  //normal view
   c2->Print("c2.gif");
}
