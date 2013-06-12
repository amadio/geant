#include <TMXsec.h>
#include <TPXsec.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TPartIndex.h>
#include <TString.h>
#include <TObjString.h>
#include <TFile.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TText.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TH1F.h>

ClassImp(TMXsec)

//___________________________________________________________________
TMXsec::TMXsec():
   fMat(0),
   fAtcm3(0),
   fEmin(0),
   fEmax(0),
   fNEbins(0),
   fElDelta(0),
   fNpart(0),
   fPXsec(0),
   fCuts(0)
{
}

//___________________________________________________________________
TMXsec::TMXsec(Int_t z, Int_t a, Float_t emin, Float_t emax, Int_t nen, Int_t np):
   TNamed(TPartIndex::MatSymb(z),TPartIndex::MatName(z)),
   fMat(z*10000+a*10),
   fAtcm3(TMath::Na()*1e-24/TPartIndex::I()->WMat(z)),
   fEmin(emin),
   fEmax(emax),
   fNEbins(nen),
   fElDelta(TMath::Log(fEmax/fEmin)/fNEbins),
   fNpart(np),
   fPXsec(new TPXsec[fNpart]),
   fCuts(new Double_t[fNpart])
{
}

//___________________________________________________________________
TMXsec::~TMXsec() {
   delete [] fPXsec;
   delete [] fCuts;
}

//___________________________________________________________________
Bool_t TMXsec::AddPart(Int_t kpart, Int_t pdg, Int_t nen, Int_t nxsec, Float_t emin, Float_t emax) {
   return fPXsec[kpart].SetPart(pdg,nen,nxsec,emin,emax);
}

//___________________________________________________________________
Bool_t TMXsec::AddPartMS(Int_t kpart, const Float_t angle[], const Float_t ansig[],
			 const Float_t length[], const Float_t lensig[]) {
   return fPXsec[kpart].SetPartMS(angle, ansig, length, lensig);
}

//___________________________________________________________________
Bool_t TMXsec::AddPartXS(Int_t kpart, const Float_t xsec[], const Short_t dict[]) {
   return fPXsec[kpart].SetPartXS(xsec,dict);
}

//___________________________________________________________________
Bool_t TMXsec::AddPartIon(Int_t kpart, const Float_t dedx[]) {
   return fPXsec[kpart].SetPartIon(dedx);
}

//___________________________________________________________________
Bool_t TMXsec::Finalise() {
   Bool_t retok = kTRUE;
   for(Int_t i=0; i<fNpart; ++i) 
      retok &= fPXsec[i].Finalise();
   return retok;
}

//___________________________________________________________________
Float_t TMXsec::XSPDG(Int_t pdg, Short_t rcode, Float_t en) const {
   for(Int_t i=0; i<fNpart; ++i) 
      if(pdg == fPXsec[i].PDG()) 
	 return fPXsec[i].XS(TPartIndex::I()->ProcIndex(rcode),en);
   return -1;
}

//___________________________________________________________________
Float_t TMXsec::XS(Int_t pindex, Short_t rindex, Float_t en) const {
   return fPXsec[pindex].XS(rindex,en);
}

//___________________________________________________________________
Float_t TMXsec::DEdxPDG(Int_t pdg, Float_t en) const {
   for(Int_t i=0; i<fNpart; ++i) 
      if(pdg == fPXsec[i].PDG()) 
	 return fPXsec[i].DEdx(en);
   return -1;
}

//___________________________________________________________________
Float_t TMXsec::DEdx(Int_t pindex, Float_t en) const {
   return fPXsec[pindex].DEdx(en);
}

//___________________________________________________________________
Bool_t TMXsec::MSPDG(Int_t pdg, Float_t en, Float_t &ang, Float_t &asig, 
		  Float_t &len, Float_t &lsig) const {
   for(Int_t i=0; i<fNpart; ++i) 
      if(pdg == fPXsec[i].PDG()) 
	 return fPXsec[i].MS(en,ang,asig,len,lsig);
   return kFALSE;
}

//___________________________________________________________________
Bool_t TMXsec::MS(Int_t pindex, Float_t en, Float_t &ang, Float_t &asig, 
		  Float_t &len, Float_t &lsig) const {
   return fPXsec[pindex].MS(en,ang,asig,len,lsig);
}

//___________________________________________________________________
void TMXsec::DumpPointers() const {
   printf("Material %d emin %f emax %f NEbins %d ElDelta %f Npart %d\n",
	  fMat,fEmin,fEmax,fNEbins,fElDelta,fNpart);
   for(Int_t i=0; i<fNpart; ++i) 
      fPXsec[i].Dump();
}

//___________________________________________________________________
TGraph* TMXsec::XSGraph(const char* part, const char *reac, 
			Float_t emin, Float_t emax, Int_t nbin) const 
{
   Char_t title[200];
   const Double_t delta = TMath::Exp(TMath::Log(emax/emin)/(nbin-1));
   Float_t *xsec = new Float_t[nbin];
   Float_t *energy = new Float_t[nbin];
   Double_t en=emin;
   Int_t rcode = TPartIndex::I()->Rcode(reac);
   if(rcode<0) {
      Error("XSGraph","Unknown reaction %s\n",reac);
      return 0;
   }
   Int_t pdg = TPartIndex::I()->PDG(part);
   if(pdg == -12345678) {
      Error("XSGraph","Unknown particle %s\n",part);
      return 0;
   }
   for(Int_t i=0; i<nbin; ++i) {
      energy[i] = en;
      xsec[i] = XSPDG(pdg,rcode,en);
      en*=delta;
   }
   TGraph *tg = new TGraph(nbin,energy,xsec);
   memset(title,0,200);
   snprintf(title,199,"%s %s on %s",part,reac,GetTitle());
   tg->SetTitle(title);
   delete [] xsec;
   delete [] energy;
   return tg;
}

//___________________________________________________________________
TGraph* TMXsec::DEdxGraph(const char* part, 
			  Float_t emin, Float_t emax, Int_t nbin) const 
{
   Char_t title[200];
   const Double_t delta = TMath::Exp(TMath::Log(emax/emin)/(nbin-1));
   Float_t *dedx = new Float_t[nbin];
   Float_t *energy = new Float_t[nbin];
   Double_t en=emin;
   Int_t pdg = TPartIndex::I()->PDG(part);
   if(pdg == -12345678) {
      Error("XSGraph","Unknown particle %s\n",part);
      return 0;
   }
   for(Int_t i=0; i<nbin; ++i) {
      energy[i] = en;
      dedx[i] = DEdxPDG(pdg,en);
      en*=delta;
   }
   TGraph *tg = new TGraph(nbin,energy,dedx);
   memset(title,0,200);
   snprintf(title,199,"%s dEdx on %s",part,GetTitle());
   tg->SetTitle(title);
   delete [] dedx;
   delete [] energy;
   return tg;
}

//___________________________________________________________________
Float_t TMXsec::LambdaPDG(Int_t pdg, Double_t en) const {
   for(Int_t i=0; i<fNpart; ++i) 
      if(pdg == fPXsec[i].PDG()) 
	 return fPXsec[i].XS(TPartIndex::I()->NReac()-1,en);
   return TMath::Limits<Float_t>::Max();
}

//___________________________________________________________________
Float_t TMXsec::Lambda(Int_t pindex, Double_t en) const {
   return fPXsec[pindex].XS(TPartIndex::I()->NReac()-1,en);
}

//___________________________________________________________________
TGraph* TMXsec::MSGraph(const char* part, const char* what,
			  Float_t emin, Float_t emax, Int_t nbin) const 
{
   Char_t title[200];
   const Char_t *whatname[4] = {"MSangle", "MSangle_sig", "MSCorr", "MSCorr_sig"};
   const Double_t delta = TMath::Exp(TMath::Log(emax/emin)/(nbin-1));
   Float_t *mscat = new Float_t[nbin];
   Float_t *energy = new Float_t[nbin];
   Double_t en=emin;
   Int_t pdg = TPartIndex::I()->PDG(part);
   if(pdg == -12345678) {
      Error("XSGraph","Unknown particle %s\n",part);
      return 0;
   }
   Int_t iopt=4;
   while(iopt--) if(!strcmp(whatname[iopt],what)) break;
   if(iopt<0)  {
      Error("MSGraph","Unknown parameter %s\nShould be one of %s %s %s %s\n",
	    what, whatname[0],whatname[1],whatname[2],whatname[3]);
      return 0;
   }
   for(Int_t i=0; i<nbin; ++i) {
      Float_t answ[4];
      energy[i] = en;
      MSPDG(pdg,en,answ[0],answ[1],answ[2],answ[3]);
      mscat[i] = answ[iopt];
      en*=delta;
   }
   TGraph *tg = new TGraph(nbin,energy,mscat);
   memset(title,0,200);
   snprintf(title,199,"%s %s on %s",part,whatname[iopt],GetTitle());
   tg->SetTitle(title);
   delete [] mscat;
   delete [] energy;
   return tg;
}

//___________________________________________________________________
void TMXsec::Draw(Option_t *option)
{
   // Draw cross sections and other physics quantities for this material
   //
   TString help("Method to draw the cross sections. The format is:\n");
   help+="particle,reaction|reaction|reaction,emin,emax,nbin:<options for draw>\n";
   help+="Available reactions are:\n";
   help+="Total,Transport,MultScatt,Ionisation,Decay,inElastic,Elastic,Capture,Brehms\n";
   help+="PairProd,dEdx,MSang,MSang_sig,MSCorr,MSCorr_sig\n";
   help+="the option All can be given to draw all reaction cross sections for a given particle\n";
   help+="the option \"same\" superimposes the plot on the previous one";

   // CoulombScatt,Photoel,Compton,Conversion,Capture,Killer,dEdx,MSangle,MSlength
   const EColor col[14] = {kBlack,kGray,kRed,kGreen,kBlue, kMagenta, kCyan,
		       kOrange, kSpring, kTeal, kAzure, kViolet, kPink };
   Char_t title[200];
   Char_t gtitle[200];

   TString opt = option;
   
   static Int_t isame=0;

   TObjArray *sections = opt.Tokenize(":");
   TString sec1 = ((TObjString*) sections->At(0))->GetName();
   TString sec2;
   if(sections->GetEntries()>1) sec2 = ((TObjString*) sections->At(1))->GetName();
   Bool_t same = sec2.Contains("same");

   TObjArray *token = sec1.Tokenize(",");
   Int_t narg = token->GetEntries();
   if(narg<2) {
      Info("Draw","%s",help.Data());
      return;
   }
   const Char_t *part = ((TObjString *) token->At(0))->GetName();
   TString reactions = ((TObjString*) token->At(1))->GetName();

   Float_t emin=1e-3;
   if(narg>2) sscanf(((TObjString*) token->At(2))->GetName(),"%f",&emin);
   Float_t emax=1e6;
   if(narg>3) sscanf(((TObjString*) token->At(3))->GetName(),"%f",&emax);
   Int_t nbin=100;
   if(narg>4) sscanf(((TObjString*) token->At(4))->GetName(),"%d",&nbin);
   if(gFile) gFile->Get("PartIndex");
   Int_t pdg = TPartIndex::I()->PDG(part);
   if(pdg == -12345678) {
      Error("Draw","Particle %s does not exist\n",part);
      TPartIndex::I()->Print("particles");
      return;
   }
   if(reactions.Contains("All") || reactions.Contains("*")) {
      TString allrea="";
      for(Int_t i=0; i<TPartIndex::I()->NReac()-1; ++i) {
	 if(XSPDG(pdg,TPartIndex::I()->ProcCode(i),emin)>=0) allrea=allrea+TPartIndex::I()->ReacName(i)+"|";
      }
      allrea+=TPartIndex::I()->ReacName(TPartIndex::I()->NReac()-1);
      reactions.ReplaceAll("All",allrea);
      reactions.ReplaceAll("*",allrea);
   }
   TObjArray *rnames = reactions.Tokenize("|");
   Int_t nreac = rnames->GetEntries();
   snprintf(gtitle,199,"%s %s on %s",part,reactions.ReplaceAll("|",",").Data(),GetTitle());
   TMultiGraph *tmg= new TMultiGraph("G5",gtitle);
   TLine **line = new TLine*[nreac];
   TText **text = new TText*[nreac];
   Float_t lstartx = 0.7;

   TCanvas *tc=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("G5canvas");
   if(!tc) {
      tc = new TCanvas("G5canvas",gtitle,600,400);      
      tc->SetLogx();
      tc->SetLogy();
      tc->SetGrid();
      same = kFALSE;
      sec2.ReplaceAll("same","");
   }

   if (same) {
      sec2.ReplaceAll("same","C");
      ++isame;
   } else {
      isame=0;
      tc->Clear();
      tc->SetTitle(gtitle);
   }

   if(isame == 2 || isame == 3) lstartx = 0.3;
   const Float_t lendx = lstartx+0.05;
   const Float_t tstartx = lstartx+0.07;
   Float_t lstarty = 0.85;
   if(isame==1 || isame==3) lstarty = 0.4;
   const Float_t lstepy = 0.03;
   Int_t coff = isame;
   Char_t ytitle[50];
   for(Int_t j=0; j<nreac; ++j) {
      const Char_t *reac = ((TObjString*) rnames->At(j))->GetName();
      TGraph *tg;
      if(TString(reac).BeginsWith("MS")) {
	 tg = MSGraph(part,reac,emin,emax,nbin);
	 snprintf(title,199,"%s %s on %s",part,reac,GetTitle());
	 tg->SetName(title);
	 tg->SetTitle(title);
	 tg->SetLineColor(col[(coff+j)%14]);
	 tmg->Add(tg,sec2.Data());
	 line[j] = new TLine(lstartx,lstarty-lstepy*j,lendx,lstarty-lstepy*j);
	 line[j]->SetLineColor(col[(coff+j)%14]);
	 text[j]=new TText(tstartx,lstarty-lstepy*j,reac);
	 if(TString(reac).BeginsWith("MSangle")) snprintf(ytitle,49,"Radians");
	 else snprintf(ytitle,49,"Relative Step Correction");
      } else if(!strcmp(reac,"dEdx")) {
	 snprintf(ytitle,49,"MeV/cm");
	 tg = DEdxGraph(part, emin, emax, nbin);
	 snprintf(title,199,"%s dEdx on %s",part,GetTitle());
	 tg->SetName(title);
	 tg->SetTitle(title);
	 tg->SetLineColor(col[(coff+j)%14]);
	 tmg->Add(tg,sec2.Data());
	 line[j] = new TLine(lstartx,lstarty-lstepy*j,lendx,lstarty-lstepy*j);
	 line[j]->SetLineColor(col[(coff+j)%14]);
	 text[j]=new TText(tstartx,lstarty-lstepy*j,reac);
      } else {
	 snprintf(ytitle,49,"barn");
	 if(TPartIndex::I()->Rcode(reac)<0) {
	    Error("Draw","Reaction %s does not exist\n",reac);
	    TPartIndex::I()->Print("reactions");
	    printf("dEdx, MSangle, MSangle_sig, MSCorr, MSCorr_sig\n");
	    return;
	 }
	 tg = XSGraph(part, reac, emin, emax, nbin);
	 snprintf(title,199,"%s %s on %s",part,reac,GetTitle());
	 tg->SetName(title);
	 tg->SetTitle(title);
	 tg->SetLineColor(col[(coff+j)%14]);
	 tmg->Add(tg,sec2.Data());
	 line[j] = new TLine(lstartx,lstarty-lstepy*j,lendx,lstarty-lstepy*j);
	 line[j]->SetLineColor(col[(coff+j)%14]);
	 text[j]=new TText(tstartx,lstarty-lstepy*j,reac);
      }
   }
   const Char_t *gopt=0;
   if(strlen(sec2.Data())) gopt = sec2.Data();
   else gopt = "AC";
   tmg->SetMinimum(1e-6);
   tmg->Draw(gopt);
   ((TAxis*) tmg->GetHistogram()->GetXaxis())->SetTitle("Energy (GeV)");
   ((TAxis*) tmg->GetHistogram()->GetYaxis())->SetTitle(ytitle);
   TText **ptext = new TText*[nreac];
   Char_t string[100]={"\0"};
   if(same) {
      snprintf(string,99,"%s on %s",part,GetTitle());
   }
   for(Int_t j=0;j<nreac;++j) {
      line[j]->SetBit(TLine::kLineNDC);
      line[j]->Draw();
      line[j]->SetLineColor(col[(coff+j)%14]);
      text[j]->SetNDC(kTRUE);
      text[j]->SetTextSize(0.03);
      text[j]->SetTextAlign(12);
      text[j]->Draw();
      if(same) {
	 ptext[j]=new TText(lstartx*0.95,lstarty-lstepy*j,string);
	 ptext[j]->SetTextSize(0.03);
	 ptext[j]->SetNDC(kTRUE);
	 ptext[j]->SetTextAlign(32);
	 ptext[j]->Draw();
      }
   }
}

