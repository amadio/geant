#include <TAxis.h>
#include <TCanvas.h>
#include <TEXsec.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLine.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPXsec.h>
#include <TPartIndex.h>
#include <TROOT.h>
#include <TString.h>
#include <TText.h>

ClassImp(TEXsec)

const char* TEXsec::fEleSymbol[NELEM]={"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg",
				"Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V",
				"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
				"Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh",
				"Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba",
				"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho",
				"Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt",
				"Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac",
				"Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
				"Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg",
				"Cn","Uut","Fl","Uup","Lv","Uus","Uuo"};
const char* TEXsec::fEleName[NELEM]={"Hydrogen","Helium","Lithium","Beryllium","Boron",
				"Carbon","Nitrogen","Oxygen","Fluorine","Neon",
				"Sodium","Magnesium","Aluminium","Silicon","Phosphorus",
				"Sulfur","Chlorine","Argon","Potassium","Calcium",
				"Scandium","Titanium","Vanadium","Chromium","Manganese",
				"Iron","Cobalt","Nickel","Copper","Zinc","Gallium",
				"Germanium","Arsenic","Selenium","Bromine","Krypton",
				"Rubidium","Strontium","Yttrium","Zirconium","Niobium",
				"Molybdenum","Technetium","Ruthenium","Rhodium","Palladium",
				"Silver","Cadmium","Indium","Tin","Antimony","Tellurium",
				"Iodine","Xenon","Caesium","Barium","Lanthanum","Cerium",
				"Praseodymium","Neodymium","Promethium","Samarium",
				"Europium","Gadolinium","Terbium","Dysprosium","Holmium",
				"Erbium","Thulium","Ytterbium","Lutetium","Hafnium",
				"Tantalum","Tungsten","Rhenium","Osmium","Iridium",
				"Platinum","Gold","Mercury","Thallium","Lead","Bismuth",
				"Polonium","Astatine","Radon","Francium","Radium","Actinium",
				"Thorium","Protactinium","Uranium","Neptunium","Plutonium",
				"Americium","Curium","Berkelium","Californium","Einsteinium",
				"Fermium","Mendelevium","Nobelium","Lawrencium","Rutherfordium",
				"Dubnium","Seaborgium","Bohrium","Hassium","Meitnerium",
				"Darmstadtium","Roentgenium","Copernicium","Ununtrium",
				"Flerovium","Ununpentium","Livermorium","Ununseptium",
				"Ununoctium"};

const Float_t TEXsec::fWElem[NELEM]={1.008,4.0026,6.94,9.0122,10.81,12.011,14.007,15.999,
					18.998,20.180,22.990,24.305,26.982,28.085,30.974,32.06,
					35.45,39.948,39.098,40.078,44.956,47.867,50.942,51.996,
					54.938,55.845,58.933,58.693,63.546,65.38,69.723,72.63,
					74.922,78.96,79.904,83.798,85.468,87.62,88.906,91.224,
					92.906,95.96,97.91,101.07,102.91,106.42,107.87,112.41,
					114.82,118.71,121.76,127.60,126.90,131.29,132.91,137.33,
					138.91,140.12,140.91,144.24,144.91,150.36,151.96,157.25,
					158.93,162.50,164.93,167.26,168.93,173.05,174.97,178.49,
					180.95,183.84,186.21,190.23,192.22,195.08,196.97,200.59,
					204.38,207.2,208.98,208.98,209.99,222.02,223.02,226.03,
					227.03,232.04,231.04,238.03,237.05,244.06,243.06,247.07,
					247.07,251.08,252.08,257.10,258.10,259.10,262.11,265.12,
					268.13,271.13,270,277.15,276.15,281.16,280.16,285.17,
					284.18,289.19,288.19,293,294,294};

TEXsec* TEXsec::fElements[NELEM]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
				  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
				  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

Int_t TEXsec::fNLdElems=0; 


//___________________________________________________________________
TEXsec::TEXsec():
   fEle(0),
   fAtcm3(0),
   fEmin(0),
   fEmax(0),
   fNEbins(0),
   fElDelta(0),
   fNRpart(0),
   fPXsec(0),
   fCuts(0)
{
}

//___________________________________________________________________
TEXsec::TEXsec(Int_t z, Int_t a, Float_t emin, Float_t emax, Int_t nen, Int_t np):
   TNamed(fEleSymbol[z-1],fEleName[z-1]),
   fEle(z*10000+a*10),
   fAtcm3(TMath::Na()*1e-24/fWElem[z-1]),
   fEmin(emin),
   fEmax(emax),
   fNEbins(nen),
   fElDelta(TMath::Log(fEmax/fEmin)/fNEbins),
   fNRpart(np),
   fPXsec(new TPXsec[fNRpart]),
   fCuts(new Double_t[fNRpart])
{
   memset(fCuts,0,fNRpart*sizeof(Double_t));
}

//___________________________________________________________________
TEXsec::~TEXsec() {
   delete [] fPXsec;
   delete [] fCuts;
}

//___________________________________________________________________
Bool_t TEXsec::AddPart(Int_t kpart, Int_t pdg, Int_t nen, Int_t nxsec, Float_t emin, Float_t emax) {
   return fPXsec[kpart].SetPart(pdg,nen,nxsec,emin,emax);
}

//___________________________________________________________________
Bool_t TEXsec::AddPart(Int_t kpart, Int_t pdg, Int_t nxsec) {
   return fPXsec[kpart].SetPart(pdg,nxsec);
}

//___________________________________________________________________
Bool_t TEXsec::AddPartMS(Int_t kpart, const Float_t angle[], const Float_t ansig[],
			 const Float_t length[], const Float_t lensig[]) {
   return fPXsec[kpart].SetPartMS(angle, ansig, length, lensig);
}

//___________________________________________________________________
Bool_t TEXsec::AddPartXS(Int_t kpart, const Float_t xsec[], const Short_t dict[]) {
   return fPXsec[kpart].SetPartXS(xsec,dict);
}

//___________________________________________________________________
Bool_t TEXsec::AddPartIon(Int_t kpart, const Float_t dedx[]) {
   return fPXsec[kpart].SetPartIon(dedx);
}

//___________________________________________________________________
Float_t TEXsec::XSPDG(Int_t pdg, Short_t rcode, Float_t en) const {
   for(Int_t i=0; i<fNRpart; ++i) 
      if(pdg == fPXsec[i].PDG()) 
	 return fPXsec[i].XS(TPartIndex::I()->ProcIndex(rcode),en);
   return -1;
}

//___________________________________________________________________
Float_t TEXsec::XS(Int_t pindex, Short_t rindex, Float_t en) const {
   return fPXsec[pindex].XS(rindex,en);
}

//___________________________________________________________________
Float_t TEXsec::DEdxPDG(Int_t pdg, Float_t en) const {
   for(Int_t i=0; i<fNRpart; ++i) 
      if(pdg == fPXsec[i].PDG()) 
	 return fPXsec[i].DEdx(en);
   return -1;
}

//___________________________________________________________________
Float_t TEXsec::DEdx(Int_t pindex, Float_t en) const {
   return fPXsec[pindex].DEdx(en);
}

//___________________________________________________________________
Bool_t TEXsec::MSPDG(Int_t pdg, Float_t en, Float_t &ang, Float_t &asig, 
		  Float_t &len, Float_t &lsig) const {
   for(Int_t i=0; i<fNRpart; ++i) 
      if(pdg == fPXsec[i].PDG()) 
	 return fPXsec[i].MS(en,ang,asig,len,lsig);
   return kFALSE;
}

//___________________________________________________________________
Bool_t TEXsec::MS(Int_t pindex, Float_t en, Float_t &ang, Float_t &asig, 
		  Float_t &len, Float_t &lsig) const {
   return fPXsec[pindex].MS(en,ang,asig,len,lsig);
}

//___________________________________________________________________
void TEXsec::DumpPointers() const {
   printf("Material %d emin %f emax %f NEbins %d ElDelta %f Npart %d\n",
	  fEle,fEmin,fEmax,fNEbins,fElDelta,fNRpart);
   for(Int_t i=0; i<fNRpart; ++i) 
      fPXsec[i].Dump();
}

//___________________________________________________________________
TGraph* TEXsec::XSGraph(const char* part, const char *reac, 
			Float_t emin, Float_t emax, Int_t nbin) const 
{
   Char_t title[200];
   const Double_t delta = TMath::Exp(TMath::Log(emax/emin)/(nbin-1));
   Float_t *xsec = new Float_t[nbin];
   Float_t *energy = new Float_t[nbin];
   Double_t en=emin;
   Int_t pindex = TPartIndex::I()->PartIndex(part);
   if(pindex < 0) {
      Error("XSGraph","Unknown particle %s\n",part);
      return 0;
   }
   Int_t proc = TPartIndex::I()->ProcIndex(reac);
   for(Int_t i=0; i<nbin; ++i) {
      energy[i] = en;
      xsec[i] = XS(pindex,proc,en);
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
TGraph* TEXsec::DEdxGraph(const char* part, 
			  Float_t emin, Float_t emax, Int_t nbin) const 
{
   Char_t title[200];
   const Double_t delta = TMath::Exp(TMath::Log(emax/emin)/(nbin-1));
   Float_t *dedx = new Float_t[nbin];
   Float_t *energy = new Float_t[nbin];
   Double_t en=emin;
   Int_t pindex = TPartIndex::I()->PartIndex(part);
   if(pindex < 0) {
      Error("DEdxGraph","Unknown particle %s\n",part);
      return 0;
   }
   for(Int_t i=0; i<nbin; ++i) {
      energy[i] = en;
      dedx[i] = DEdxPDG(pindex,en);
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
Float_t TEXsec::LambdaPDG(Int_t pdg, Double_t en) const {
   Double_t xs=0;
   for(Int_t i=0; i<fNRpart; ++i) 
      if(pdg == fPXsec[i].PDG()) {
	 xs = fPXsec[i].XS(TPartIndex::I()->NProc()-1,en);
	 break;
      }
   return xs?1./(fAtcm3*xs):TMath::Limits<Float_t>::Max();
}

//___________________________________________________________________
Float_t TEXsec::Lambda(Int_t pindex, Double_t en) const {
   Double_t xs=0;
   xs = fPXsec[pindex].XS(TPartIndex::I()->NProc()-1,en);
   return xs?1./(fAtcm3*xs):TMath::Limits<Float_t>::Max();
}

//___________________________________________________________________
Int_t TEXsec::SampleReac(Int_t pindex, Double_t en) const {
   return fPXsec[pindex].SampleReac(en);
}

//___________________________________________________________________
Int_t TEXsec::SampleReacPDG(Int_t pdg, Double_t en) const {
   for(Int_t i=0; i<fNRpart; ++i) 
      if(pdg == fPXsec[i].PDG()) {
	 return fPXsec[i].SampleReac(en);
      }
   return -1;
}

//___________________________________________________________________
TGraph* TEXsec::MSGraph(const char* part, const char* what,
			  Float_t emin, Float_t emax, Int_t nbin) const 
{
   Char_t title[200];
   const Char_t *whatname[4] = {"MSangle", "MSangle_sig", "MSCorr", "MSCorr_sig"};
   const Double_t delta = TMath::Exp(TMath::Log(emax/emin)/(nbin-1));
   Float_t *mscat = new Float_t[nbin];
   Float_t *energy = new Float_t[nbin];
   Double_t en=emin;
   Int_t pindex = TPartIndex::I()->PartIndex(part);
   if(pindex < 0) {
      Error("MSGraph","Unknown particle %s\n",part);
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
      MS(pindex,en,answ[0],answ[1],answ[2],answ[3]);
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
TEXsec *TEXsec::GetElement(Int_t z, Int_t a, TFile* f) {
   Int_t ecode = z*10000+a*10;
   for(Int_t el=0; el<fNLdElems; ++el) 
      if(ecode == fElements[el]->Ele()) return fElements[el];
   
   TFile *ff=gFile;
   if(f) ff=f;
   fElements[fNLdElems] = (TEXsec *) ff->Get(fEleSymbol[z-1]);
   return fElements[fNLdElems++];
}

//___________________________________________________________________
void TEXsec::Draw(Option_t *option)
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
   Int_t pindex = TPartIndex::I()->PartIndex(part);
   printf("pdg %d pindex %d\n",pdg,pindex);
   if(pindex < 0) {
      Error("Draw","Unknown particle %s\n",part);
      return;
   }
   if(reactions.Contains("All") || reactions.Contains("*")) {
      TString allrea="";
      for(Int_t i=0; i<TPartIndex::I()->NProc()-1; ++i) {
	 if(XS(pindex,i,emin)>=0) allrea=allrea+TPartIndex::I()->ProcName(i)+"|";
      }
      allrea+=TPartIndex::I()->ProcName(TPartIndex::I()->NProc()-1);
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
   //   if(!tc) tc = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
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
	 if(TPartIndex::I()->ProcIndex(reac)<0) {
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
   if(!same) {
      ((TAxis*) tmg->GetHistogram()->GetXaxis())->SetTitle("Energy (GeV)");
      ((TAxis*) tmg->GetHistogram()->GetYaxis())->SetTitle(ytitle);
   }
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

