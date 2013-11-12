#include "TAxis.h"
#include "TCanvas.h"
#include "TEXsec.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include <TObjArray.h>
#include <TObjString.h>
#include <TPXsec.h>
#include <TPartIndex.h>
#include <TROOT.h>
#include <TString.h>
#include <TText.h>
#include <TFrame.h>
#include <TGFrame.h>
#include <TGComboBox.h>
#include <TGListBox.h>
#include <TGLabel.h>
#include <TRootEmbeddedCanvas.h>

ClassImp(TEXsec)


TEXsec* TEXsec::fElements[NELEM]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

Int_t TEXsec::fNLdElems=0;

TGMainFrame         *TEXsec::fMain=0;
TGHorizontalFrame   *TEXsec::fSecond=0;
TRootEmbeddedCanvas *TEXsec::fCanvas=0;
TGListBox           *TEXsec::fReactionBox=0;
TGListBox           *TEXsec::fParticleBox=0;

//___________________________________________________________________
TEXsec::TEXsec():
fEle(0),
fDens(0),
fAtcm3(0),
fEmin(0),
fEmax(0),
fNEbins(0),
fEilDelta(0),
fEGrid(TPartIndex::I()->EGrid()),
fNRpart(0),
fPXsec(0)
{
}

//___________________________________________________________________
TEXsec::TEXsec(Int_t z, Int_t a, Float_t dens, Int_t np):
TNamed(TPartIndex::I()->EleSymb(z),TPartIndex::I()->EleName(z)),
fEle(z*10000+a*10),
fDens(dens),
fAtcm3(fDens*TMath::Na()*1e-24/TPartIndex::I()->WEle(z)),
fEmin(TPartIndex::I()->Emin()),
fEmax(TPartIndex::I()->Emax()),
fNEbins(TPartIndex::I()->NEbins()),
fEilDelta(TPartIndex::I()->EilDelta()),
fEGrid(TPartIndex::I()->EGrid()),
fNRpart(np),
fPXsec(new TPXsec[fNRpart])
{
    memset(fCuts,0,4*sizeof(Float_t));
}

//___________________________________________________________________
TEXsec::~TEXsec() {
    delete [] fPXsec;
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
Bool_t TEXsec::AddPartXS(Int_t kpart, const Float_t xsec[], const Int_t dict[]) {
    return fPXsec[kpart].SetPartXS(xsec,dict);
}

//___________________________________________________________________
Bool_t TEXsec::AddPartIon(Int_t kpart, const Float_t dedx[]) {
    return fPXsec[kpart].SetPartIon(dedx);
}

//___________________________________________________________________
Float_t TEXsec::XS(Int_t pindex, Int_t rindex, Float_t en) const {
    return fPXsec[pindex].XS(rindex,en);
}

//___________________________________________________________________
Float_t TEXsec::DEdx(Int_t pindex, Float_t en) const {
    return fPXsec[pindex].DEdx(en);
}

//___________________________________________________________________
Bool_t TEXsec::MS(Int_t pindex, Float_t en, Float_t &ang, Float_t &asig,
                  Float_t &len, Float_t &lsig) const {
    return fPXsec[pindex].MS(en,ang,asig,len,lsig);
}

//___________________________________________________________________
void TEXsec::DumpPointers() const {
    printf("Material %d emin %f emax %f NEbins %d EilDelta %f Npart %d\n",
           fEle,fEmin,fEmax,fNEbins,fEilDelta,fNRpart);
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
        dedx[i] = DEdx(pindex,en);
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
Float_t TEXsec::Lambda(Int_t pindex, Double_t en) const {
  Double_t xs = fPXsec[pindex].XS(TPartIndex::I()->NProc()-1,en);
  return xs?1./(fAtcm3*xs):TMath::Limits<Float_t>::Max();
}

//___________________________________________________________________
Bool_t TEXsec::Lambda_v(Int_t npart, const Int_t pindex[], const Double_t en[], Double_t lam[]) const {
  const Int_t itot = TPartIndex::I()->NProc()-1;
  for(Int_t ip=0; ip<npart; ++ip) {
    Double_t xs = fPXsec[pindex[ip]].XS(itot,en[ip]);
    lam[ip] = xs?1./(fAtcm3*xs):TMath::Limits<Float_t>::Max();
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t TEXsec::Lambda_v(Int_t npart, Int_t pindex, const Double_t en[], Double_t lam[]) const {
  const Int_t itot = TPartIndex::I()->NProc()-1;
  return fPXsec[pindex].XS_v(npart, itot,en,lam);
}

//___________________________________________________________________
Int_t TEXsec::SampleReac(Int_t pindex, Double_t en) const {
    return fPXsec[pindex].SampleReac(en);
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
    //   printf("Getting Element %d %d %d\n",z,a,fNLdElems);
    Int_t ecode = z*10000+a*10;
    for(Int_t el=0; el<fNLdElems; ++el)
        if(ecode == fElements[el]->Ele())
            return fElements[el];
    
    // Element not found in memory, getting it from file
    TFile *ff=gFile;
    if(f) ff=f;
    if(!ff) ::Fatal("TEXsec::GetElement","No file open!");
    fElements[fNLdElems] = (TEXsec *) ff->Get(TPartIndex::I()->EleSymb(z));
    if(!fElements[fNLdElems]) {
        ::Fatal("GetElement","Element z %d a %d not found",z,a);
        return 0; // just to make the compiler happy
    } else {
        // We loaded the element, however we have to see whether
        // the energy grid is the right one
        if(FloatDiff(TPartIndex::I()->Emin(),fElements[fNLdElems]->Emin(),1e-7) ||
           FloatDiff(TPartIndex::I()->Emax(),fElements[fNLdElems]->Emax(),1e-7) ||
           TPartIndex::I()->NEbins() != fElements[fNLdElems]->NEbins())
            // we have to resize the energy grid of the element
            fElements[fNLdElems]->Resample();
        return fElements[fNLdElems++];
    }
}

//___________________________________________________________________
Bool_t TEXsec::Prune()
{
    for(Int_t ip=0; ip<fNRpart; ++ip)
        fPXsec[ip].Prune();
    return kTRUE;
}

//___________________________________________________________________
Bool_t TEXsec::Resample()
{
    for(Int_t ip=0; ip<fNRpart; ++ip)
        fPXsec[ip].Resample();
    fEmin = TPartIndex::I()->Emin();
    fEmax = TPartIndex::I()->Emax();
    fNEbins = TPartIndex::I()->NEbins();
    fEGrid = TPartIndex::I()->EGrid();
    return kTRUE;
}

//___________________________________________________________________
void TEXsec::Draw(Option_t *option)  // mode=0->terminal, mode=1->viewer
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
    
    Float_t emin=TPartIndex::I()->Emin();
    if(narg>2) sscanf(((TObjString*) token->At(2))->GetName(),"%f",&emin);
    Float_t emax=TPartIndex::I()->Emax();
    if(narg>3) sscanf(((TObjString*) token->At(3))->GetName(),"%f",&emax);
    Int_t nbin=100;
    if(narg>4) sscanf(((TObjString*) token->At(4))->GetName(),"%d",&nbin);
    Int_t mode=0;
    if(narg>5) sscanf(((TObjString*) token->At(5))->GetName(),"%d",&mode);
    if(gFile) gFile->Get("PartIndex");
    Int_t pindex = TPartIndex::I()->PartIndex(part);
    if(pindex < 0) {
        Error("Draw","Unknown particle %s\n",part);
        return;
    } else if (pindex > TPartIndex::I()->NPartReac()) {
        Error("Draw","No reaction for particle %s\n",part);
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
    
    TCanvas *tc;
    if (mode){
        tc = fCanvas->GetCanvas();
        tc->SetLogx();
        tc->SetLogy();
        tc->SetGrid();
        same = kFALSE;
    }else{
        tc=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("G5canvas");
        //   if(!tc) tc = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
        if(!tc) {
            tc = new TCanvas("G5canvas",gtitle,600,400);
            tc->SetLogx();
            tc->SetLogy();
            tc->SetGrid();
            same = kFALSE;
            sec2.ReplaceAll("same","");
        }
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
    Double_t ymin=1e10;
    for(Int_t j=0; j<nreac; ++j) {
        const Char_t *reac = ((TObjString*) rnames->At(j))->GetName();
        TGraph *tg;
        if(TString(reac).BeginsWith("MS")) {
            tg = MSGraph(part,reac,emin,emax,nbin);
            for(Int_t i=0; i<tg->GetN(); ++i) {
                Double_t y=tg->GetY()[i];
                if(y>0 && y<ymin) ymin=y;
            }
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
            snprintf(ytitle,49,"GeV/cm / g/cm^3");
            tg = DEdxGraph(part, emin, emax, nbin);
            for(Int_t i=0; i<tg->GetN(); ++i) {
                Double_t y=tg->GetY()[i];
                if(y>0 && y<ymin) ymin=y;
            }
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
            for(Int_t i=0; i<tg->GetN(); ++i) {
                Double_t y=tg->GetY()[i];
                if(y>0 && y<ymin) ymin=y;
            }
            // a x-sec less than 1nb makes little sense...
            ymin = ymin<2e-9?2e-9:ymin;
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
    tmg->SetMinimum(0.5*ymin);
    tmg->Draw(gopt);
    if(!same) {
        ((TAxis*) tmg->GetHistogram()->GetXaxis())->SetTitle("Energy (GeV)");
        ((TAxis*) tmg->GetHistogram()->GetYaxis())->SetTitle(ytitle);
        ((TAxis*) tmg->GetHistogram()->GetXaxis())->CenterTitle();
        ((TAxis*) tmg->GetHistogram()->GetYaxis())->CenterTitle();
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
    tc->Update();
}

//___________________________________________________________________

//  All the methods below are used by the Viewer
void TEXsec::Viewer()
{
    if (fMain==NULL) {   // if the viewer ins't already active
        fMain = new TGMainFrame(gClient->GetRoot(),1050,520);
        fMain->Connect("CloseWindow()", "TEXsec", this, "ResetFrame()");
        fSecond = new TGHorizontalFrame(fMain, 200, 40);   //contains everything
        fCanvas = new TRootEmbeddedCanvas(0,fSecond,900, 520);   // contains the graphs
        
        TGVerticalFrame *fVertical = new TGVerticalFrame(fSecond, 200, 520);   // Create a vertical frame containing the lists & the buttons
        TGHorizontalFrame *fSelect = new TGHorizontalFrame(fVertical, 200, 40);   //contains the select all & deselect all buttons
        
        TGLabel *lParticle = new TGLabel(fVertical, "Particles : "); // label for the particle selector
        lParticle->SetTextFont("verdana");
        fVertical->AddFrame(lParticle, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2)); 
        
        fParticleBox = new TGListBox(fVertical);
        fParticleBox->Connect("Selected(Int_t)","TEXsec",this,"UpdateReactions()");   // when the user changes the particle with the mouse
        
        for (Int_t i=0; i<TPartIndex::I()->NPartReac(); i++) {
            fParticleBox->AddEntry(TPartIndex::I()->PartName(i), i);  // adding all the particles of the list
        }
        
        fParticleBox->SortByName();
        fParticleBox->Select(8);
        fParticleBox->Resize(150, 160);
        fVertical->AddFrame(fParticleBox,new TGLayoutHints(kLHintsCenterX, 5, 5, 5, 5));
        
        TGLabel *lSpace = new TGLabel(fVertical, "         ");
        fVertical->AddFrame(lSpace, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2));
        
        TGLabel *lReaction = new TGLabel(fVertical, "Reactions : "); // label for the reaction selector
        lReaction->SetTextFont("verdana");
        fVertical->AddFrame(lReaction, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2));
        
        fReactionBox = new TGListBox(fVertical, 90);
        fReactionBox->Connect("Selected(Int_t)","TEXsec",this,"PreDraw()");
        this->UpdateReactions();
        fReactionBox->Resize(150, 160);
        fReactionBox->SetMultipleSelections(kTRUE);
        fVertical->AddFrame(fReactionBox,new TGLayoutHints(kLHintsTop, 5, 5, 5, 5));
        
        TGTextButton *select = new TGTextButton(fSelect,"&Select All");      // adding the 3 buttons
        select->Connect("Clicked()","TEXsec",this,"SelectAll()");
        fSelect->AddFrame(select, new TGLayoutHints(kLHintsLeft, 0,5,3,4));
        TGTextButton *deselect = new TGTextButton(fSelect,"Deselect &All");
        deselect->Connect("Clicked()","TEXsec",this,"DeselectAll()");
        fSelect->AddFrame(deselect, new TGLayoutHints(kLHintsRight, 5,5,3,4));
        fVertical->AddFrame(fSelect, new TGLayoutHints(kLHintsLeft, 5, 5, 2, 2));
        
    
        TGTextButton *exit = new TGTextButton(fVertical,"&Exit", "gApplication->Terminate(0)");
        fVertical->AddFrame(exit, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5,5,3,4));
        
        fSecond->AddFrame(fVertical, new TGLayoutHints(kLHintsCenterY,2,2,2,2));
        fSecond->AddFrame(fCanvas, new TGLayoutHints(kLHintsRight |kLHintsExpandY | kLHintsExpandX,0,0,1,1));
        fMain->AddFrame(fSecond, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX, 0,0,1,1));
        
        TString name="ROOT XSec Viewer - ";
        name+=GetTitle();
        fMain->SetWindowName(name);
        fMain->MapSubwindows();
        fMain->MapWindow();
        fMain->Layout();
        
    }
    

    
}
//___________________________________________________________________
void TEXsec::SelectAll()
{
    
    for (Int_t i=0; i<fReactionBox->GetNumberOfEntries(); i++) {
        if(!fReactionBox->GetEntry(i)->IsActive()) {
            fReactionBox->GetEntry(i)->Toggle();
        }
    }
    fReactionBox->Layout();
    PreDraw();
}

//___________________________________________________________________
void TEXsec::DeselectAll()
{
    for (Int_t i=0; i<fReactionBox->GetNumberOfEntries(); i++) {
        if(fReactionBox->GetEntry(i)->IsActive()) {
            fReactionBox->GetEntry(i)->Toggle();
        }
    }
    fReactionBox->Layout();
    PreDraw();
}

//___________________________________________________________________
void TEXsec::UpdateReactions()   // when the user changes the selected particle
{
    fReactionBox->RemoveAll();
    if (fParticleBox->GetSelected()>-1) {  // if a particle is selected
        Int_t nb=0;
        for(Int_t i=0; i<TPartIndex::I()->NProc()-1; ++i) {
            if(XS(TPartIndex::I()->PartIndex(((TGTextLBEntry *)fParticleBox->GetSelectedEntry())->GetTitle()),i,0.0000000001)>=0) {
                fReactionBox->AddEntry(TPartIndex::I()->ProcName(i), nb);   // adding the possible reactions for the chosen particle
                fReactionBox->GetEntry(nb)->Toggle();
                nb++;
            }
        }
    }

    fReactionBox->Layout();
    
    PreDraw();
}

//___________________________________________________________________
void TEXsec::PreDraw()   // preparation of Draw() in the viewer, generation of the option
{
    if (fParticleBox->GetSelected()>-1) {  // if a particle is selected
        TString opt = ((TGTextLBEntry *)fParticleBox->GetSelectedEntry())->GetTitle();
        
        opt+=",Total|";
        for (Int_t i=0; i<fReactionBox->GetNumberOfEntries(); i++) {
            if(fReactionBox->GetEntry(i)->IsActive()) {
                opt+=((TGTextLBEntry *)fReactionBox->GetEntry(i))->GetTitle();
                opt+="|";
            }
        }
        opt+=",0.000000001,1000000000,100,1";
        Draw(opt);
    }
}

//___________________________________________________________________
void TEXsec::ResetFrame()
{
   delete fMain;
   fMain=NULL;
}



