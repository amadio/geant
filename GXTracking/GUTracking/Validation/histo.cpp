/*
=======
Author: M. Bandieramonte
 
 NB: 
 1. How to compile:
 g++ histo.cpp -o histo `root-config --cflags --glibs`

 2. Then execute the file and give the input interactively
     
     Once defined the projectile input energy (i.e. 100 MeV), the program is looking for the following input files:
     geant4_100MeV.root
     scalar_100MeV.root
     vector_100MeV.root
=======
*/

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <TVector3.h>
#include <cmath>
#include <unistd.h>
#include <TMarker.h>
#include <TGraph2D.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TAttText.h>

#include <stdio.h>
#include <string.h>

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TH3D.h"
#include "TF1.h"
#include "TImage.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <stdlib.h>
#include "TGaxis.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TLatex.h"

constexpr double electron_mass_c2 = 0.510998910 * 1.;
constexpr double inv_electron_mass_c2 = 1.0/electron_mass_c2;

using namespace ROOT::Math;

//_________
void MSaveBigPDF(double scale=5) {
    TCanvas* old_canv = gPad->GetCanvas();
    
    gROOT->SetBatch(kTRUE);
    gROOT->ForceStyle(kTRUE);
    
    Int_t orig_msz = gStyle->GetMarkerSize();
    Int_t orig_mst = gStyle->GetMarkerStyle();
    Int_t orig_lt  = gStyle->GetLineWidth();
    
    gStyle->SetMarkerSize(1.0+scale/5);
    gStyle->SetMarkerStyle(20);
    gStyle->SetLineWidth(orig_lt*scale);
   
    TString filename = old_canv->GetName();
    filename += ".png";
   
    
    Int_t old_width  = old_canv->GetWindowWidth();
    Int_t old_height = old_canv->GetWindowHeight();
    
    Int_t new_width = old_width * scale;
    Int_t new_height= old_height* scale;
    
    TCanvas* temp_canvas = new TCanvas("temp", "", new_width, new_height);
    old_canv->DrawClonePad();
    
    temp_canvas->SetLineWidth(orig_lt*40);
    temp_canvas->Draw();
    temp_canvas->SaveAs(filename);
    temp_canvas->Close();
    
    gStyle->SetMarkerSize(orig_msz);
    gStyle->SetMarkerStyle(orig_mst);
    gStyle->SetLineWidth(orig_lt);
    
    gROOT->ForceStyle(kFALSE);
    gROOT->SetBatch(kFALSE);
    
    std::cout<<"Saving the image as: "<<filename<<"\n";
    
    return;
}


//_________
void drawtext(double x, double y, const char *s)
{
    TLatex *t = new TLatex(x,y,Form("#chi^{2}: %s",s));
    t->SetTextFont(4);
    //t->SetTextAlign(12);
    t->SetNDC();
    t->SetTextColor(6);
    t->SetTextSize(0.048);
    t->Draw();
}


//______________________________________________________________________________
void validatePdf(TH1F* eOutScalar, TH1F* eOutVector, TH1F* eOutG4, char* energy)
{
    
    Int_t entries=eOutScalar->GetSize()-2; //NB: In Root the size of an histogram is nBin+2
    //eOutScalar->GetEntries();
    //std::cout<<"N of histogram entries: "<<eOutScalar->GetEntries()<<"\n";
    
    
    Double_t xScalar[entries], yScalar[entries], zScalar[entries];
    Double_t xVector[entries], yVector[entries], zVector[entries];
    Double_t xGeant4[entries], yGeant4[entries], zGeant4[entries];
    
    TString pdfFileName=Form("pdf_%sMeV.root",energy);
    TFile *fPfd = new TFile(pdfFileName,"w");
    TGraph *pdfGraph= (TGraph*)fPfd->Get("Pdf2.0"); //Get the first pdf, for fMinX value
    
    
    TString c2Name=Form("Pdf %s MeV",energy);
    
    TCanvas *c2 = new TCanvas("c2",c2Name,200,10,700,500);
    pdfGraph->SetLineColor(kBlue);
    //pdfGraph->SetTitle("pdf; mcs; dm");
    
    pdfGraph->Draw();
    c2->Update();
    
    
    //Scale the histograms
    double norm=1.;
    TH1F *eOutScalarScaled = (TH1F*)(eOutScalar->Clone("eOutScalarScaled"));
    norm= eOutScalarScaled->GetEntries();
    eOutScalarScaled->Scale(1/norm);
    
    
    TH1F *eOutVectorScaled = (TH1F*)(eOutVector->Clone("eOutVectorScaled"));
    norm= eOutVectorScaled->GetEntries();
    eOutVectorScaled->Scale(1/norm);
    
    
    TH1F *eOutGeant4Scaled = (TH1F*)(eOutG4->Clone("eOutGeant4Scaled"));
    norm= eOutGeant4Scaled->GetEntries();
    eOutGeant4Scaled->Scale(1/norm);

    
    for(int j = 0; j < entries ; ++j){
        yScalar[j] = eOutScalarScaled->GetBinContent(j+1);
        pdfGraph->GetPoint(j, xScalar[j],zScalar[j]); //to map
        
        yVector[j] = eOutVectorScaled->GetBinContent(j+1);
        pdfGraph->GetPoint(j, xVector[j],zVector[j]); //to map
        
        yGeant4[j] = eOutGeant4Scaled->GetBinContent(j+1);
        pdfGraph->GetPoint(j, xGeant4[j],zGeant4[j]); //to map
    }
    
    
    //Read the pdf File for 500MeV
    TCanvas *c3 = new TCanvas("c3","Graph comparison",200,10,700,500);
    pdfGraph->SetLineColor(kBlue);
    pdfGraph->Draw();
    c3->Update();
    
    
    
    //Create the Graph with x-values taken from the pdfGraph and the y-values taken from the simulation histogram
    TGraph *myCopyVectorGr=new TGraph(entries,xVector,yVector);
    //myCopyVectorGr->SetLineColor(kYellow+10);
    //myCopyVectorGr->SetLineStyle(2);
    //myCopyVectorGr->SetLineWidth(2);
    //myCopyVectorGr->Draw("LP");
    
    //TMarker *mark=new TMarker();
    //mark->SetMarkerStyle(2);
    //mark->SetMarkerSize(1);
    
    TGraph *myCopyScalarGr=new TGraph(entries,xScalar,yScalar);
    myCopyScalarGr->SetLineColor(kYellow+10);
    //myCopyScalarGr->SetLineStyle(4);
    //myCopyScalarGr->SetLineWidth(2);
    
    //myCopyScalarGr->SetMarkerStyle(20);
    myCopyScalarGr->Draw("LP");
    
    //TMarker *mark2=new TMarker();
    //mark2->SetMarkerStyle(4);
    
    
    TGraph *myCopyG4Gr=new TGraph(entries,xGeant4,yGeant4);
    myCopyG4Gr->SetLineColor(kRed);
    //myCopyG4Gr->SetLineStyle(5);
    //myCopyG4Gr->SetLineWidth(2);
    myCopyG4Gr->Draw("LP");
    
    TLegend *leg = new TLegend(0.4,0.6,0.9,0.9);
    leg->SetHeader("Legend");
    //leg->AddEntry(my,"Histogram filled with random numbers","f");
    //leg->AddEntry(my,"Function abs(#frac{sin(x)}{x})","l");
    //leg->AddEntry("gr","Graph with error bars","lep");
    TString legName=Form("Pdf for E_in=%s MeV",energy);
    leg->AddEntry(pdfGraph,legName,"l");
    leg->AddEntry(myCopyScalarGr,"E_out Histogram values scaled - Scalar","l");
    //leg->AddEntry(myCopyVectorGr,"E_out Histogram values scaled - Vector","l");
    leg->AddEntry(myCopyG4Gr,"E_out Histogram values scaled - Geant4","l");
    
    leg->Draw();
    c3->Update();
    c3->SaveAs("pdfValidation.pdf");
    
    
    ////Calculate chi-square
    std::cout<<"#E_out Histogram entries: "<<myCopyScalarGr->GetN()<<"#Pdf entries: "<<pdfGraph->GetN()<<"\n";
    double chiSquare_vector=0,chiSquare_scalar=0, chiSquare_geant4=0, expX,expY ;
    for (int i=0; i<pdfGraph->GetN(); i++)
    {
        pdfGraph->GetPoint(i, expX, expY);
        //myCopyScalarGr->GetPoint(i, obsX, obsY);
        //std::cout<<"expX: "<<expX<<" expY: "<<expY<<" obsX: "<<obsX<<" obsY: "<<obsY<<std::endl;
        chiSquare_scalar+=(yScalar[i]-expY)*(yScalar[i]-expY)/expY;
        chiSquare_vector+=(yVector[i]-expY)*(yVector[i]-expY)/expY;
        chiSquare_geant4+=(yGeant4[i]-expY)*(yGeant4[i]-expY)/expY;
    }
    std::cout<<"chiSquare scalar: "<<chiSquare_scalar<<" chiSquare-reduced: "<<chiSquare_scalar/(eOutScalarScaled->GetSize()-2)<<std::endl;
    std::cout<<"chiSquare vector: "<<chiSquare_vector<<" chiSquare-reduced: "<<chiSquare_vector/(eOutVectorScaled->GetSize()-2)<<std::endl;
    std::cout<<"chiSquare geant4: "<<chiSquare_geant4<<" chiSquare-reduced: "<<chiSquare_geant4/(eOutGeant4Scaled->GetSize()-2)<<std::endl;
    
}

//______________________________________________________________________________
void chiSquare_pdf(TH1F* eOutScalar, TH1F* eOutVector, TH1F* eOutG4, int energy)
{
    
    Int_t entries=eOutScalar->GetSize()-2; //NB: In Root the size of an histogram is nBin+2
    //eOutScalar->GetEntries();
    //std::cout<<"N of histogram entries: "<<eOutScalar->GetEntries()<<"\n";
    
    
    Double_t xScalar[entries], yScalar[entries], zScalar[entries];
    Double_t xVector[entries], yVector[entries], zVector[entries];
    Double_t xGeant4[entries], yGeant4[entries], zGeant4[entries];
    Double_t pdf[entries], x[entries];
    TGraph *pdfGraph;
    
    double logxmin = log(1);
    double dx = (log(10000) - logxmin)/99;
    
    //// pdf calculation
    double energy0 = exp(logxmin + dx*energy);
    
    double ymin = energy0/(1+2.0*energy0*inv_electron_mass_c2); //MeV
    double dy = (energy0 - ymin)/(1000);
    double yo = ymin + 0.5*dy;
    
    
    double sum = 0.;
    //double integralPdf=0;
    
    for(int j = 0; j < entries ; ++j) {
        //for each output energy bin
        double energy1 = yo + dy*j;
        
        double *grej=new double();
        
        
        
        double E0_m = energy0/0.510998910 ;
        double epsilon = energy1/energy0;
        
        double onecost = (1.- epsilon)/(epsilon*E0_m);
        double sint2   = onecost*(2.-onecost);
        double greject = 1. - epsilon*sint2/(1.+ epsilon*epsilon);
        double xsec = (epsilon + 1./epsilon)*greject;
        
        x[j]=energy1;
        pdf[j] = xsec;
        sum += xsec;
    }
    std::cout<<"xsec sum: "<<sum<<" 1/sum: "<<1/sum<<" \n";
    
    sum=1./sum;
    for(int j = 0; j < entries ; ++j)
        pdf[j]*=sum;
    
    pdfGraph = new TGraph(entries,x,pdf);
    
    TCanvas *c1 = new TCanvas("c1","pdf",200,10,700,500);
    pdfGraph->SetLineColor(kBlue);
    pdfGraph->Draw();
    c1->Update();

    
    
    
    //TString pdfFileName=Form("pdf_%sMeV.root",energy);
    //TFile *fPfd = new TFile(pdfFileName,"w");
    //TGraph *pdfGraph= (TGraph*)fPfd->Get("Pdf2.0"); //Get the first pdf, for fMinX value
    

    //Scale the histograms
    double norm=1.;
    TH1F *eOutScalarScaled = (TH1F*)(eOutScalar->Clone("eOutScalarScaled"));
    norm= eOutScalarScaled->GetEntries();
    eOutScalarScaled->Scale(1/norm);
    
    
    TH1F *eOutVectorScaled = (TH1F*)(eOutVector->Clone("eOutVectorScaled"));
    norm= eOutVectorScaled->GetEntries();
    eOutVectorScaled->Scale(1/norm);
    
    
    TH1F *eOutGeant4Scaled = (TH1F*)(eOutG4->Clone("eOutGeant4Scaled"));
    norm= eOutGeant4Scaled->GetEntries();
    eOutGeant4Scaled->Scale(1/norm);
    
    
    for(int j = 0; j < entries ; ++j){
        yScalar[j] = eOutScalarScaled->GetBinContent(j+1);
        //pdfGraph->GetPoint(j, xScalar[j],zScalar[j]); //to map
        
        yVector[j] = eOutVectorScaled->GetBinContent(j+1);
        //pdfGraph->GetPoint(j, xVector[j],zVector[j]); //to map
        
        yGeant4[j] = eOutGeant4Scaled->GetBinContent(j+1);
        //pdfGraph->GetPoint(j, xGeant4[j],zGeant4[j]); //to map
    }
    
    ////Calculate chi-square
    //std::cout<<"#E_out Histogram entries: "<<myCopyScalarGr->GetN()<<"\n";
    double chiSquare_vector=0,chiSquare_scalar=0, chiSquare_geant4=0, expX,expY ;
    for (int i=0; i<entries; i++)
    {
        chiSquare_scalar+=(yScalar[i]-pdf[i])*(yScalar[i]-pdf[i])/pdf[i];
        chiSquare_vector+=(yVector[i]-pdf[i])*(yVector[i]-pdf[i])/pdf[i];
        chiSquare_geant4+=(yGeant4[i]-pdf[i])*(yGeant4[i]-pdf[i])/pdf[i];
    }
    std::cout<<"chiSquare scalar: "<<chiSquare_scalar<<" chiSquare-reduced: "<<chiSquare_scalar/(eOutScalarScaled->GetSize()-2)<<std::endl;
    std::cout<<"chiSquare vector: "<<chiSquare_vector<<" chiSquare-reduced: "<<chiSquare_vector/(eOutVectorScaled->GetSize()-2)<<std::endl;
    std::cout<<"chiSquare geant4: "<<chiSquare_geant4<<" chiSquare-reduced: "<<chiSquare_geant4/(eOutGeant4Scaled->GetSize()-2)<<std::endl;
    
    //return pdfGraph;
}


//______________________________________________________________________________
double chiSquare(TH1F* eOutG4, TH1F* eOutScalar)
{
    
    int entries=eOutScalar->GetSize()-2; //NB: In Root the size of an histogram is nBin+2
    double chiSquare =0. ,chiSquare_scaled =0. , exp=0., obs=0. ;
 
    //calculate the chi-square of the scaled histo
    double norm=1.;
    TH1F *eOutScalarScaled = (TH1F*)(eOutScalar->Clone("eOutScalarScaled"));
    norm= eOutScalarScaled->GetEntries();
    eOutScalarScaled->Scale(1/norm);
    
    TH1F *eOutGeant4Scaled = (TH1F*)(eOutG4->Clone("eOutGeant4Scaled"));
    norm= eOutGeant4Scaled->GetEntries();
    eOutGeant4Scaled->Scale(1/norm);

    
    for (int i=0; i<entries; i++)
    {
        exp = eOutG4->GetBinContent(i); //to verify
        obs = eOutScalar->GetBinContent(i);
        
        //std::cout<<"exp:  "<<exp<<" obs:  "<<obs<<std::endl;
        if(exp!=0)
        chiSquare+=((obs-exp)*(obs-exp)/exp);
        
        exp = eOutGeant4Scaled->GetBinContent(i); //to verify
        obs = eOutScalarScaled->GetBinContent(i);
        if(exp!=0)
            chiSquare_scaled+=((obs-exp)*(obs-exp)/exp);
        
    }
    std::cout<<"chiSquare "<<chiSquare<<std::endl;
    std::cout<<"chiSquare_scaled "<<chiSquare_scaled<<std::endl;
    
    return chiSquare_scaled;
}

//______________________________________________________________________________
void genAllHisto(const char *process, const char* energy)
{
    
    TString processName= process;
    //TString histogramName= variable;
    TString geant4RootFileName= Form("geant4_%sMeV.root",energy);
    TString scalarRootFileName= Form("scalar_%sMeV.root",energy);
    TString vectorRootFileName= Form("vector_%sMeV.root",energy);
    
    TString g4HistoName=Form("%s/geant4",process);
    TString gvHistoScalarName=Form("%s/geantVscalar",process);
    TString gvHistoVectorName=Form("%s/geantVvector",process);
    
    TString histPath= Form("%s/%s",process, "EnergyOut1");
    TString histName= Form("%s/%s/%sMeV",process, "EnergyOut1", energy);
    //TString histoFileNamePdf= Form("%s-%s-%sMeV.pdf",process, variable, energy);
    TString histoFileNameEps= Form("%s-%sMeV.eps",process, energy);
    
    //Geant 4
    TFile *fG4 = new TFile(geant4RootFileName,"w");
    TH1F *eOutG4 = (TH1F*)fG4->Get(histPath);
    if(eOutG4) eOutG4->SetName(" geant4 ");
    else std::cout<<"eOutG4 is null\n";
    
    //GeantV scalar
    TFile *fScalar = new TFile(scalarRootFileName,"w");
    TH1F *eOutScalar = (TH1F*)fScalar->Get(histPath);
    if(eOutScalar) eOutScalar->SetName("geantVscalar");
    else std::cout<<"eOutScalar is null\n";
    
    //GeantV vector
    TFile *fVector = new TFile(vectorRootFileName,"w");
    TH1F *eOutVector = (TH1F*)fVector->Get(histPath);
    if(eOutVector) eOutVector->SetName("geantVvector");
    else std::cout<<"eOutVector is null\n";
    
    TCanvas *MyC1 = new TCanvas(process,process,800,200,1200,1000);
    //TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);

    MyC1->Divide(2,2);
    
    //1
    MyC1->cd(1);
    
    

    
    if(!strcmp(process,"KleinNishina"))
    {
        eOutG4->GetXaxis()->SetTitle("Scattered Photon Energy [MeV]");
        eOutG4->GetYaxis()->SetTitle("dN/dE");
    } else
    {
        eOutG4->GetXaxis()->SetTitle("EnergyOut1");
        eOutG4->GetYaxis()->SetTitle("dN/dE");
    }
    
    gStyle->SetEndErrorSize(3);
    gStyle->SetErrorX(1.);
    //eOutG4->SetMarkerStyle(20);
    //he->Draw("E1");
    
    eOutG4->Draw("E");
    MyC1->Update();
    
    TPaveStats *ps1 = (TPaveStats*)eOutG4->GetListOfFunctions()->FindObject("stats");
    ps1->SetX1NDC(0.4); ps1->SetX2NDC(0.55);
    MyC1->Modified();
    MyC1->Update();
    eOutScalar->SetLineColor(kYellow+10);
    //eOutScalar->Draw("][sames");
    eOutScalar->Draw("sames");
    
    MyC1->Update();
    
    TPaveStats *ps2 = (TPaveStats*)eOutScalar->GetListOfFunctions()->FindObject("stats");
    ps2->SetX1NDC(0.6); ps2->SetX2NDC(0.75);
    ps2->SetTextColor(kYellow+10);
    ps2->SetOptStat(1111);
    
    eOutVector->SetLineColor(kMagenta);
    //eOutVector->Draw("][sames");
    eOutVector->Draw("sames");
    
    
    double chiS=chiSquare(eOutG4,eOutScalar);
    TString chiSquareString= Form("%f",chiS);
    drawtext(.1, 0.92, chiSquareString);
    
    MyC1->Update();
    TPaveStats *ps3 = (TPaveStats*)eOutVector->GetListOfFunctions()->FindObject("stats");
    ps3->SetX1NDC(0.8); ps3->SetX2NDC(0.95);
    ps3->SetTextColor(kMagenta);
    TPaveText *t = new TPaveText(0.0, 0.95, 0.3, 1.0, "brNDC"); // left-up
    //t->AddText(histName);
    

    //2
    MyC1->cd(2);
    histPath= Form("%s/%s",process, "EnergyOut2");
    histName= Form("%s/%s/%sMeV",process, "EnergyOut2", energy);
    eOutG4 = (TH1F*)fG4->Get(histPath);
    eOutScalar = (TH1F*)fScalar->Get(histPath);
    eOutVector = (TH1F*)fVector->Get(histPath);
    eOutG4->SetName(" geant4 ");
    eOutScalar->SetName("geantVscalar");
    eOutVector->SetName("geantVvector");
    eOutScalar->SetLineColor(kYellow+10);
    eOutVector->SetLineColor(kMagenta);

    
    if(!strcmp(process,"KleinNishina"))
    {
        eOutG4->GetXaxis()->SetTitle("Electron Energy [MeV]");
        eOutG4->GetYaxis()->SetTitle("dN/dE");
    } else
    {
        eOutG4->GetXaxis()->SetTitle("EnergyOut2");
        eOutG4->GetYaxis()->SetTitle("dN/dE");
    }
    
    
    eOutG4->Draw("E");
    eOutScalar->Draw("sames");
    eOutVector->Draw("sames");
    MyC1->Update();
    TPaveStats *psEnOut1 = (TPaveStats*)eOutG4->GetListOfFunctions()->FindObject("stats");
    psEnOut1->SetX1NDC(0.4); psEnOut1->SetX2NDC(0.55);
    //psEnOut1->SetTextColor(kMagenta);

    TPaveStats *psEnOut2 = (TPaveStats*)eOutScalar->GetListOfFunctions()->FindObject("stats");
    psEnOut2->SetX1NDC(0.6); psEnOut2->SetX2NDC(0.75);
    psEnOut2->SetTextColor(kYellow+10);
    psEnOut2->SetOptStat(1111);

    chiS=chiSquare(eOutG4,eOutScalar);
    chiSquareString= Form("%f",chiS);
    drawtext(.1, 0.92, chiSquareString);
    
    TPaveStats *psEnOut3 = (TPaveStats*)eOutVector->GetListOfFunctions()->FindObject("stats");
    psEnOut3->SetX1NDC(0.8); psEnOut3->SetX2NDC(0.95);
    psEnOut3->SetTextColor(kMagenta);
    
    
    //3
    MyC1->cd(3);
    histPath= Form("%s/%s",process, "AngleOut1");
    histName= Form("%s/%s/%sMeV",process, "AngleOut1", energy);
    eOutG4 = (TH1F*)fG4->Get(histPath);
    eOutScalar = (TH1F*)fScalar->Get(histPath);
    eOutVector = (TH1F*)fVector->Get(histPath);
    eOutG4->SetName(" geant4 ");
    eOutScalar->SetName("geantVscalar");
    eOutVector->SetName("geantVvector");
    
    if(!strcmp(process,"KleinNishina"))
    {
        eOutG4->GetXaxis()->SetTitle("Scattered Photon Angle");
        //eOutG4->GetYaxis()->SetTitle("dN/dE");
    } else
    {
        eOutG4->GetXaxis()->SetTitle("AngleOut1");
        //eOutG4->GetYaxis()->SetTitle("dN/dE");
    }
    
    //eOutG4->GetXaxis()->SetTitle("AngleOut1");
    
    eOutG4->Draw("E");
    eOutScalar->SetLineColor(kYellow+10);
    eOutVector->SetLineColor(kMagenta);
    eOutScalar->Draw("sames");
    eOutVector->Draw("sames");
    MyC1->Update();
    
    chiS=chiSquare(eOutG4,eOutScalar);
    chiSquareString= Form("%f",chiS);
    drawtext(.1, 0.92, chiSquareString);
    
    TPaveStats *psAOut1 = (TPaveStats*)eOutG4->GetListOfFunctions()->FindObject("stats");
    psAOut1->SetX1NDC(0.4); psAOut1->SetX2NDC(0.55);

    TPaveStats *psAOut2 = (TPaveStats*)eOutScalar->GetListOfFunctions()->FindObject("stats");
    psAOut2->SetX1NDC(0.6); psAOut2->SetX2NDC(0.75);
    psAOut2->SetTextColor(kYellow+10);
    psAOut2->SetOptStat(1111);
    
    TPaveStats *psAOut3 = (TPaveStats*)eOutVector->GetListOfFunctions()->FindObject("stats");
    psAOut3->SetX1NDC(0.8); psAOut3->SetX2NDC(0.95);
    psAOut3->SetTextColor(kMagenta);
    
    //4
    MyC1->cd(4);
    histPath= Form("%s/%s",process, "AngleOut2");
    histName= Form("%s/%s/%sMeV",process, "AngleOut2", energy);
    eOutG4 = (TH1F*)fG4->Get(histPath);
    eOutScalar = (TH1F*)fScalar->Get(histPath);
    eOutVector = (TH1F*)fVector->Get(histPath);
    eOutG4->SetName(" geant4 ");
    eOutScalar->SetName("geantVscalar");
    eOutVector->SetName("geantVvector");
    
    

    if(!strcmp(process,"KleinNishina"))
    {
        eOutG4->GetXaxis()->SetTitle("Electron Angle");
        //eOutG4->GetYaxis()->SetTitle("dN/dE");
    } else
    {
        eOutG4->GetXaxis()->SetTitle("AngleOut2");
        //eOutG4->GetYaxis()->SetTitle("dN/dE");
    }
    //eOutG4->GetXaxis()->SetTitle("AngleOut2");
    eOutG4->Draw("E");
    eOutScalar->SetLineColor(kYellow+10);
    eOutVector->SetLineColor(kMagenta);
    eOutScalar->Draw("sames");
    eOutVector->Draw("sames");
    MyC1->Update();
    
    chiS=chiSquare(eOutG4,eOutScalar);
    chiSquareString= Form("%f",chiS);
    drawtext(.1, 0.92, chiSquareString);
    
    TPaveStats *psAngleOut1 = (TPaveStats*)eOutG4->GetListOfFunctions()->FindObject("stats");
    psAngleOut1->SetX1NDC(0.4); psAngleOut1->SetX2NDC(0.55);
    
    TPaveStats *psAngleOut2 = (TPaveStats*)eOutScalar->GetListOfFunctions()->FindObject("stats");
    psAngleOut2->SetX1NDC(0.6); psAngleOut2->SetX2NDC(0.75);
    psAngleOut2->SetTextColor(kYellow+10);
    psAngleOut2->SetOptStat(1111);
    
    TPaveStats *psAngleOut3 = (TPaveStats*)eOutVector->GetListOfFunctions()->FindObject("stats");
    psAngleOut3->SetX1NDC(0.8); psAngleOut3->SetX2NDC(0.95);
    psAngleOut3->SetTextColor(kMagenta);

    //MSaveBigPDF();
    MyC1->SaveAs(histoFileNameEps);

}


//______________________________________________________________________________
void genSingleHisto(const char *process, const char *variable, const char* energy)
{
    TString processName= process;
    TString histogramName= variable;
    
    TString histPath= Form("%s/%s",process, variable);
    TString histName= Form("%s/%s/%sMeV",process, variable, energy);
    TString histoFileNamePdf= Form("%s-%s-%sMeV.pdf",process, variable, energy);
    TString histoFileNameEps= Form("%s-%s-%sMeV.eps",process, variable, energy);
    
    TString g4HistoName=Form("%s/geant4",process);
    TString gvHistoScalarName=Form("%s/geantVscalar",process);
    TString gvHistoVectorName=Form("%s/geantVvector",process);
    
    
    TString geant4RootFileName= Form("geant4_%sMeV.root",energy);
    TString scalarRootFileName= Form("scalar_%sMeV.root",energy);
    TString vectorRootFileName= Form("vector_%sMeV.root",energy);
    //TString pdfFileName=Form("pdf_%sMeV.root",energy);
    
    //Geant 4
    TFile *fG4 = new TFile(geant4RootFileName,"w");
    //fG4->ls();
    //std::cout<<"Path to the histogram:"<< histPath<<std::endl;
    TH1F *eOutG4 = (TH1F*)fG4->Get(histPath);
    if(eOutG4) eOutG4->SetName("geant4");
    else std::cout<<"eOutG4 is null\n";
    
    //GeantV scalar
    TFile *fScalar = new TFile(scalarRootFileName,"w");
    TH1F *eOutScalar = (TH1F*)fScalar->Get(histPath);
    if(eOutScalar) eOutScalar->SetName("geantVscalar");
    else std::cout<<"eOutScalar is null\n";
    
    //GeantV vector
    TFile *fVector = new TFile(vectorRootFileName,"w");
    TH1F *eOutVector = (TH1F*)fVector->Get(histPath);
    if(eOutVector) eOutVector->SetName("geantVvector");
    else std::cout<<"eOutVector is null\n";
    
    //int intenergy= atoi(energy);
    //TGraph my=chiSquare_pdf(eOutG4, eOutScalar, eOutVector, intenergy);
    double chiS= chiSquare(eOutG4, eOutScalar);
    
    std::cout<<"chiS: "<<chiS<<"\n";
    
    TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);
    

    if(!strcmp(process,"KleinNishina"))
       {
           std::cout<<"esatto\n";
           if(!strcmp(variable,"EnergyIn"))
           {
               eOutG4->GetXaxis()->SetTitle("EnergyIn [MeV]");
               eOutG4->GetYaxis()->SetTitle("dN/dE");
               
           }else
               if(!strcmp(variable,"EnergyOut1"))
               {
                   eOutG4->GetXaxis()->SetTitle("Scattered Photon Energy [MeV]");
                   eOutG4->GetYaxis()->SetTitle("dN/dE");
                   
               }else
                   if(!strcmp(variable,"EnergyOut2"))
                   {
                       eOutG4->GetXaxis()->SetTitle("Electron Energy [MeV]");
                       eOutG4->GetYaxis()->SetTitle("dN/dE");
                       
                   }else
                       //angle of the scatterred photon
                       if(!strcmp(variable,"AngleOut1"))
                       {
                           eOutG4->GetXaxis()->SetTitle("Scattered Photon Angle");
                           //eOutG4->GetYaxis()->SetTitle("dN/dsin(Theta)");
                           
                       }else
                           if(!strcmp(variable,"AngleOut2"))
                           {
                               eOutG4->GetXaxis()->SetTitle("Electron Angle");
                               //eOutG4->GetYaxis()->SetTitle("dN/dsin(Theta)");
                           }
       } else
       {
           std::cout<<"sbagliato\n";
           if(!strcmp(variable,"EnergyIn"))
           {
               eOutG4->GetXaxis()->SetTitle("EnergyIn [MeV]");
               eOutG4->GetYaxis()->SetTitle("dN/dE");
               
           }else
               if(!strcmp(variable,"EnergyOut1"))
               {
                   eOutG4->GetXaxis()->SetTitle("EnergyOut1 [MeV]");
                   eOutG4->GetYaxis()->SetTitle("dN/dE");
                   
               }else
                   if(!strcmp(variable,"EnergyOut2"))
                   {
                       eOutG4->GetXaxis()->SetTitle("EnergyOut2 [MeV]");
                       eOutG4->GetYaxis()->SetTitle("dN/dE");
                       
                   }else
                       
                       if(!strcmp(variable,"AngleOut1"))
                       {
                           eOutG4->GetXaxis()->SetTitle("AngleOut1");
                           //eOutG4->GetYaxis()->SetTitle("dN/dsin(Theta)");
                           
                       }else
                           if(!strcmp(variable,"AngleOut2"))
                           {
                               eOutG4->GetXaxis()->SetTitle("AngleOut2");
                               //eOutG4->GetYaxis()->SetTitle("dN/dsin(Theta)");
                           }

       
       }
    
    eOutG4->Draw("E");
    //eOutG4->SetLineColor(kGreen);
    c1->Update();

    
    TPaveStats *ps1 = (TPaveStats*)eOutG4->GetListOfFunctions()->FindObject("stats");
    //ps1->SetTextSize(0.005);
    ps1->SetX1NDC(0.4); ps1->SetX2NDC(0.55);
    //ps1->SetX1NDC(0.6); ps1->SetX2NDC(0.75);
    
    c1->Modified();
    c1->Update();
    
    eOutScalar->SetLineColor(kYellow+10);
    eOutScalar->Draw("sames");
    
    c1->Update();
    
    TPaveStats *ps2 = (TPaveStats*)eOutScalar->GetListOfFunctions()->FindObject("stats");
    ps2->SetX1NDC(0.6); ps2->SetX2NDC(0.75);
    //ps2->SetX1NDC(0.8); ps2->SetX2NDC(0.95);
    ps2->SetTextColor(kYellow+10);
    ps2->SetOptStat(1111);
    
    
    eOutVector->SetLineColor(kMagenta);
    eOutVector->Draw("sames");
    
    c1->Update();
    
    
    TPaveStats *ps3 = (TPaveStats*)eOutVector->GetListOfFunctions()->FindObject("stats");
    ps3->SetX1NDC(0.8); ps3->SetX2NDC(0.95);
    ps3->SetTextColor(kMagenta);
    
    TPaveText *t = new TPaveText(0.0, 0.95, 0.3, 1.0, "brNDC"); // left-up
    t->AddText(histName);
    
    gROOT->SetStyle("Plain");
    //gStyle->SetOptStat(0000);
    gStyle->SetOptFit(1111);
    gStyle->SetOptTitle(0);
    t->SetBorderSize(0);
    t->SetFillColor(gStyle->GetTitleFillColor());
    t->Draw();
    
    TString chiSquareString= Form("%f",chiS);
    drawtext(.1, 0.92, chiSquareString);
    
    c1->Modified();
    c1->Update();

    c1->SaveAs(histoFileNameEps);
    
    //validatePdf(eOutScalar, eOutVector, eOutG4, energy);
    /*c1->Destructor();
    if(c1!=NULL)
    {
        c1->Clear();
        c1->Closed();
        c1 = NULL;
    }*/
}

//______________________________________________________________________________
void scaleHisto(TH1F* eOutScalar, TH1F* eOutVector, TH1F* eOutG4)
{
    
    //Scale the histograms and plot them
    TCanvas *c0 = new TCanvas("c0","Scaled histograms",200,10,700,500);
    
    double norm=1.;
    TH1F *eOutScalarScaled = (TH1F*)(eOutScalar->Clone("eOutScalarScaled"));
    norm= eOutScalarScaled->GetEntries();
    eOutScalarScaled->Scale(1/norm);
    eOutScalarScaled->SetLineColor(kMagenta);
    eOutScalarScaled->Draw("][sames");
    c0->Update();
    
    
    TH1F *eOutVectorScaled = (TH1F*)(eOutVector->Clone("eOutVectorScaled"));
    norm= eOutVectorScaled->GetEntries();
    eOutVectorScaled->Scale(1/norm);
    eOutVectorScaled->SetLineColor(kBlue);
    eOutVectorScaled->Draw("][sames");
    c0->Update();
    
    TH1F *eOutGeant4Scaled = (TH1F*)(eOutG4->Clone("eOutGeant4Scaled"));
    norm= eOutGeant4Scaled->GetEntries();
    eOutGeant4Scaled->Scale(1/norm);
    eOutGeant4Scaled->SetLineColor(kGreen+2);
    eOutGeant4Scaled->Draw("][sames");
    c0->Update();
    
}





//______________________________________________________________________________
int main( int argc,  char *argv[])
{
	TApplication theApp("App",&argc,argv);
    
    
    const char* GUPhysicsModelName[5] = {
        "KleinNishina",  //Compton - gamma
        "BetheHeitler",  //Conversion - pair production - gamma
        "SauterGavrila", // Photo-Electric Effect - gamma
        "MollerBhabha",  // Ionization -electron
        "SeltzerBerger"  // Bremsstrahlung
    };
    
    const char* GUPhysicsModelVariable[5] = {
        "EnergyIn",
        "EnergyOut1",
        "EnergyOut2",
        "AngleOut1",
        "AngleOut2"
    };
    

    TString modelName;
    TString variableName;
    TString  in1, in2, in3;
    const char* phisicsModel, *variable, *energy;
    TString myInput;
    std::cout<<"Starting.\n";
    
    do{
    
        
        std::cout<<"Press 'D' for the default execution, press 'C' for the customized execution\n";
    
        std::cin>>myInput;
        if(myInput.EqualTo("D")|| myInput.EqualTo("d"))
        {
            std::cout<<"Insert the energy of the projectile [Mev]\n";
            std::cin>>in3;
            energy=in3.Data();
            std::cout << "Starting default mode.\n";
            for(int i=0; i<5 ; i++)
                genAllHisto(GUPhysicsModelName[i], energy);
            //std::cout << "Run completed.\n";
            //theApp.Run();
            //return 0;
            
        }
        else
        {
            std::cout<<"PhysicsModel: KleinNishina, BetheHeitler, SauterGavrila, MollerBhabha, SeltzerBerger.\n";
            std::cin>>in1;
            phisicsModel=in1.Data();
            std::cout<<"Variable: EnergyIn, EnergyOut1, EnergyOut2, AngleOut1, AngleOut2, or all.\n";
            std::cin>>in2;
            variable=in2.Data();
            
            std::cout<<"Energy of the projectile [Mev]\n";
            std::cin>>in3;
            energy=in3.Data();
            if(in2.EqualTo("all")||in2.EqualTo("All"))
                genAllHisto(phisicsModel,energy);
            else
                genSingleHisto(phisicsModel, variable, energy);
        }
        std::cout<<"Press 'c' to continue, 'e' to exit\n";
        std::cin>>myInput;
    } while(myInput.EqualTo("c")|| myInput.EqualTo("C"));
    
    
    /*
    if (argc<2)
    {
        std::cout << "Usage " << argv[0] << "\n <PHYSICSMODEL_NAME> \n: KleinNishina  \n: BetheHeitler \n: SauterGavrila \n: MollerBhabha \n: SeltzerBerger\n \n <HISTOGRAM_NAME> \n: EnergyIn \n: EnergyOut1 \n: EnergyOut2 \n: AngleOut1 \n: AngleOut2\n \n<ENERGY> [MeV] \n";
        std::cout << "Default mode usage: "<< argv[0] <<" <ENERGY> [MeV] \n";
        return 0;
    }
    
    int argIndex = 1;
    bool default=true;
    
    while (argIndex < argc){
        
        if (!strcmp(argv[argIndex],"-default")) {
            ++argIndex;
            
            ++argIndex;
        }

        else if (!strcmp(argv[argIndex],"-model")) {
            ++argIndex;
            modelName= theApp.Argv(argIndex);
            ++argIndex;
        }
        
        else if (!strcmp(argv[argIndex],"-histo")) {
            ++argIndex;
            variableName= theApp.Argv(argIndex);
            ++argIndex;
        }
        else if (!strcmp(argv[argIndex],"-energy")) {
            ++argIndex;
            energy= theApp.Argv(argIndex);
            ++argIndex;
        }
        
    }
    
    
    
    if(argc==2)
    {
        std::cout << "Starting default mode.\n";
        for(int i=0; i<5 ; i++)
            for(int j=0; j<5; j++)
            {
                genSingleHisto(GUPhysicsModelName[i], GUPhysicsModelVariable[j], theApp.Argv(1), i+j);
                std::cout << "Faccio altro, magari .\n";
            }
        std::cout << "Run completed.\n";
        theApp.Run();
        return 0;
    }
    
    
    else if (argc<4)
    {
        std::cout << "Usage " << argv[0] << "\n <PHYSICSMODEL_NAME> \n: KleinNishina  \n: BetheHeitler \n: SauterGavrila \n: MollerBhabha \n: SeltzerBerger\n \n <HISTOGRAM_NAME> \n: EnergyIn \n: EnergyOut1 \n: EnergyOut2 \n: AngleOut1 \n: AngleOut2\n \n<ENERGY> [MeV] \n";
        return 0;
    }
    

    else
    if( !(strcmp(argv[1], "KleinNishina")==0) && !(strcmp(argv[1], "BetheHeitler")==0) && !(strcmp(argv[1], "SauterGravila")==0) && !(strcmp(argv[1], "MollerBhabha")==0) && !(strcmp(argv[1], "SeltzerBerger")==0))
    {
        std::cout<<"Wrong PHYSICSMODEL_NAME.\n";
        std::cout << "Usage " << argv[0] << "\n <PHYSICSMODEL_NAME> \n: KleinNishina  \n: BetheHeitler \n: SauterGavrila \n: MollerBhabha \n: SeltzerBerger\n \n <HISTOGRAM_NAME> \n: EnergyIn \n: EnergyOut1 \n: EnergyOut2 \n: AngleOut1 \n: AngleOut2\n\n <ENERGY> [MeV] \n";
        return 0;
    }
   
    else if( !(strcmp(argv[2], "EnergyIn")==0) && !(strcmp(argv[2], "EnergyOut1")==0) && !(strcmp(argv[2], "EnergyOut2")==0) && !(strcmp(argv[2], "AngleOut1")==0) && !(strcmp(argv[2], "AngleOut2")==0))
    {
        std::cout<<"Wrong HISTOGRAM_NAME.\n";
        std::cout << "Usage " << argv[0] << "\n <PHYSICSMODEL_NAME> \n: KleinNishina  \n: BetheHeitler \n: SauterGavrila \n: MollerBhabha \n: SeltzerBerger\n \n <HISTOGRAM_NAME> \n: EnergyIn \n: EnergyOut1 \n: EnergyOut2 \n: AngleOut1 \n: AngleOut2\n \n <ENERGY> [MeV] \n";
        return 0;
    }
    
    std::cout << "Starting.\n";
    genSingleHisto(theApp.Argv(1), theApp.Argv(2), theApp.Argv(3),1);
    
    TCanvas * newCanvas=new TCanvas("c1","c",200,10,700,500);

    double y = 0.95;
    for (int f = 12; f<=152; f+=10) {
        if (f!=142) drawtext(0.02,y, f,"ABCDEFGH abcdefgh 0123456789 @#$");
        else drawtext(0.02,y, f,"ABCD efgh 01234 @#$");
        y -= 0.065;
    }
    genSingleHisto("KleinNishina", "AngleOut2", "500",6);*/
    

    /*
    TString processName= theApp.Argv(1);
    TString histogramName= theApp.Argv(2);
 
    TString histPath= Form("%s/%s",theApp.Argv(1), theApp.Argv(2));
    TString histName= Form("%s/%s/%sMeV",theApp.Argv(1), theApp.Argv(2), theApp.Argv(3));
    TString histoFileName= Form("%s-%s-%sMeV.pdf",theApp.Argv(1), theApp.Argv(2), theApp.Argv(3));
    
    TString g4HistoName=Form("%s/geant4",theApp.Argv(1));
    TString gvHistoScalarName=Form("%s/geantVscalar",theApp.Argv(1));
    TString gvHistoVectorName=Form("%s/geantVvector",theApp.Argv(1));
    
    TString geant4RootFileName= Form("geant4_%sMeV.root",theApp.Argv(3));
    TString scalarRootFileName= Form("scalar_%sMeV.root",theApp.Argv(3));
    TString vectorRootFileName= Form("vector_%sMeV.root",theApp.Argv(3));
    TString pdfFileName=Form("pdf_%sMeV.root",theApp.Argv(3));
    
    //Geant 4
    
    TFile *fG4 = new TFile(geant4RootFileName,"w");
    fG4->ls();
    std::cout<<"Path to the histogram:"<< histPath<<std::endl;
    //TH1F *eOutG4 = (TH1F*)fG4->Get("KleinNishina /EnergyOut1");
    //eOutG4->SetName("KleinNishina/geant4");
    
    TH1F *eOutG4 = (TH1F*)fG4->Get(histPath);
    if(eOutG4) eOutG4->SetName("geant4");
    else std::cout<<"eOutG4 is null\n";
    
    
    //GeantV scalar
    TFile *fScalar = new TFile(scalarRootFileName,"w");
    //TH1F *eOutScalar = (TH1F*)fScalar->Get("KleinNishina /EnergyOut1");
    //eOutScalar->SetName("KleinNishina/scalar");
    
    TH1F *eOutScalar = (TH1F*)fScalar->Get(histPath);
    
    if(eOutScalar) eOutScalar->SetName("geantVscalar");
    else std::cout<<"eOutScalar is null\n";
    
    
    //GeantV vector
    TFile *fVector = new TFile(vectorRootFileName,"w");
    //TH1F *eOutScalar = (TH1F*)fScalar->Get("KleinNishina /EnergyOut1");
    //eOutScalar->SetName("KleinNishina/scalar");
    
    TH1F *eOutVector = (TH1F*)fVector->Get(histPath);
    if(eOutVector) eOutVector->SetName("geantVvector");
    else std::cout<<"eOutVector is null\n";
    
    */
    
    
    
    

    /*
    //Plot the not scaled histograms
    TCanvas *c1 = new TCanvas("c1","Physics validation",200,10,700,500);
    eOutG4->Draw();
    c1->Update();
    
    TPaveStats *ps1 = (TPaveStats*)eOutG4->GetListOfFunctions()->FindObject("stats");
    ps1->SetX1NDC(0.4); ps1->SetX2NDC(0.55);
    //ps1->SetX1NDC(0.5); ps1->SetX2NDC(0.7);
    
    c1->Modified();
    c1->Update();

    eOutScalar->SetLineColor(kYellow+10);
    eOutScalar->Draw("][sames");
    c1->Update();
    
    TPaveStats *ps2 = (TPaveStats*)eOutScalar->GetListOfFunctions()->FindObject("stats");
    ps2->SetX1NDC(0.6); ps2->SetX2NDC(0.75);
    //ps2->SetX1NDC(0.75); ps2->SetX2NDC(0.95);
    ps2->SetTextColor(kYellow+10);
    
    eOutVector->SetLineColor(kMagenta);
    eOutVector->Draw("][sames");
    c1->Update();
    
    TPaveStats *ps3 = (TPaveStats*)eOutVector->GetListOfFunctions()->FindObject("stats");
    ps3->SetX1NDC(0.8); ps3->SetX2NDC(0.95);
    ps3->SetTextColor(kMagenta);

    TPaveText *t = new TPaveText(0.0, 0.95, 0.3, 1.0, "brNDC"); // left-up
    t->AddText(histName);
    
    gROOT->SetStyle("Plain");
    //gStyle->SetOptStat(0000);
    gStyle->SetOptFit(1111);
    gStyle->SetOptTitle(0);
    t->SetBorderSize(0);
    t->SetFillColor(gStyle->GetTitleFillColor());
    t->Draw();
    
    c1->Modified();
    c1->Update();
    */
    
    
    
    /*
    
    //PDF validation
    Int_t entries=eOutScalar->GetSize()-2; //NB: In Root the size of an histogram is nBin+2
    //eOutScalar->GetEntries();
    //std::cout<<"N of histogram entries: "<<eOutScalar->GetEntries()<<"\n";
    
    
    Double_t xScalar[entries], yScalar[entries], zScalar[entries];
    Double_t xVector[entries], yVector[entries], zVector[entries];
    Double_t xGeant4[entries], yGeant4[entries], zGeant4[entries];
    
    TFile *fPfd = new TFile(pdfFileName,"w");
    TGraph *pdfGraph= (TGraph*)fPfd->Get("Pdf2.0"); //Get the first pdf, for fMinX value
    
    
    TString c2Name=Form("Pdf %s MeV",theApp.Argv(3));
    
    TCanvas *c2 = new TCanvas("c2",c2Name,200,10,700,500);
    pdfGraph->SetLineColor(kBlue);
    //pdfGraph->SetTitle("pdf; mcs; dm");
 
    pdfGraph->Draw();
    c2->Update();

    
    
    for(int j = 0; j < entries ; ++j){
        yScalar[j] = eOutScalarScaled->GetBinContent(j+1);
        pdfGraph->GetPoint(j, xScalar[j],zScalar[j]); //to map
        
        yVector[j] = eOutVectorScaled->GetBinContent(j+1);
        pdfGraph->GetPoint(j, xVector[j],zVector[j]); //to map
        
        yGeant4[j] = eOutGeant4Scaled->GetBinContent(j+1);
        pdfGraph->GetPoint(j, xGeant4[j],zGeant4[j]); //to map
     }
    
    
    //Read the pdf File for 500MeV
    TCanvas *c3 = new TCanvas("c3","Graph comparison",200,10,700,500);
    pdfGraph->SetLineColor(kBlue);
    pdfGraph->Draw();
    c3->Update();

    
    
    //Create the Graph with x-values taken from the pdfGraph and the y-values taken from the simulation histogram
    TGraph *myCopyVectorGr=new TGraph(entries,xVector,yVector);
    //myCopyVectorGr->SetLineColor(kYellow+10);
    //myCopyVectorGr->SetLineStyle(2);
    //myCopyVectorGr->SetLineWidth(2);
    //myCopyVectorGr->Draw("LP");
    
    //TMarker *mark=new TMarker();
    //mark->SetMarkerStyle(2);
    //mark->SetMarkerSize(1);
    
    TGraph *myCopyScalarGr=new TGraph(entries,xScalar,yScalar);
    myCopyScalarGr->SetLineColor(kYellow+10);
    //myCopyScalarGr->SetLineStyle(4);
    //myCopyScalarGr->SetLineWidth(2);
    
    //myCopyScalarGr->SetMarkerStyle(20);
    myCopyScalarGr->Draw("LP");
    
    //TMarker *mark2=new TMarker();
    //mark2->SetMarkerStyle(4);

    
    TGraph *myCopyG4Gr=new TGraph(entries,xGeant4,yGeant4);
    myCopyG4Gr->SetLineColor(kRed);
    //myCopyG4Gr->SetLineStyle(5);
    //myCopyG4Gr->SetLineWidth(2);
    myCopyG4Gr->Draw("LP");
    
    TLegend *leg = new TLegend(0.4,0.6,0.9,0.9);
    leg->SetHeader("Legend");
    //leg->AddEntry(my,"Histogram filled with random numbers","f");
    //leg->AddEntry(my,"Function abs(#frac{sin(x)}{x})","l");
    //leg->AddEntry("gr","Graph with error bars","lep");
    TString legName=Form("Pdf for E_in=%s MeV",theApp.Argv(3));
    leg->AddEntry(pdfGraph,legName,"l");
    leg->AddEntry(myCopyScalarGr,"E_out Histogram values scaled - Scalar","l");
    //leg->AddEntry(myCopyVectorGr,"E_out Histogram values scaled - Vector","l");
    leg->AddEntry(myCopyG4Gr,"E_out Histogram values scaled - Geant4","l");
    
    leg->Draw();
    c3->Update();
    std::cout<<"5\n";
    
    ////Calculate chi-square
    std::cout<<"#E_out Histogram entries: "<<myCopyScalarGr->GetN()<<"#Pdf entries: "<<pdfGraph->GetN()<<"\n";
    double chiSquare_vector=0,chiSquare_scalar=0, chiSquare_geant4=0, expX,expY ;
    for (int i=0; i<pdfGraph->GetN(); i++)
    {
        pdfGraph->GetPoint(i, expX, expY);
        //myCopyScalarGr->GetPoint(i, obsX, obsY);
        //std::cout<<"expX: "<<expX<<" expY: "<<expY<<" obsX: "<<obsX<<" obsY: "<<obsY<<std::endl;
        chiSquare_scalar+=(yScalar[i]-expY)*(yScalar[i]-expY)/expY;
        chiSquare_vector+=(yVector[i]-expY)*(yVector[i]-expY)/expY;
        chiSquare_geant4+=(yGeant4[i]-expY)*(yGeant4[i]-expY)/expY;
    }
    std::cout<<"chiSquare scalar: "<<chiSquare_scalar<<" chiSquare-reduced: "<<chiSquare_scalar/(eOutScalarScaled->GetSize()-2)<<std::endl;
    std::cout<<"chiSquare vector: "<<chiSquare_vector<<" chiSquare-reduced: "<<chiSquare_vector/(eOutVectorScaled->GetSize()-2)<<std::endl;
    std::cout<<"chiSquare geant4: "<<chiSquare_geant4<<" chiSquare-reduced: "<<chiSquare_geant4/(eOutGeant4Scaled->GetSize()-2)<<std::endl;
    
    
    //gr[i]->SetTitle(graphName);
    //gr[i]->GetXaxis()->SetTitle(" EnergyOut [MeV]");
    //gr[i]->GetYaxis()->SetTitle("Cross section");
    //myCopyGr->Write();
    c1->SaveAs(histoFileName);
     */
  	std::cout << "\nRun completed.\n";
    //theApp.Run();
	return 0;
}
