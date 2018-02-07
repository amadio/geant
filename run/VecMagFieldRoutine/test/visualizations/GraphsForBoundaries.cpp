#include "iostream"
#include "MagField.h"
#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"
#include <vector>

//For plotting graphs using root 
#include <TGraph.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TMultiGraph.h>

using namespace std;


int main(int argc, char *argv[]){
    MagField m1;
    m1.ReadVectorData("/home/ananya/Work/MagFieldRoutine/cms2015.txt");

    const double kRDiff = 50;
    const double kZDiff = 200;
    const double kRMax= 9000;
    const double kZMax= 16000;

    //For storing values in order to plot graphs later 
    vector<double> allR, allZ;
    vector<double> z1Bz, z2Bz, z3Bz, z4Bz, z5Bz;
    vector<double> z1Br, z2Br, z3Br, z4Br, z5Br;
    vector<double> z1Bphi, z2Bphi, z3Bphi, z4Bphi, z5Bphi;
    vector<double> r1Bz, r2Bz, r3Bz, r4Bz, r5Bz;
    vector<double> r1Br, r2Br, r3Br, r4Br, r5Br;
    vector<double> r1Bphi, r2Bphi, r3Bphi, r4Bphi, r5Bphi;

    //Z boundaries 
    for (double r = 0.; r <= kRMax; r =r+kRDiff){
        
        //(r,0,z) == (r,z)
        // z = -16k, -8k, 0, 8k, 16k
        vecgeom::Vector3D<double> pos1(r,0.,-kZMax), pos2(r,0.,-kZMax*0.5), pos3(r,0.,0.), pos4(r,0.,kZMax*0.5), pos5(r,0.,kZMax);
        vecgeom::Vector3D<double> rzField1, rzField2, rzField3, rzField4, rzField5;
        m1.GetFieldValueTest(pos1, rzField1);
        m1.GetFieldValueTest(pos2, rzField2);
        m1.GetFieldValueTest(pos3, rzField3);
        m1.GetFieldValueTest(pos4, rzField4);
        m1.GetFieldValueTest(pos5, rzField5);
        allR.push_back(pos1.x()); z1Br.push_back(rzField1.x()); z1Bphi.push_back(rzField1.y()); z1Bz.push_back(rzField1.z()); 
        z2Br.push_back(rzField2.x()); z2Bphi.push_back(rzField2.y()); z2Bz.push_back(rzField2.z()); 
        z3Br.push_back(rzField3.x()); z3Bphi.push_back(rzField3.y()); z3Bz.push_back(rzField3.z()); 
        z4Br.push_back(rzField4.x()); z4Bphi.push_back(rzField4.y()); z4Bz.push_back(rzField4.z()); 
        z5Br.push_back(rzField5.x()); z5Bphi.push_back(rzField5.y()); z5Bz.push_back(rzField5.z()); 
    }
    

    //r-boundaries 
    for (double z = -kZMax; z <= kZMax; z = z+kZDiff){

        //r = 0, 2250,4500,6750,9000
        vecgeom::Vector3D<double> pos1(0.,0.,z), pos2(kRMax*0.25,0.,z), pos3(kRMax*0.5,0.,z), pos4(kRMax*0.75,0.,z), pos5(kRMax,0.,z);
        vecgeom::Vector3D<double> rzField1, rzField2, rzField3, rzField4, rzField5;
        m1.GetFieldValueTest(pos1, rzField1);
        m1.GetFieldValueTest(pos2, rzField2);
        m1.GetFieldValueTest(pos3, rzField3);
        m1.GetFieldValueTest(pos4, rzField4);
        m1.GetFieldValueTest(pos5, rzField5);
        allZ.push_back(pos1.z()); r1Br.push_back(rzField1.x()); r1Bphi.push_back(rzField1.y()); r1Bz.push_back(rzField1.z()); 
        r2Br.push_back(rzField2.x()); r2Bphi.push_back(rzField2.y()); r2Bz.push_back(rzField2.z()); 
        r3Br.push_back(rzField3.x()); r3Bphi.push_back(rzField3.y()); r3Bz.push_back(rzField3.z()); 
        r4Br.push_back(rzField4.x()); r4Bphi.push_back(rzField4.y()); r4Bz.push_back(rzField4.z()); 
        r5Br.push_back(rzField5.x()); r5Bphi.push_back(rzField5.y()); r5Bz.push_back(rzField5.z()); 

    }

    

    int flag = 6;//to decide which graphs to display
    TApplication *app = new TApplication("App", &argc, argv);
    int style1, style2, style3, style4, style5;
    style1 = style2 = style3 = style4 = style5 = 7;
    int color5= 5; //yellow in place of previous black
    int color4= 28; //trying dark brown

    
    TCanvas *c1 = new TCanvas("c1","Magnetic Field on Boundaries ",200,10,700,500);
    TMultiGraph *mg = new TMultiGraph("mg", "Magnetic Field on Boundaries");
    //c1->Divide(1,2);

    if(flag == 1){
        TGraph *gr1 = new TGraph( allR.size(), &(allR[0]), &(z1Bz[0]) );
        //gr1->SetLineColor(2);
        // gr->SetLineWidth(4);
        gr1->SetMarkerColor(4); //Blue 
        gr1->SetMarkerStyle(style1);
        //gr1->SetTitle("z = -16000");
        gr1->GetXaxis()->SetTitle("radius (mm)");
        gr1->GetYaxis()->SetTitle("Bz (tesla)");
        gr1->GetYaxis()->SetRangeUser(-2.0, 4.2);

        TGraph *gr2 = new TGraph( allR.size(), &(allR[0]), &(z2Bz[0]) );
        gr2->SetMarkerColor(2); //red 
        gr2->SetMarkerStyle(style2);
        
        TGraph *gr3 = new TGraph( allR.size(), &(allR[0]), &(z3Bz[0]) );
        gr3->SetMarkerColor(3); //green
        gr3->SetMarkerStyle(style3);
        
        TGraph *gr4 = new TGraph( allR.size(), &(allR[0]), &(z4Bz[0]) );
        gr4->SetMarkerColor(color4); //some shade of red
        gr4->SetMarkerStyle(style4);
        
        TGraph *gr5 = new TGraph( allR.size(), &(allR[0]), &(z5Bz[0]) );
        gr5->SetMarkerColor(color5); //black
        gr5->SetMarkerStyle(style5);

        // gr1->Draw("ALP");
        // gr2->Draw("LP");
        // gr3->Draw("LP");
        // gr4->Draw("LP");
        // gr5->Draw("LP");

        gr1->SetTitle("z = -16000");
        gr2->SetTitle("z= - 8000");
        gr3->SetTitle("z = 0");
        gr4->SetTitle("z = 8000");
        gr5->SetTitle("z=16000");


        mg->Add(gr2);
        mg->Add(gr3);
        mg->Add(gr4);
        mg->Add(gr5);
        gr1->Draw("ALP");
        mg->Draw("LP");
        c1->BuildLegend();

    }
    else if(flag == 2){
        TGraph *gr1 = new TGraph( allR.size(), &(allR[0]), &(z1Br[0]) );
        gr1->SetMarkerColor(4); //Blue 
        gr1->SetMarkerStyle(style1);
        //gr1->SetTitle("Z boundaries");
        gr1->GetXaxis()->SetTitle("radius (mm)");
        gr1->GetYaxis()->SetTitle("Br (tesla)");
        gr1->GetYaxis()->SetRangeUser(-0.5,0.5);

        TGraph *gr2 = new TGraph( allR.size(), &(allR[0]), &(z2Br[0]) );
        gr2->SetMarkerColor(2); //red 
        gr2->SetMarkerStyle(style2);

        TGraph *gr3 = new TGraph( allR.size(), &(allR[0]), &(z3Br[0]) );
        gr3->SetMarkerColor(3); //green
        gr3->SetMarkerStyle(style3);

        TGraph *gr4 = new TGraph( allR.size(), &(allR[0]), &(z4Br[0]) );
        gr4->SetMarkerColor(color4); //some shade of red
        gr4->SetMarkerStyle(style4);

        TGraph *gr5 = new TGraph( allR.size(), &(allR[0]), &(z5Br[0]) );
        gr5->SetMarkerColor(color5); //black
        gr5->SetMarkerStyle(style5);

        // gr1->Draw("ALP");
        // gr2->Draw("LP");
        // gr3->Draw("LP");
        // gr4->Draw("LP");
        // gr5->Draw("LP");

        gr1->SetTitle("z = -16000");
        gr2->SetTitle("z= - 8000");
        gr3->SetTitle("z = 0");
        gr4->SetTitle("z = 8000");
        gr5->SetTitle("z=16000");

        mg->Add(gr2);
        mg->Add(gr3);
        mg->Add(gr4);
        mg->Add(gr5);
        gr1->Draw("ALP");
        mg->Draw("LP");
        c1->BuildLegend();
    }
    else if(flag == 3){
        TGraph *gr1 = new TGraph( allR.size(), &(allR[0]), &(z1Bphi[0]) );
        gr1->SetMarkerColor(4); //Blue 
        gr1->SetMarkerStyle(style1);
        //gr1->SetTitle("Z boundaries");
        gr1->GetXaxis()->SetTitle("radius (mm)");
        gr1->GetYaxis()->SetTitle("Bphi (tesla)");
        gr1->GetYaxis()->SetRangeUser(-0.01,0.07);

        TGraph *gr2 = new TGraph( allR.size(), &(allR[0]), &(z2Bphi[0]) );
        gr2->SetMarkerColor(2); //red 
        gr2->SetMarkerStyle(style2);

        TGraph *gr3 = new TGraph( allR.size(), &(allR[0]), &(z3Bphi[0]) );
        gr3->SetMarkerColor(3); //green
        gr3->SetMarkerStyle(style3);

        TGraph *gr4 = new TGraph( allR.size(), &(allR[0]), &(z4Bphi[0]) );
        gr4->SetMarkerColor(color4); //some shade of red
        gr4->SetMarkerStyle(style4);

        TGraph *gr5 = new TGraph( allR.size(), &(allR[0]), &(z5Bphi[0]) );
        gr5->SetMarkerColor(color5); //black
        gr5->SetMarkerStyle(style5);

        // gr1->Draw("ALP");
        // gr2->Draw("LP");
        // gr3->Draw("LP");
        // gr4->Draw("LP");
        // gr5->Draw("LP");

        gr1->SetTitle("z = -16000");
        gr2->SetTitle("z= - 8000");
        gr3->SetTitle("z = 0");
        gr4->SetTitle("z = 8000");
        gr5->SetTitle("z=16000");

        mg->Add(gr2);
        mg->Add(gr3);
        mg->Add(gr4);
        mg->Add(gr5);
        gr1->Draw("ALP");
        mg->Draw("LP");
        c1->BuildLegend();
    }

    else if(flag==4){
        TGraph *gr1 = new TGraph( allZ.size(), &(allZ[0]), &(r1Bphi[0]) );
        gr1->SetMarkerColor(4); //Blue 
        gr1->SetMarkerStyle(style1);
        gr1->SetTitle("R boundaries");
        gr1->GetXaxis()->SetTitle("z (mm)");
        gr1->GetYaxis()->SetTitle("Bphi (tesla)");
        gr1->GetYaxis()->SetRangeUser(-0.8,0.8);

        TGraph *gr2 = new TGraph( allZ.size(), &(allZ[0]), &(r2Bphi[0]) );
        gr2->SetMarkerColor(2); //red 
        gr2->SetMarkerStyle(style2);

        TGraph *gr3 = new TGraph( allZ.size(), &(allZ[0]), &(r3Bphi[0]) );
        gr3->SetMarkerColor(3); //green
        gr3->SetMarkerStyle(style3);

        TGraph *gr4 = new TGraph( allZ.size(), &(allZ[0]), &(r4Bphi[0]) );
        gr4->SetMarkerColor(color4); //some shade of red
        gr4->SetMarkerStyle(style4);

        TGraph *gr5 = new TGraph( allZ.size(), &(allZ[0]), &(r5Bphi[0]) );
        gr5->SetMarkerColor(color5); //black
        gr5->SetMarkerStyle(style5);

        // gr1->Draw("ALP");
        // gr2->Draw("LP");
        // gr3->Draw("LP");
        // gr4->Draw("LP");
        // gr5->Draw("LP");

        gr1->SetTitle("r = 0");
        gr2->SetTitle("r = 2250");
        gr3->SetTitle("r = 4500");
        gr4->SetTitle("r = 6750");
        gr5->SetTitle("r = 9000");  

        mg->Add(gr2);
        mg->Add(gr3);
        mg->Add(gr4);
        mg->Add(gr5);
        gr1->Draw("ALP");
        mg->Draw("LP");
        c1->BuildLegend();  
    }

    else if(flag==5){
        TGraph *gr1 = new TGraph( allZ.size(), &(allZ[0]), &(r1Br[0]) );
        gr1->SetMarkerColor(4); //Blue 
        gr1->SetMarkerStyle(style1);
        gr1->SetTitle("R boundaries");
        gr1->GetXaxis()->SetTitle("z (mm)");
        gr1->GetYaxis()->SetTitle("Br (tesla)");
        gr1->GetYaxis()->SetRangeUser(-2.5,2.5);

        TGraph *gr2 = new TGraph( allZ.size(), &(allZ[0]), &(r2Br[0]) );
        gr2->SetMarkerColor(2); //red 
        gr2->SetMarkerStyle(style2);

        TGraph *gr3 = new TGraph( allZ.size(), &(allZ[0]), &(r3Br[0]) );
        gr3->SetMarkerColor(3); //green
        gr3->SetMarkerStyle(style3);

        TGraph *gr4 = new TGraph( allZ.size(), &(allZ[0]), &(r4Br[0]) );
        gr4->SetMarkerColor(color4); //some shade of red
        gr4->SetMarkerStyle(style4);

        TGraph *gr5 = new TGraph( allZ.size(), &(allZ[0]), &(r5Br[0]) );
        gr5->SetMarkerColor(color5); //black
        gr5->SetMarkerStyle(style5);

        // gr1->Draw("ALP");
        // gr2->Draw("LP");
        // gr3->Draw("LP");
        // gr4->Draw("LP");
        // gr5->Draw("LP");     

        gr1->SetTitle("r = 0");
        gr2->SetTitle("r = 2250");
        gr3->SetTitle("r = 4500");
        gr4->SetTitle("r = 6750");
        gr5->SetTitle("r = 9000");  

        mg->Add(gr2);
        mg->Add(gr3);
        mg->Add(gr4);
        mg->Add(gr5);
        gr1->Draw("ALP");
        mg->Draw("LP");
        c1->BuildLegend();
    }

    else if(flag==6){
        TGraph *gr1 = new TGraph( allZ.size(), &(allZ[0]), &(r1Bz[0]) );
        gr1->SetMarkerColor(4); //Blue 
        gr1->SetMarkerStyle(style1);
        gr1->SetTitle("R boundaries");
        gr1->GetXaxis()->SetTitle("z (mm)");
        gr1->GetYaxis()->SetTitle("Bz (tesla)");
        gr1->GetYaxis()->SetRangeUser(-2,4);

        TGraph *gr2 = new TGraph( allZ.size(), &(allZ[0]), &(r2Bz[0]) );
        gr2->SetMarkerColor(2); //red 
        gr2->SetMarkerStyle(style2);

        TGraph *gr3 = new TGraph( allZ.size(), &(allZ[0]), &(r3Bz[0]) );
        gr3->SetMarkerColor(3); //green
        gr3->SetMarkerStyle(style3);

        TGraph *gr4 = new TGraph( allZ.size(), &(allZ[0]), &(r4Bz[0]) );
        gr4->SetMarkerColor(color4); //some shade of red
        gr4->SetMarkerStyle(style4);

        TGraph *gr5 = new TGraph( allZ.size(), &(allZ[0]), &(r5Bz[0]) );
        gr5->SetMarkerColor(color5); //black
        gr5->SetMarkerStyle(style5);

        // gr1->Draw("ALP");
        // gr2->Draw("LP");
        // gr3->Draw("LP");
        // gr4->Draw("LP");
        // gr5->Draw("LP");     
        gr1->SetTitle("r = 0");
        gr2->SetTitle("r = 2250");
        gr3->SetTitle("r = 4500");
        gr4->SetTitle("r = 6750");
        gr5->SetTitle("r = 9000");  

        mg->Add(gr2);
        mg->Add(gr3);
        mg->Add(gr4);
        mg->Add(gr5);
        gr1->Draw("ALP");
        mg->Draw("LP");
        c1->BuildLegend();
    }

    app->Run();

}