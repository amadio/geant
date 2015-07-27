#include <stdlib.h>
#include "TNudyElementRN.h"

ClassImp(TNudyElementRN)

double TNudyElementRN::fCCodeRange[26]={2.0e+32,3.0e+7,1.1e+6,1.4e+5,3.4e+4,
				    8.6e+3, 3e+3, 1.3e+3, 6e+2,2.9e+2,
				    1.6e+2,83.4,43,23.5,12,6.2,3.5,
				    1.8,0.9,0.5,2.3e-1,1.0e-1,4.6e-2,
					  1.4e-2,8.2e-4,0.0};
int TNudyElementRN::fCCodeColor[26][3]={{ 50,101,200},{ 50,122,200},{ 50,144,200},{ 50,165,200},{ 50,187,200},
					  { 50,208,200},{ 50,200,200},{ 61,200,179},{ 72,200,158},{ 84,200,136},
					  { 95,200,115},{107,200,93 },{118,200,72 },{220,220,1  },{220,206,15 },
					  {220,178,43 },{220,163,58 },{220,149,72 },{220,135,86 },{220,180,120},
					  {220,166,134},{220,152,148},{220,138,162},{220,123,177},{220,109,191},{220, 95,205}};

TNudyElementRN::TNudyElementRN()
{
  fEle = 0;
  fSize = 10;
  fScale = 1;
  fPadding = 1;
  fCoSize = 26;
}
TNudyElementRN::TNudyElementRN(TGeoElementRN* elem, Float_t x, Float_t y)
{
  fEle = elem;
  fSize = 10;
  fScale = 1;
  fPadding = 1;
  fCoSize = 26;
  fBox = NULL;
  fX = x;
  fY = y;
}
Color_t TNudyElementRN::GetColor(double halfLife)
{
  int index;

  for(index=0; fCCodeRange[index] > halfLife && index < fCoSize; index++);
  if(halfLife <= 0)
    return TColor::GetColor(0,0,0);
  return TColor::GetColor(fCCodeColor[index][0],fCCodeColor[index][1],fCCodeColor[index][2]);
}
void TNudyElementRN::Draw(Option_t* option){
  fBox = new TBox(fX+fPadding,fY+fPadding,fX+fSize-fPadding,fY+fSize-fPadding);
  fBox->SetFillColor(GetColor(fEle->HalfLife()));
  fBox->SetLineColor(1);
  fBox->SetLineWidth(2);
  /*  char* name;
  strcpy(name,fEle->GetName());
  cout<<name<<endl;
  int a,z;
  sscanf(name,"%d-%s-%d",&a,name,&z);
  TString nname = *name;*/
  TString name = fEle->GetName();
  //name = name.SubString("-[a-z-A-Z]+-");
  TString toolTip = name+"\nA : "+ TString(Form("%d",fEle->AtomicNo())) + "\n"
    +"ENDF Code : " + TString(Form("%d",fEle->ENDFCode()))+"\n"
    +"Z : "+TString(Form("%d",fEle->MassNo()))+"\n"
    +"HalfLife : "+TString(Form("%e",fEle->HalfLife()));
  fBox->Draw();
  fBox->SetToolTipText(toolTip);
}
void TNudyElementRN::Move(Float_t x, Float_t y)
{
  fX = x;
  fY = y;
  fBox->SetX1((fX+fPadding)*fScale);
  fBox->SetY1((fY+fPadding)*fScale);
  fBox->SetX2((fX+fSize-fPadding)*fScale);
  fBox->SetY2((fY+fSize-fPadding)*fScale);
}
 
 
