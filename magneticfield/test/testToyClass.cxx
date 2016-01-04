#include <iostream>
#include "ToyClass.h"
#include <stdlib.h>
#include "ToyClass1.h"
#include "ToyClassTemplate.h"

using namespace std;

int main(int argc, char *argv[]){

  // int flag = atoi(argv[1]);

  int n = 16;
  double yIn[n], yOut[n];
  srand (time(NULL));
  // srand(flag);
  // srand(9); //passes without sigsev

  for (int i = 0; i < n; ++i)
  {
    yIn[i] = (float) rand()/(RAND_MAX) ;
  }

  
/*  ToyClass* objToyClass = new ToyClass(n);
  int temp = objToyClass->ToyMethod(yIn, yOut);
  cout<<temp<<endl;
  delete objToyClass; */

  ToyClass1* objToyClass1 = new ToyClass1(n);
  objToyClass1->ToyMethod(yIn, yOut);
  delete objToyClass1; 


/*  ToyClassTemplate* objToyClass1 = new ToyClassTemplate(n, 1 ); //4);
  objToyClass1->ToyMethod(yIn, yOut);
  delete objToyClass1; 
*/

/*  if(flag)
  {
    ToyClass objToyClass(n);
    int temp = objToyClass.ToyMethod(yIn, yOut);
    cout<<temp<<endl;
  }
  else
  {
    ToyClass* objToyClass = new ToyClass(n);
    int temp = objToyClass->ToyMethod(yIn, yOut);
    cout<<temp<<endl;
    delete objToyClass; 
  }*/

  // cout<<"\n----Output obtained is: "<<yOut<<endl;


  return 0;
}