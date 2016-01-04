#include <fstream>
#include "iostream"
#include "stdio.h"
#include <sstream>
#include <vector>
using namespace std;

int main(){

	cout<<"Test1"<<endl;

	string line;
	int row,col;
	int no_rows=29141;
	int no_cols=6;
	vector<double> radius,phi,z,Bz,Br,Bphi;
	string s1,s2,s3,s4,s5,s0;
	double d1,d2,d3,d4,d5,d0;
	//ifstream pFile("cms2015.ascii");
	ifstream pFile("cms2015.txt");
	if (pFile.is_open()){
		row=0;
		// getline() returns the stream. testing the stream with while returns error such as EOF
		while(getline(pFile,line)){
			 // so here we know that the read was a success and that line has valid data
			stringstream ss(line);

			//parsing all the parts. s0's store the string names which are of no use to us. 
			ss>> s0>> d1>> s1>> d0>> s2>> d2>> s3>> d3>> s4>> d4>> s5>> d5;
			radius.push_back(d1);
			phi.push_back(d0);
			z.push_back(d2);
			Bz.push_back(d3);
			Br.push_back(d4);
			Bphi.push_back(d5);

			row++;

		}
		pFile.close();
	}
	else
		cout<<"Unable to open file";


	// for (int i = 0; i < radius.size(); ++i)
	// {
	// 	/* code */
	// 	cout<<radius[i]<<" "<<z[i]<<" "<<Bz[i]<< " "<<Bphi[i]<<" "<<Br[i]<<endl;
	// }
	
	return 0;
}