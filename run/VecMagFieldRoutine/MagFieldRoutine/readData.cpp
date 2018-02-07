#include <fstream>
#include <iostream>
#include "stdio.h"
#include <sstream>
#include <vector>
using namespace std;

int main(){
	string line;
	int row,col;
	int no_rows=29141;
	int no_cols=6;
	double test_data[no_rows][no_cols];
	string s1,s2,s3,s4,s5,s0;
	//ifstream pFile("cms2015.ascii");
	ifstream pFile("cms2015.txt");
	if (pFile.is_open()){
		row=0;
		// getline() returns the stream. testing the stream with while returns error such as EOF
		while(getline(pFile,line)){
			 // so here we know that the read was a success and that line has valid data
			stringstream ss(line);

			//parsing all the parts. s0's store the string names which are of no use to us. 
			ss>> s0>>test_data[row][0] >>s1 >>test_data[row][1] >>s2 >> test_data[row][2] >>s3 >> test_data[row][3] >>s4 >> test_data[row][4]>> s5 >> test_data[row][5];

			row++;

		}
		pFile.close();
	}
	else
		cout<<"Unable to open file";


	// for(int i=0; i<no_rows/2;i++){
	// 	for(int j=0;j<no_cols;j++){
	// 		cout<<test_data[i][j]<<" ";
	// 	}

	// 	cout<<"\n";
	// }
	
	return 0;
}