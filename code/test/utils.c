#include <stdio.h>  
#include <stdlib.h>  
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TDirectoryFile.h"
using namespace std;

extern "C"{
const char * string_conventor(char* s, int len){
	string sss(s,len);
	cout<<sss<<endl;
	return sss.c_str();
}
}
