#include <stdio.h>  
#include <stdlib.h>  
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TDirectoryFile.h"
using namespace std;

extern "C"{
void reasdd(string sss){
	//string sss(s,len);
	TString infile = sss;
	cout<<"hhh"<<endl;
    TFile *f1 = new TFile(infile);
    TDirectoryFile *td1 = (TDirectoryFile*)f1->Get("demo");
    TTree *t1 = (TTree*)td1->Get("ntree");

    double ht, theWeight, r = 2.5;
    t1->SetBranchAddress("ht",&ht);
    Long64_t nentries1 = t1->GetEntries();

    for(int i1=0; i1<nentries1;i1++){
        t1->GetEntry(i1);
		//printf("ht= %f\n",ht);
        cout<<ht<<endl;
        //h1->Fill(ht_600to800, 1.65/nentries1);
        
    }
    //return r;
}
}
