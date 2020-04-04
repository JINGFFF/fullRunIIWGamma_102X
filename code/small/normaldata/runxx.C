#include "xx.h"
#include "xx.C"

int main(int argc, char** argv) {

//gROOT->LoadMacro("xx.C");
TString in = argv[1];
//TString dir= "/home/pku/pengj/VBSWGamma/94X2016/data/94x_control_smuH_v1/191217_082115/0000/";
TString dir= "/home/pku/pengj/VBSWGamma/94X2016/data/";  // -loose-iso/";
TString outdir = "/home/pku/pengj/VBSWGamma/94X2016/fake_lepton_rate/data/";
ifstream infile(in);
string buffer; 
TString infilename;

int k=1;

while (k>0){
getline (infile, buffer) ;
infilename = buffer;
if(infilename.Contains(".root")==0) {k=-2; continue;}
TString outname= outdir + "out"+infilename;

cout<<outname<<endl;

TFile *file1 =new TFile(dir+infilename);
TDirectory * dir1 = (TDirectory*)file1->Get("treeDumper");
TTree *tree1 = (TTree*) dir1->Get("PKUCandidates");
xx m1(tree1,outname);
cout<<outname<<endl;
m1.Loop();
m1.endJob();
}
return 0; 



}

