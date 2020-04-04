#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TTree.h"
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
//barrel
TH1D* hist_barrel_pt_25_30 	= new TH1D("hist_barrel_pt_25_30", "hist_barrel_pt_25_30", 128, 0.00, 0.04);
TH1D* hist_barrel_pt_30_40 	= new TH1D("hist_barrel_pt_30_40", "hist_barrel_pt_30_40", 128, 0.00, 0.04);
TH1D* hist_barrel_pt_40_50 	= new TH1D("hist_barrel_pt_40_50", "hist_barrel_pt_40_50", 128, 0.00, 0.04);
TH1D* hist_barrel_pt_50_70 	= new TH1D("hist_barrel_pt_50_70", "hist_barrel_pt_50_70", 128, 0.00, 0.04);
TH1D* hist_barrel_pt_70_100 	= new TH1D("hist_barrel_pt_70_100", "hist_barrel_pt_70_100", 128, 0.00, 0.04);
TH1D* hist_barrel_pt_100_135 	= new TH1D("hist_barrel_pt_100_135", "hist_barrel_pt_100_135", 128, 0.00, 0.04);
TH1D* hist_barrel_pt_135_400 	= new TH1D("hist_barrel_pt_135_400", "hist_barrel_pt_135_400", 128, 0.00, 0.04);
//endcap
TH1D* hist_endcap_pt_25_30      = new TH1D("hist_endcap_pt_25_30", "hist_endcap_pt_25_30", 128, 0.00, 0.04);
TH1D* hist_endcap_pt_30_40      = new TH1D("hist_endcap_pt_30_40", "hist_endcap_pt_30_40", 128, 0.00, 0.04);
TH1D* hist_endcap_pt_40_50      = new TH1D("hist_endcap_pt_40_50", "hist_endcap_pt_40_50", 128, 0.00, 0.04);
TH1D* hist_endcap_pt_50_70      = new TH1D("hist_endcap_pt_50_70", "hist_endcap_pt_50_70", 128, 0.00, 0.04);
TH1D* hist_endcap_pt_70_100     = new TH1D("hist_endcap_pt_70_100", "hist_endcap_pt_70_100", 128, 0.00, 0.04);
TH1D* hist_endcap_pt_100_135    = new TH1D("hist_endcap_pt_100_135", "hist_endcap_pt_100_135", 128, 0.00, 0.04);
TH1D* hist_endcap_pt_135_400    = new TH1D("hist_endcap_pt_135_400", "hist_endcap_pt_135_400", 128, 0.00, 0.04);

void templateMaker(){
    TTreeReader fReader("demo", ipt);
    TTreeReaderValue<Double_t> scalef              = {fReader, "scalef"};
    TTreeReaderValue<Int_t>    nevent              = {fReader, "nevent"};
    TTreeReaderValue<Int_t>    run                 = {fReader, "run"};
    TTreeReaderValue<Int_t>    ls                  = {fReader, "ls"};
    TTreeReaderValue<Int_t>    nVtx                = {fReader, "nVtx"};
    TTreeReaderValue<Double_t> theWeight           = {fReader, "theWeight"};
    TTreeReaderValue<Double_t> npT                 = {fReader, "npT"};
    TTreeReaderValue<Int_t>    lep                 = {fReader, "lep"};
    TTreeReaderValue<Double_t> ptVlepJEC           = {fReader, "ptVlepJEC"};
    TTreeReaderValue<Double_t> yVlepJEC            = {fReader, "yVlepJEC"};
    TTreeReaderValue<Double_t> phiVlepJEC          = {fReader, "phiVlepJEC"};
    TTreeReaderValue<Double_t> mtVlepJECnew        = {fReader, "mtVlepJECnew"};
    TTreeReaderValue<Int_t>    ngoodeles          = {fReader, "ngoodeles"};
    TTreeReaderValue<Int_t>    ngoodmus           = {fReader, "ngoodmus"};
    TTreeReaderValue<Int_t>    nlooseeles          = {fReader, "nlooseeles"};
    TTreeReaderValue<Int_t>    nloosemus           = {fReader, "nloosemus"};
    TTreeReaderArray<Double_t> photon_pt           = {fReader, "photon_pt"};
    TTreeReaderArray<Double_t> photon_eta          = {fReader, "photon_eta"};
    TTreeReaderArray<Double_t> photon_phi          = {fReader, "photon_phi"};
    TTreeReaderArray<Bool_t>   photon_ppsv         = {fReader, "photon_ppsv"};
    TTreeReaderArray<Double_t> photon_hoe          = {fReader, "photon_hoe"};
    TTreeReaderArray<Double_t> photon_sieie        = {fReader, "photon_sieie"};
    TTreeReaderArray<Double_t> photon_sieie2       = {fReader, "photon_sieie2"};
    TTreeReaderArray<Double_t> photon_chiso        = {fReader, "photon_chiso"};
    TTreeReaderArray<Double_t> photon_nhiso        = {fReader, "photon_nhiso"};
    TTreeReaderArray<Double_t> photon_phoiso       = {fReader, "photon_phoiso"};
    TTreeReaderArray<Int_t>    photon_isprompt     = {fReader, "photon_isprompt"};
    TTreeReaderArray<Double_t> photon_drla         = {fReader, "photon_drla"};
    TTreeReaderArray<Double_t> photon_mla          = {fReader, "photon_mla"};
    TTreeReaderArray<Double_t> photonsc_eta        = {fReader, "photonsc_eta"};
    TTreeReaderArray<Double_t> photonsc_phi        = {fReader, "photonsc_phi"};
    TTreeReaderValue<Double_t> ptlep1              = {fReader, "ptlep1"};
    TTreeReaderValue<Double_t> etalep1             = {fReader, "etalep1"};
    TTreeReaderValue<Double_t> philep1             = {fReader, "philep1"};
    TTreeReaderValue<Int_t>    isprompt            = {fReader, "isprompt"};
    TTreeReaderValue<Int_t>    ispromptLep         = {fReader, "ispromptLep"};
    TTreeReaderValue<Double_t> MET_et              = {fReader, "MET_et"};
    TTreeReaderValue<Double_t> MET_phi             = {fReader, "MET_phi"};
    TTreeReaderValue<Int_t>    HLT_Ele1            = {fReader, "HLT_Ele1"};
    TTreeReaderValue<Int_t>    HLT_Ele2            = {fReader, "HLT_Ele2"};
    TTreeReaderValue<Int_t>    HLT_Mu1             = {fReader, "HLT_Mu1"};
    TTreeReaderValue<Int_t>    HLT_Mu2             = {fReader, "HLT_Mu2"};
    TTreeReaderValue<Int_t>    HLT_Mu3             = {fReader, "HLT_Mu3"};
    TTreeReaderValue<Double_t> lumiWeight          = {fReader, "lumiWeight"};
    TTreeReaderValue<Double_t> pileupWeight        = {fReader, "pileupWeight"};
    TTreeReaderValue<Double_t> photonsceta         = {fReader, "photonsceta"};
    TTreeReaderValue<Double_t> photonscphi         = {fReader, "photonscphi"};
    TTreeReaderValue<Double_t> Mla                 = {fReader, "Mla"};
    TTreeReaderValue<Bool_t>   photonhaspixelseed  = {fReader, "photonhaspixelseed"};
    TTreeReaderValue<Bool_t>   photonpasseleveto   = {fReader, "photonpasseleveto"};
    TTreeReaderValue<Double_t> photonet            = {fReader, "photonet"};
    TTreeReaderValue<Double_t> photoneta           = {fReader, "photoneta"};
    TTreeReaderValue<Double_t> photonphi           = {fReader, "photonphi"};
    TTreeReaderValue<Double_t> photone             = {fReader, "photone"};
    TTreeReaderValue<Double_t> drla                = {fReader, "drla"};
    TTreeReaderValue<Double_t> jet1pt              = {fReader, "jet1pt"};
    TTreeReaderValue<Double_t> jet1eta             = {fReader, "jet1eta"};
    TTreeReaderValue<Double_t> jet1phi             = {fReader, "jet1phi"};
    TTreeReaderValue<Double_t> jet1e               = {fReader, "jet1e"};
    TTreeReaderValue<Double_t> jet1icsv            = {fReader, "jet1icsv"};
    TTreeReaderValue<Double_t> jet2pt              = {fReader, "jet2pt"};
    TTreeReaderValue<Double_t> jet2eta             = {fReader, "jet2eta"};
    TTreeReaderValue<Double_t> jet2phi             = {fReader, "jet2phi"};
    TTreeReaderValue<Double_t> jet2e               = {fReader, "jet2e"};
    TTreeReaderValue<Double_t> jet2icsv            = {fReader, "jet2icsv"};
    TTreeReaderValue<Double_t> drj1a               = {fReader, "drj1a"};
    TTreeReaderValue<Double_t> drj2a               = {fReader, "drj2a"};
    TTreeReaderValue<Double_t> drj1l               = {fReader, "drj1l"};
    TTreeReaderValue<Double_t> drj2l               = {fReader, "drj2l"};
    TTreeReaderValue<Double_t> Mjj                 = {fReader, "Mjj"};
    TTreeReaderValue<Double_t> deltaeta            = {fReader, "deltaeta"};
    TTreeReaderValue<Double_t> zepp                = {fReader, "zepp"};
    TTreeReaderValue<Double_t> j1metPhi            = {fReader, "j1metPhi"};
    TTreeReaderValue<Double_t> j2metPhi            = {fReader, "j2metPhi"};
    TTreeReaderValue<Double_t> Dphiwajj            = {fReader, "



}

int main(int argc, char** argv){
    startTime = clock();
    TString type= argv[1];
    if(!(type == "MuonChannel_data" || type == "MuonChannel_fake" type == "MuonChannel_true"|| type == "ElectronChannel_data" || type == "ElectronChannel_fake" type == "ElectronChannel_true")){
        cout<<endl<<"Please input 'MuonChannel_data', 'MuonChannel_fake', 'MuonChannel_true', 'ElectronChannel_data', 'ElectronChannel_fake', or  'ElectronChannel_true' as the type!"<<endl<<endl;
        return 0;
    }

    cout<<"==========================================="<<endl;
    cout<<"   Make "<<type<<" template   "<<endl;
    cout<<"==========================================="<<endl<<endl;

    TString input_file = "input.txt";
    ifstream in(input_file);

    string line;
    int ww=0;
    TString muon, electron, mc;

    while (getline (in, line)){
        TString test = line;
        if(test.Contains("#")) continue;
        if (ww==0) muon 	= line;
        if (ww==1) electron 	= line;
        if (ww==2) mc 		= line;
        ww++;
    }
    //file in
    TFile* file_muon		= TFile::Open(muon);
    TFile* file_electronon	= TFile::Open(electron);
    TFile* file_mc    		= TFile::Open(mc);


}















