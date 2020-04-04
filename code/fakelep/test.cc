#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
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

Double_t LUMI=-35.857;
clock_t startTime,endTime;
TH2D* t_l;
TH2D* c_l;
TH2D* t_wz;
TH2D* c_wz;

void count(TString type, TFile* ipt,TString filename, TString str) {
    TTreeReader fReader("demo", ipt);
    TTreeReaderValue<Double_t> scalef              = {fReader, "scalef"};
    TTreeReaderValue<Int_t>    nevent              = {fReader, "nevent"};
    TTreeReaderValue<Int_t>    run                 = {fReader, "run"};
    TTreeReaderValue<Int_t>    ls                  = {fReader, "ls"};
    TTreeReaderValue<Double_t> theWeight           = {fReader, "theWeight"};        
    TTreeReaderValue<Int_t>    lep                 = {fReader, "lep"};
    TTreeReaderValue<Int_t>    jet1hf_orig         = {fReader, "jet1hf_orig"};
    TTreeReaderValue<Int_t>    jet1pf_orig         = {fReader, "jet1pf_orig"};
    TTreeReaderValue<Int_t>    jet2hf_orig         = {fReader, "jet2hf_orig"};
    TTreeReaderValue<Int_t>    jet2pf_orig         = {fReader, "jet2pf_orig"};
    TTreeReaderValue<Double_t> jet1pt_orig         = {fReader, "jet1pt_orig"};
    TTreeReaderValue<Double_t> jet1eta_orig        = {fReader, "jet1eta_orig"};
    TTreeReaderValue<Double_t> jet1phi_orig        = {fReader, "jet1phi_orig"};
    TTreeReaderValue<Double_t> jet1e_orig          = {fReader, "jet1e_orig"};
    TTreeReaderValue<Double_t> jet2pt_orig         = {fReader, "jet2pt_orig"};
    TTreeReaderValue<Double_t> jet2eta_orig        = {fReader, "jet2eta_orig"};
    TTreeReaderValue<Double_t> jet2phi_orig        = {fReader, "jet2phi_orig"};
    TTreeReaderValue<Double_t> jet2e_orig          = {fReader, "jet2e_orig"};
    TTreeReaderValue<Double_t> jet1csv_orig        = {fReader, "jet1csv_orig"};
    TTreeReaderValue<Double_t> jet2csv_orig        = {fReader, "jet2csv_orig"};
    TTreeReaderValue<Double_t> jet1icsv_orig       = {fReader, "jet1icsv_orig"};
    TTreeReaderValue<Double_t> jet2icsv_orig       = {fReader, "jet2icsv_orig"};
    TTreeReaderValue<Double_t> drj1a_orig          = {fReader, "drj1a_orig"};
    TTreeReaderValue<Double_t> drj2a_orig          = {fReader, "drj2a_orig"};
    TTreeReaderValue<Double_t> drj1l_orig          = {fReader, "drj1l_orig"};
    TTreeReaderValue<Double_t> drj2l_orig          = {fReader, "drj2l_orig"};
    TTreeReaderValue<Int_t>    nlooseeles          = {fReader, "nlooseeles"};
    TTreeReaderValue<Int_t>    nloosemus           = {fReader, "nloosemus"};
    TTreeReaderValue<Int_t>    ngoodeles          = {fReader, "ngoodeles"};
    TTreeReaderValue<Int_t>    ngoodmus           = {fReader, "ngoodmus"};
    TTreeReaderValue<Double_t> ptlep1              = {fReader, "ptlep1"};
    TTreeReaderValue<Double_t> etalep1             = {fReader, "etalep1"};
    TTreeReaderValue<Int_t>    HLT_Ele1            = {fReader, "HLT_Ele1"};
    TTreeReaderValue<Int_t>    HLT_Ele2            = {fReader, "HLT_Ele2"};
    TTreeReaderValue<Int_t>    HLT_Mu1             = {fReader, "HLT_Mu1"};
    TTreeReaderValue<Int_t>    HLT_Mu2             = {fReader, "HLT_Mu2"};
    TTreeReaderValue<Int_t>    HLT_Mu3             = {fReader, "HLT_Mu3"};
    TTreeReaderValue<Double_t> lumiWeight          = {fReader, "lumiWeight"};
    TTreeReaderValue<Double_t> pileupWeight        = {fReader, "pileupWeight"};
    TTreeReaderValue<Double_t> MET_et              = {fReader, "MET_et"};
    TTreeReaderValue<Double_t> mtVlepJECnew        = {fReader, "mtVlepJECnew"};

    Long64_t                   maxEntries          = fReader.GetEntries(false);
    Bool_t muon_cut, electron_cut, cut;
    int i = 0;
    cout<<maxEntries<<"   "<<filename<<endl;

    while (fReader.Next()) {
	if (i % 10000 == 0) {
            std::cout << "Analyzing event " << i <<"  "<<setw(6)<<i*100./maxEntries<<"%"<<"\r" << std::flush;
        }
	i++;
	//if (i>10000) break;	
	//cut definition
        muon_cut = *HLT_Mu2 == 1 &&*lep == 13 &&*ngoodmus == 1 &&*ngoodeles == 0 &&*mtVlepJECnew < 20 &&*MET_et < 30 &&((*jet1pt_orig > 20 && *drj1l_orig > 0.3) || (*jet2pt_orig > 20 && *drj2l_orig > 0.3));
        electron_cut = *HLT_Ele2 == 1 &&*lep == 11 &&*ngoodmus == 0 &&*ngoodeles == 1 &&*mtVlepJECnew < 20 &&*MET_et < 30 &&((*jet1pt_orig > 30 && *drj1l_orig > 0.3) || (*jet2pt_orig > 30 && *drj2l_orig > 0.3));
      	if(type == "muon") cut = muon_cut;
	if(type == "electron") cut = electron_cut;

	//fill histgram
	double ptl1 = *ptlep1;
	if(ptl1 >=50) ptl1 = 45.;
	if(cut){
	    //cout<<"pass"<<endl;
	    if(str == "tight lepton"){
		//cout<<"fill"<<endl;
		t_l->Fill(fabs(*etalep1), ptl1);
	    }

	    else if(str == "loose lepton"){
		c_l->Fill(fabs(*etalep1), ptl1);
	    }

            else if(str == "tight wz"){
                t_wz->Fill(fabs(*etalep1), ptl1, (*scalef)*(*lumiWeight) * (*pileupWeight));
            }

            else if(str == "loose wz"){
                c_wz->Fill(fabs(*etalep1), ptl1, (*scalef)* (*lumiWeight) * (*pileupWeight));
            }
        }  
    }
    cout<<endl<<endl;
}

int main(int argc, char** argv){
    startTime = clock();
    TString type= argv[1];
    if(!(type == "muon" || type == "electron")){
	cout<<endl<<"Please input 'muon' or 'electron' as the type!"<<endl<<endl;
	return 0;
    }

    cout<<"==========================================="<<endl;
    cout<<"   Calculate "<<type<<" fake rate   "<<endl;
    cout<<"==========================================="<<endl<<endl;

    TString input_file_tight_data, input_file_loose_data, input_file_tight_mc, input_file_loose_mc;
    if(type == "muon") {
	input_file_tight_data = "tight_muon.txt";
	input_file_loose_data = "loose_muon.txt";
    }
    if(type == "electron"){ 
	input_file_tight_data = "tight_electron.txt";
        input_file_loose_data = "loose_electron.txt";
    }
    input_file_tight_mc = "tight_wz.txt"; 
    input_file_loose_mc = "loose_wz.txt";
   
    ifstream in_t_data(input_file_tight_data);
    ifstream in_l_data(input_file_loose_data);
    ifstream in_t_mc(input_file_tight_mc);
    ifstream in_l_mc(input_file_loose_mc);


    //file out
    TFile* output = new TFile("fake_" + type + ".root", "RECREATE");

    //4 2D histgrams
    Int_t n_pt_bin, n_eta_bin;
    if(type == "electron"){
	n_pt_bin   = 2;
	n_eta_bin  = 5;
	Double_t pt_bin[3]  = {30, 40, 50};
	Double_t eta_bin[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};

    	t_l = new TH2D("count_tight_" + type, "count_tight_" + type + ";|#eta|;p_{T}(GeV)",n_eta_bin, eta_bin, n_pt_bin, pt_bin);
    	c_l = new TH2D("count_control_" + type, "count_control_" + type + ";|#eta|;p_{T}(GeV)",n_eta_bin, eta_bin, n_pt_bin, pt_bin);
    	t_wz = new TH2D("count_tight_wz", "count_tight_wz ;|#eta|;p_{T}(GeV)",n_eta_bin, eta_bin, n_pt_bin, pt_bin);
    	c_wz = new TH2D("count_control_wz", "count_control_wz ;|#eta|;p_{T}(GeV)",n_eta_bin, eta_bin, n_pt_bin, pt_bin);
    }

    if(type == "muon"){
	n_pt_bin   = 3;
	n_eta_bin  = 5;
	Double_t pt_bin[4]  = {25, 30, 40, 50};
	Double_t eta_bin[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};

    	t_l = new TH2D("count_tight_" + type, "count_tight_" + type + ";|#eta|;p_{T}(GeV)",n_eta_bin, eta_bin, n_pt_bin, pt_bin);
    	c_l = new TH2D("count_control_" + type, "count_control_" + type + ";|#eta|;p_{T}(GeV)",n_eta_bin, eta_bin, n_pt_bin, pt_bin);
    	t_wz = new TH2D("count_tight_wz", "count_tight_wz ;|#eta|;p_{T}(GeV)",n_eta_bin, eta_bin, n_pt_bin, pt_bin);
    	c_wz = new TH2D("count_control_wz_", "count_control_wz ;|#eta|;p_{T}(GeV)",n_eta_bin, eta_bin, n_pt_bin, pt_bin);
    } 


    string line1, line2, line3, line4;
    int ww=0;
    TString tight_lepton, control_lepton, tight_wz, control_wz;

    //tight lepton
    cout<<type<<" : tight lepton "<<endl;
    while (getline (in_t_data, line1)){
        //cout << line1 << endl;
        TString test1 = line1;
        if(test1.Contains("#")) continue;
        tight_lepton = line1;
        TFile* file_tight_lepton = TFile::Open(tight_lepton);
	count(type, file_tight_lepton, test1, "tight lepton");

       
    }

    //loose lepton
    cout<<type<<" : loose lepton "<<endl;
    while (getline (in_l_data, line2)){
        //cout << line2 << endl;
        TString test2 = line2;
        if(test2.Contains("#")) continue;
        control_lepton = line2;
        TFile* file_control_lepton = TFile::Open(control_lepton);
        count(type, file_control_lepton, test2, "loose lepton");

    }

    //tight wz
    cout<<type<<" : tight wz "<<endl;
    while (getline (in_t_mc, line3)){
        //cout << line3 << endl;
        TString test3 = line3;
        if(test3.Contains("#")) continue;
        tight_wz = line3;
        TFile* file_tight_wz = TFile::Open(tight_wz);
        count(type, file_tight_wz, test3, "tight wz");

    }

    //loose wz
    cout<<type<<" : loose wz "<<endl;
    while (getline (in_l_mc, line4)){
        //cout << line4 << endl;
        TString test4 = line4;
        if(test4.Contains("#")) continue;
        control_wz = line4;
        TFile* file_control_wz = TFile::Open(control_wz);
        count(type, file_control_wz, test4, "loose wz");

    }

    //TFile* file_tight_lepton = TFile::Open(tight_lepton);
    //TFile* file_control_lepton = TFile::Open(control_lepton);
    //TFile* file_tight_wz = TFile::Open(tight_wz);
    //TFile* file_control_wz = TFile::Open(control_wz);

    //count(type, file_tight_lepton, "tight lepton");
    //count(type, file_control_lepton, "loose lepton");
    //count(type, file_tight_wz, "tight wz");
    //count(type, file_control_wz, "loose wz");

    //calculate the fake rate
    TH2D* fakerate;
    TH2D* fenzi;
    TH2D* fenmu;
    TH2D* weight; 

    fenzi = (TH2D*)t_l->Clone();
    fenzi->Add(t_wz,LUMI);
    fenzi->SetName("numerator");
    fenzi->SetTitle("numerator");
     
    fenmu = (TH2D*)t_l->Clone();
    fenmu->Add(c_l);
    fenmu->Add(t_wz,LUMI);
    fenmu->Add(c_wz,LUMI);
    fenmu->SetName("denominator");
    fenmu->SetTitle("denominator");

    fakerate = (TH2D*)fenzi->Clone();
    fakerate->Divide(fenmu);
    fakerate->SetName("fakerate");
    fakerate->GetXaxis()->SetTitle("|#eta|");
    fakerate->GetYaxis()->SetTitle("p_{T}(GeV)");
    fakerate->SetTitle("fake_rate");

    weight = (TH2D*)fakerate->Clone();
    weight->SetName("weight");
    weight->GetXaxis()->SetTitle("|#eta|");
    weight->GetYaxis()->SetTitle("p_{T}(GeV)");
    weight->SetTitle("fake_lepton_weight");

    for (int i = 0; i < n_pt_bin; i++) {
        for (int j = 0; j < n_eta_bin; j++) {
            Double_t tmp     = fakerate->GetBinContent(j + 1, n_pt_bin - i);
            Double_t wgt_tmp = tmp / (1 - tmp);
            weight->SetBinContent(j + 1, n_pt_bin - i, wgt_tmp);
            cout<<tmp<<"   ";
        }
        cout<<endl;
    }
    cout<<endl;

    output->cd();
    t_l->Write();
    c_l->Write();
    t_wz->Write();
    c_wz->Write();
    fenzi->Write();
    fenmu->Write();
    fakerate->Write();
    weight->Write();    
    output->Close();

    endTime = clock();
    cout<<"The run time is: "<<(double)(endTime - startTime)/CLOCKS_PER_SEC<<"s"<<endl<<endl;
    
}
