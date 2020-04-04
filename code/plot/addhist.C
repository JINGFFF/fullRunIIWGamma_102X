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
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>

using namespace std;

void addhist(int addZA, int addVV, int addTTA, int addSTop, int addWGJJ, int addWGJets, int addDATA, int addFakePhoton, int addFakeLepton, int addDoubleFake){

    TString xtitle;
    TString string[36] = {"barrel_nVtx","barrel_ptlep1","barrel_etalep1","barrel_mtVlepJECnew","barrel_ptVlepJEC","barrel_photonet","barrel_photoneta","barrel_photonsceta","barrel_jet1pt","barrel_jet1eta","barrel_jet2pt","barrel_jet2eta","barrel_Mjj","barrel_zepp","barrel_deltaeta","barrel_MET_et","barrel_Dphiwajj","barrel_Mla","endcap_nVtx","endcap_ptlep1","endcap_etalep1","endcap_mtVlepJECnew","endcap_ptVlepJEC","endcap_photonet","endcap_photoneta","endcap_photonsceta","endcap_jet1pt","endcap_jet1eta","endcap_jet2pt","endcap_jet2eta","endcap_Mjj","endcap_zepp","endcap_deltaeta","endcap_MET_et","endcap_Dphiwajj","endcap_Mla"};

    //file prepare
    //MC: ZA
    TFile * za1 = TFile::Open("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/hist_store/hist_ZG.root");
    //MC: VV
    TFile * vv1 = TFile::Open("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/hist_store/hist_VV.root");
    //MC: TTA
    TFile * tta1 = TFile::Open("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/hist_store/hist_TTA.root");
    //MC: STop
    TFile * st1 = TFile::Open("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/hist_store/hist_ST.root");
    //MC: WGJJ (signal)
    TFile * wgjj1 = TFile::Open("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/hist_store/hist_WAJJ.root");
    //MC: WGJets
    TFile * wgjets1 = TFile::Open("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/hist_store/hist_WA.root");

    //DATA
    TFile * data1 = TFile::Open("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/hist_store/hist_data_mu.root");
    //fake photon data 
    TFile * data2 = TFile::Open("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/hist_store/hist_data_fakephoton_mu.root");
    //fake lepton data
    TFile * data3 = TFile::Open("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/hist_store/hist_data_fake_mu.root");
    //double fake data
    TFile * data4 = TFile::Open("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/hist_store/hist_data_doublefake_mu.root");

    //output file to save result (a lot of canvas)
    TFile* output = new TFile("DATAvsMC/result.root","RECREATE");

    //Draw Hist
    for(int i = 1; i<=36; i++){
	//set hist x title
	if(string[i-1].Contains("nVtx"))                xtitle = "nVtx";
    	if(string[i-1].Contains("ptlep1"))              xtitle = "ptlep1 [GeV]";
	if(string[i-1].Contains("etalep1"))             xtitle = "etalep1";
        if(string[i-1].Contains("mtVlepJECnew"))        xtitle = "mtVlepJECnew [GeV]";
        if(string[i-1].Contains("ptVlepJEC"))           xtitle = "ptVlepJEC [GeV]";
        if(string[i-1].Contains("photonet"))            xtitle = "photonet [GeV]";
        if(string[i-1].Contains("photoneta"))           xtitle = "photoneta";
        if(string[i-1].Contains("photonsceta"))         xtitle = "photonsceta";
        if(string[i-1].Contains("jet1pt"))              xtitle = "jet1pt [GeV]";
        if(string[i-1].Contains("jet1eta"))             xtitle = "jet1eta";
        if(string[i-1].Contains("jet2pt"))              xtitle = "jet2pt [GeV]";
        if(string[i-1].Contains("jet2eta"))             xtitle = "jet2eta";
        if(string[i-1].Contains("Mjj"))                 xtitle = "Mjj [GeV]";
        if(string[i-1].Contains("zepp"))                xtitle = "zepp";
        if(string[i-1].Contains("MET_et"))              xtitle = "MET_et [GeV]";
        if(string[i-1].Contains("deltaeta"))            xtitle = "deltaeta";
        if(string[i-1].Contains("Dphiwajj"))            xtitle = "Dphiwajj";
        if(string[i-1].Contains("Mla"))                 xtitle = "Mla [GeV]";

	//set canvas
	TCanvas* cv = new TCanvas("cv_" + string[i-1], "cv_" + string[i-1], 700, 700);
        TPad* fPads1 = new TPad("pad1", "", 0.00, 0.25, 1.00, 1.00);
        TPad* fPads2 = new TPad("pad2", "", 0.00, 0.00, 1.00, 0.25);
        fPads1->SetFillColor(0);
        fPads1->SetLineColor(0);
        fPads2->SetFillColor(0);
        fPads2->SetLineColor(0);
        fPads1->SetBottomMargin(0);
        fPads2->SetTopMargin(0);
        fPads2->SetBottomMargin(0.4);
        fPads1->Draw();
        fPads2->Draw();

	//set stack
	THStack* hs = new THStack(string[i-1], string[i-1]);


	//setting for hist
	//MC: ZA
	TH1F* ZA1 = (TH1F*)za1->Get(string[i-1]);
        TH1F* ZA = (TH1F*)ZA1->Clone();
        ZA->SetFillColor(600);
        ZA->SetLineColor(0);
        ZA->SetLineWidth(0);
	//MC: VV
	TH1F* VV1 = (TH1F*)vv1->Get(string[i-1]);
	TH1F* VV = (TH1F*)VV1->Clone();
        VV->SetFillColor(616);
        VV->SetLineColor(0);
        VV->SetLineWidth(0);
	//MC: TTA
	TH1F* TTA1 = (TH1F*)tta1->Get(string[i-1]);
	TH1F* TTA = (TH1F*)TTA1->Clone();
        TTA->SetFillColor(416);
        TTA->SetLineColor(0);
        TTA->SetLineWidth(0);
	//MC: STop
	TH1F* ST1 = (TH1F*)st1->Get(string[i-1]);
	TH1F* ST = (TH1F*)ST1->Clone();
	ST->SetFillColor(797);
        ST->SetLineColor(0);
        ST->SetLineWidth(0);
	//MC: WGJJ (signal)
	TH1F* WGJJ1 = (TH1F*)wgjj1->Get(string[i-1]);
	TH1F* WGJJ = (TH1F*)WGJJ1->Clone();
	WGJJ->SetFillColor(432);
        WGJJ->SetLineColor(0);
       	WGJJ->SetLineWidth(0);
	//MC: WGJets
	TH1F* WGJETS1 = (TH1F*)wgjets1->Get(string[i-1]);
	TH1F* WGJETS = (TH1F*)WGJETS1->Clone();       
	WGJETS->SetFillColor(632);
        WGJETS->SetLineColor(0);
       	WGJETS->SetLineWidth(0);
	//DATA: fake photon
        TH1F* DATA2 = (TH1F*)data2->Get(string[i-1]);
        TH1F* fp = (TH1F*)DATA2->Clone();
        fp->SetFillColor(400);
        fp->SetLineColor(0);
        fp->SetLineWidth(0);
	//DATA: fake lepton
        TH1F* DATA3 = (TH1F*)data3->Get(string[i-1]);
        TH1F* fl = (TH1F*)DATA3->Clone();
        fl->SetFillColor(419);
        fl->SetLineColor(0);
        fl->SetLineWidth(0);
        //DATA: double fake
        TH1F* DATA4 = (TH1F*)data4->Get(string[i-1]);
        TH1F* df = (TH1F*)DATA4->Clone();
        df->SetFillColor(620);
        df->SetLineColor(0);
        df->SetLineWidth(0);

	//fp->Add(df,-1);
	//fl->Add(df,-1);	

	//DATA
        TH1F* DATA1 = (TH1F*)data1->Get(string[i-1]);
        TH1F* DATA = (TH1F*)DATA1->Clone();
        DATA->SetLineColor(kBlack);
        DATA->SetLineWidth(2);
        DATA->GetXaxis()->SetTitle(xtitle);
        DATA->GetYaxis()->SetTitle("Events/bin");
        DATA->SetMarkerStyle(20);
        DATA->SetMarkerSize(0.8);
        DATA->SetStats(kFALSE);

	//add hist of MC and Data, and draw them
	TH1F* hmc = (TH1F*)ZA->Clone();
	hmc->Reset();
	if( addZA==1) 		hmc->Add(ZA,1);
        if( addVV==1) 		hmc->Add(VV,1);
        if( addTTA==1) 		hmc->Add(TTA,1);
        if( addSTop==1) 	hmc->Add(ST,1);
        if( addWGJJ==1) 	hmc->Add(WGJJ,1);
        if( addWGJets==1) 	hmc->Add(WGJETS,1);
        if( addFakePhoton==1)   hmc->Add(fp,1);
        if( addFakeLepton==1)   hmc->Add(fl,1);
        if( addDoubleFake==1)    hmc->Add(df,1);

	//stack the hist
        TH1F* h = (TH1F*)ZA->Clone();
	h->Reset();
	hs->Add(h);

	if( addZA==1)		hs->Add(ZA);
       	if( addVV==1)		hs->Add(VV);
        if( addTTA==1)		hs->Add(TTA);
        if( addSTop==1)		hs->Add(ST);
	if( addWGJJ==1)		hs->Add(WGJJ);
        if( addDoubleFake==1)   hs->Add(df);
        if( addFakeLepton==1)   hs->Add(fl);
        if( addFakePhoton==1)   hs->Add(fp);
	if( addWGJets==1)	hs->Add(WGJETS);
     


        //entries
        int data_e=DATA->GetEntries();
	TString data_entry = to_string(data_e);
        int df_e=df->GetEntries();
        TString df_entry = to_string(df_e);
        int fl_e=fl->GetEntries();
        TString fl_entry = to_string(fl_e);
        int fp_e=fp->GetEntries();
        TString fp_entry = to_string(fp_e);

        int za_e=ZA->GetEntries();
        TString za_entry = to_string(za_e);
        int vv_e=VV->GetEntries();
        TString vv_entry = to_string(vv_e);
        int tta_e=TTA->GetEntries();
        TString tta_entry = to_string(tta_e);
        int st_e=ST->GetEntries();
        TString st_entry = to_string(st_e);
        int wgjj_e=WGJJ->GetEntries();
        TString wgjj_entry = to_string(wgjj_e);
        int wgjets_e=WGJETS->GetEntries();
        TString wgjets_entry = to_string(wgjets_e);

	//set legend
	TLegend* leg = new TLegend(0.58,0.6,0.88,0.89);
	leg->SetBorderSize(0);
        if( addDATA==1)         leg->AddEntry(DATA,  "DATA		[ entries = " + data_entry + " ]","lp");
        if( addZA==1)		leg->AddEntry(ZA,    "ZA      		[ entries = " + za_entry + " ]","f");
        if( addVV==1)		leg->AddEntry(VV,    "VV              	[ entries = " + vv_entry + " ]","f");
        if( addTTA==1)		leg->AddEntry(TTA,   "TTA             	[ entries = " + tta_entry + " ]","f");
        if( addSTop==1)		leg->AddEntry(ST,    "STop            	[ entries = " + st_entry + " ]","f");
	if( addWGJJ==1)		leg->AddEntry(WGJJ,  "WGJJ            	[ entries = " + wgjj_entry + " ]","f");
	if( addWGJets==1)	leg->AddEntry(WGJETS,"WG+jets        	[ entries = " + wgjets_entry + " ]","f");
        if( addDoubleFake==1)   leg->AddEntry(df,    "Double Fake	[ entries = " + df_entry + " ]","f");
        if( addFakePhoton==1)   leg->AddEntry(fp,    "Fake Photon 	[ entries = " + fp_entry + " ]","f");
        if( addFakeLepton==1)   leg->AddEntry(fl,    "Fake Lepton 	[ entries = " + fl_entry + " ]","f");

	//DATA and MC pad
	fPads1->cd();
        hs->Draw("HIST");
      	hs->GetYaxis()->SetTitle("Events/bin");
        hs->GetYaxis()->SetTitleSize(0.04);
        hs->GetYaxis()->SetTitleOffset(1.3);

     	double maximumMC       = 1.6 * hmc->GetMaximum();
      	double maximumDATA     = 1.6 * DATA->GetMaximum();
       	double maximumForStack = (maximumMC > maximumDATA ? maximumMC : maximumDATA);
        hs->SetMaximum(maximumForStack);

	if( addDATA==1) DATA->Draw("SAME APE");
	leg->Draw();


        //DATA DIVIDE MC
        fPads2->cd();
       	fPads2->SetGridy();
	TH1F* divide = (TH1F*)DATA->Clone();
	divide->Divide(hmc);
        divide->SetTitle("");
       	divide->SetStats(kFALSE);
      	divide->SetLineColor(kBlue);
        divide->SetLineWidth(2);
        divide->GetYaxis()->SetTitle("Data/MC");
        divide->GetYaxis()->CenterTitle();
        divide->GetYaxis()->SetTitleSize(0.1);
       	divide->GetYaxis()->SetTitleOffset(0.3);
        divide->GetYaxis()->SetLabelSize(0.07);
        divide->GetXaxis()->SetTitle(xtitle);
        divide->GetXaxis()->SetLabelSize(0.07);
       	divide->GetXaxis()->SetTitleSize(0.07);
	divide->Draw("HIST P e");
	divide->Draw();

	cv->SaveAs("DATAvsMC/"+string[i-1]+".png");
	output->cd();
	cv->Write();
	//output->Write();
    }

    output->Close();



}   

int main() {

    int addZA, addVV, addTTA, addSTop, addWGJJ, addWGJets;
    int addDATA, addFakePhoton, addFakeLepton, addDoubleFake;
    //code help
    cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
    cout<<"<<                                                                                  <<"<<endl;
    cout<<"<<  This program is used to add hist, you should tell it which hist should be add.  <<"<<endl;
    cout<<"<<  The hist sets to 1, means to add.                                               <<"<<endl;
    cout<<"<<  There are two kinds of hist.                                                    <<"<<endl;
    cout<<"<<  The hists of MC are   : addZA, addVV, addTTA, addSTop, addWGJJ, addWGJets       <<"<<endl;
    cout<<"<<  The hists of DATA are : addDATA, addFakePhoton, addDoubleFake                   <<"<<endl;
    cout<<"<<                                                                                  <<"<<endl;
    cout<<"<<                                                               -- Jing Peng       <<"<<endl;
    cout<<"<<                                                                                  <<"<<endl;
    cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl<<endl;
    
    cout<<"Starting adding hists................................................................."<<endl<<endl;
    
    //for MC
    cout<<"Prepare hists of MC:"<<endl;
    cout<<"addZA        = ";
    cin>>addZA;

    cout<<"addVV        = ";
    cin>>addVV;

    cout<<"addTTA       = ";
    cin>>addTTA;

    cout<<"addSTop      = ";
    cin>>addSTop;

    cout<<"addWGJJ      = ";
    cin>>addWGJJ;

    cout<<"addWGJets    = ";
    cin>>addWGJets;
    cout<<endl;

    //for DATA
    cout<<"Prepare hists of DATA:"<<endl;
    cout<<"addDATA         = ";
    cin>>addDATA;

    cout<<"addFakePhoton   = ";
    cin>>addFakePhoton;

    cout<<"addFakeLepton   = ";
    cin>>addFakeLepton;

    cout<<"addDoubleFake   = ";
    cin>>addDoubleFake;

    cout<<endl;
    addhist(addZA, addVV, addTTA, addSTop, addWGJJ, addWGJets, addDATA, addFakePhoton, addFakeLepton, addDoubleFake);

    cout<<endl;
    cout<<"Ending everything, thanks for your using."<<endl;
    cout<<"See you next time!"<<endl;

    return 0;
}

