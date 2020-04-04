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
void hist_data_fake(){

	double fake_lepton_weight, weight;

        TFile * input = TFile::Open("/home/pku/pengj/VBSWGamma/94X2016/myplot/plotdata/fakelepton_mu.root");
        TFile* output = new TFile("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/hist_store/hist_data_fake_mu.root","RECREATE");

	//fake lepton weight
	TFile * frw = TFile::Open("/home/pku/pengj/VBSWGamma/94X2016/fake_lepton_rate/fake_muon.root");
	TH2D* h = (TH2D*)frw->Get("weight");

    	TCanvas* cb1 = new TCanvas("cb1", "cb1", 700, 700);
        TCanvas* cb2 = new TCanvas("cb2", "cb2", 700, 700);
        TCanvas* cb3 = new TCanvas("cb3", "cb3", 700, 700);
        TCanvas* cb4 = new TCanvas("cb4", "cb4", 700, 700);
        TCanvas* cb5 = new TCanvas("cb5", "cb5", 700, 700);
        TCanvas* cb6 = new TCanvas("cb6", "cb6", 700, 700);
        TCanvas* cb7 = new TCanvas("cb7", "cb7", 700, 700);
        TCanvas* cb8 = new TCanvas("cb8", "cb8", 700, 700);
        TCanvas* cb9 = new TCanvas("cb9", "cb9", 700, 700);
        TCanvas* cb10 = new TCanvas("cb10", "cb10", 700, 700);
        TCanvas* cb11 = new TCanvas("cb11", "cb11", 700, 700);
        TCanvas* cb12 = new TCanvas("cb12", "cb12", 700, 700);
        TCanvas* cb13 = new TCanvas("cb13", "cb13", 700, 700);
        TCanvas* cb14 = new TCanvas("cb14", "cb14", 700, 700);
        TCanvas* cb15 = new TCanvas("cb15", "cb15", 700, 700);
        TCanvas* cb16 = new TCanvas("cb16", "cb16", 700, 700);
        TCanvas* cb17 = new TCanvas("cb17", "cb17", 700, 700);
        TCanvas* cb18 = new TCanvas("cb18", "cb18", 700, 700);

        TCanvas* ce1 = new TCanvas("ce1", "ce1", 700, 700);
        TCanvas* ce2 = new TCanvas("ce2", "ce2", 700, 700);
        TCanvas* ce3 = new TCanvas("ce3", "ce3", 700, 700);
        TCanvas* ce4 = new TCanvas("ce4", "ce4", 700, 700);
        TCanvas* ce5 = new TCanvas("ce5", "ce5", 700, 700);
        TCanvas* ce6 = new TCanvas("ce6", "ce6", 700, 700);
        TCanvas* ce7 = new TCanvas("ce7", "ce7", 700, 700);
        TCanvas* ce8 = new TCanvas("ce8", "ce8", 700, 700);
        TCanvas* ce9 = new TCanvas("ce9", "ce9", 700, 700);
        TCanvas* ce10 = new TCanvas("ce10", "ce10", 700, 700);
        TCanvas* ce11 = new TCanvas("ce11", "ce11", 700, 700);
        TCanvas* ce12 = new TCanvas("ce12", "ce12", 700, 700);
        TCanvas* ce13 = new TCanvas("ce13", "ce13", 700, 700);
        TCanvas* ce14 = new TCanvas("ce14", "ce14", 700, 700);
        TCanvas* ce15 = new TCanvas("ce15", "ce15", 700, 700);
        TCanvas* ce16 = new TCanvas("ce16", "ce16", 700, 700);
        TCanvas* ce17 = new TCanvas("ce17", "ce17", 700, 700);
        TCanvas* ce18 = new TCanvas("ce18", "ce18", 700, 700);

	//barrel hist
	TH1F* barrel_nVtx = new TH1F("barrel_nVtx", "barrel_nVtx", 30, 0, 50);
        TH1F* barrel_ptlep1 = new TH1F("barrel_ptlep1", "barrel_ptlep1", 30, 25.0, 200.0);
        TH1F* barrel_etalep1 = new TH1F("barrel_etalep1", "barrel_etalep1", 30, -2.4, 2.4);
        TH1F* barrel_mtVlepJECnew = new TH1F("barrel_mtVlepJECnew", "barrel_mtVlepJECnew", 30, 35.0, 100.0);
        TH1F* barrel_ptVlepJEC = new TH1F("barrel_ptVlepJEC", "barrel_ptVlepJEC", 30, 0.0, 200.0);
        TH1F* barrel_photonet = new TH1F("barrel_photonet", "barrel_photonet", 30, 30.0, 200.0);
        TH1F* barrel_photoneta = new TH1F("barrel_photoneta", "barrel_photoneta", 30, -2.5, 2.5);
        TH1F* barrel_photonsceta = new TH1F("barrel_photonsceta", "barrel_photonsceta", 30, -2.5, 2.5);
        TH1F* barrel_jet1pt = new TH1F("barrel_jet1pt", "barrel_jet1pt", 30, 40.0, 300.0);
        TH1F* barrel_jet1eta = new TH1F("barrel_jet1eta", "barrel_jet1eta", 30, -5.0, 5.0);
        TH1F* barrel_jet2pt = new TH1F("barrel_jet2pt", "barrel_jet2pt", 30, 30.0, 200.0);
        TH1F* barrel_jet2eta = new TH1F("barrel_jet2eta", "barrel_jet2eta", 30, -5.0, 5.0);
        TH1F* barrel_Mjj = new TH1F("barrel_Mjj", "barrel_Mjj", 30, 200., 400.);
        TH1F* barrel_zepp = new TH1F("barrel_zepp", "barrel_zepp", 30, 0.0, 5.0);
        TH1F* barrel_deltaeta = new TH1F("barrel_deltaeta", "barrel_deltaeta", 30, 0.5, 5.0);
        TH1F* barrel_MET_et = new TH1F("barrel_MET_et", "barrel_MET_et", 30, 30.0, 200.0);
        TH1F* barrel_Dphiwajj = new TH1F("barrel_Dphiwajj", "barrel_Dphiwajj", 30, 0.0, 3.2);
        TH1F* barrel_Mla = new TH1F("barrel_Mla", "barrel_Mla", 30, 0.0, 200);
	//end cap
        TH1F* endcap_nVtx = new TH1F("endcap_nVtx", "endcap_nVtx", 30, 0, 50);
        TH1F* endcap_ptlep1 = new TH1F("endcap_ptlep1", "endcap_ptlep1", 30, 25.0, 200.0);
        TH1F* endcap_etalep1 = new TH1F("endcap_etalep1", "endcap_etalep1", 30, -2.4, 2.4);
        TH1F* endcap_mtVlepJECnew = new TH1F("endcap_mtVlepJECnew", "endcap_mtVlepJECnew", 30, 35.0, 100.0);
        TH1F* endcap_ptVlepJEC = new TH1F("endcap_ptVlepJEC", "endcap_ptVlepJEC", 30, 0.0, 200.0);
        TH1F* endcap_photonet = new TH1F("endcap_photonet", "endcap_photonet", 30, 30.0, 200.0);
        TH1F* endcap_photoneta = new TH1F("endcap_photoneta", "endcap_photoneta", 30, -2.5, 2.5);
        TH1F* endcap_photonsceta = new TH1F("endcap_photonsceta", "endcap_photonsceta", 30, -2.5, 2.5);
        TH1F* endcap_jet1pt = new TH1F("endcap_jet1pt", "endcap_jet1pt", 30, 40.0, 300.0);
        TH1F* endcap_jet1eta = new TH1F("endcap_jet1eta", "endcap_jet1eta", 30, -5.0, 5.0);
        TH1F* endcap_jet2pt = new TH1F("endcap_jet2pt", "endcap_jet2pt", 30, 30.0, 200.0);
        TH1F* endcap_jet2eta = new TH1F("endcap_jet2eta", "endcap_jet2eta", 30, -5.0, 5.0);
        TH1F* endcap_Mjj = new TH1F("endcap_Mjj", "endcap_Mjj", 30, 200., 400.);
        TH1F* endcap_zepp = new TH1F("endcap_zepp", "endcap_zepp", 30, 0.0, 5.0);
        TH1F* endcap_deltaeta = new TH1F("endcap_deltaeta", "endcap_deltaeta", 30, 0.5, 5.0);
        TH1F* endcap_MET_et = new TH1F("endcap_MET_et", "endcap_MET_et", 30, 30.0, 200.0);
        TH1F* endcap_Dphiwajj = new TH1F("endcap_Dphiwajj", "endcap_Dphiwajj", 30, 0.0, 3.2);
        TH1F* endcap_Mla = new TH1F("endcap_Mla", "endcap_Mla", 30, 0.0, 200);

        TTreeReader fReader("demo", input);
        //TTreeReaderValue<Double_t> scalef              = {fReader, "scalef"};
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
        //TTreeReaderValue<Int_t>    lep1_sign           = {fReader, "lep1_sign"};
        //TTreeReaderValue<Int_t>    muon1_trackerLayers = {fReader, "muon1_trackerLayers"};
        //TTreeReaderValue<Double_t> matchedgenMu1_pt    = {fReader, "matchedgenMu1_pt"};
        TTreeReaderValue<Double_t> photonsceta         = {fReader, "photonsceta"};
        TTreeReaderValue<Double_t> photonscphi         = {fReader, "photonscphi"};

        TTreeReaderValue<Double_t> Mla                 = {fReader, "Mla"};
        TTreeReaderValue<Bool_t>   photonhaspixelseed  = {fReader, "photonhaspixelseed"};
	TTreeReaderValue<Bool_t>   photonhaspixelseed_f  = {fReader, "photonhaspixelseed_f"};
        TTreeReaderValue<Bool_t>   photonpasseleveto   = {fReader, "photonpasseleveto"};
        TTreeReaderValue<Double_t> photonet            = {fReader, "photonet"};
        TTreeReaderValue<Double_t> photoneta           = {fReader, "photoneta"};
        TTreeReaderValue<Double_t> photonphi           = {fReader, "photonphi"};
        TTreeReaderValue<Double_t> photone             = {fReader, "photone"};
        //TTreeReaderValue<Double_t> photonsieie         = {fReader, "photonsieie"};
        //TTreeReaderValue<Double_t> photonphoiso        = {fReader, "photonphoiso"};
        //TTreeReaderValue<Double_t> photonchiso         = {fReader, "photonchiso"};
        //TTreeReaderValue<Double_t> photonnhiso         = {fReader, "photonnhiso"};
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
        TTreeReaderValue<Double_t> Dphiwajj            = {fReader, "Dphiwajj"};

        TTreeReaderValue<Double_t> j1metPhi_f            = {fReader, "j1metPhi_f"};
        TTreeReaderValue<Double_t> j2metPhi_f            = {fReader, "j2metPhi_f"};
        TTreeReaderValue<Double_t> drj1a_f               = {fReader, "drj1a_f"};
        TTreeReaderValue<Double_t> drj2a_f               = {fReader, "drj2a_f"};
        TTreeReaderValue<Double_t> drj1l_f               = {fReader, "drj1l_f"};
        TTreeReaderValue<Double_t> drj2l_f               = {fReader, "drj2l_f"};
        TTreeReaderValue<Double_t> jet1icsv_f            = {fReader, "jet1icsv_f"};
        TTreeReaderValue<Double_t> jet2icsv_f            = {fReader, "jet2icsv_f"};
        TTreeReaderValue<Double_t> deltaeta_f            = {fReader, "deltaeta_f"};
        TTreeReaderValue<Double_t> photonet_f            = {fReader, "photonet_f"};
        TTreeReaderValue<Double_t> drla_f                = {fReader, "drla"};
        //TTreeReaderValue<Double_t> MET_et_f              = {fReader, "MET_et"};

        Long64_t maxEntries = fReader.GetEntries(false);
	cout << "Number of events to be analyzed : " << maxEntries << std::endl;
	int i =0;
	double etal1, ptl1;
	while (fReader.Next()){
		if (i % 10000 == 0){
			cout<<i*100/maxEntries<<"%"<< "\r" << std::flush;
		}
		i++;

      
                if(!(*HLT_Mu2==1
                && *Mjj>200. && *Mjj<400.
                && *jet1pt>40. && abs(*jet1eta)<4.7 && *jet2pt>30. && abs(*jet2eta)<4.7
                && abs(*lep)==13 && *ptlep1>30. && abs(*etalep1)<2.4 && *ngoodmus==1 && *ngoodeles==0
                && *photonhaspixelseed==0 && *photonet>25. && *photonet < 400.
                && *MET_et>30. && *mtVlepJECnew>30
                && *drla>0.5 && *drj1l>0.5 && *drj2l>0.5 && *drj1a>0.5 && *drj2a>0.5
                && abs(*j1metPhi)>0.5 && abs(*j2metPhi)>0.5
		&& *jet1icsv < 0.8484 && *jet2icsv < 0.8484 && *jet1icsv > -10. && *jet2icsv > -10.)) continue;


		etal1 = *etalep1;
		if (*ptlep1 >= 50) ptl1 = 45;
		else ptl1 = *ptlep1;
		fake_lepton_weight = h->GetBinContent(h->GetXaxis()->FindBin(fabs(etal1)),h->GetYaxis()->FindBin(ptl1));
		weight = fake_lepton_weight;
		//barell hist fill
		if (fabs(*photonsceta) < 1.4442){
			barrel_nVtx->Fill(*nVtx,weight);
			barrel_ptlep1->Fill(*ptlep1,weight);
			barrel_etalep1->Fill(*etalep1,weight);
			barrel_mtVlepJECnew->Fill(*mtVlepJECnew,weight);
			barrel_ptVlepJEC->Fill(*ptVlepJEC,weight);
			barrel_photonet->Fill(*photonet,weight);
			barrel_photoneta->Fill(*photoneta,weight);
			barrel_photonsceta->Fill(*photonsceta,weight);
			barrel_jet1pt->Fill(*jet1pt,weight);
			barrel_jet1eta->Fill(*jet1eta,weight);
                        barrel_jet2pt->Fill(*jet2pt,weight);
                        barrel_jet2eta->Fill(*jet2eta,weight);
			barrel_Mjj->Fill(*Mjj,weight);
			barrel_zepp->Fill(*zepp,weight);
			barrel_deltaeta->Fill(*deltaeta,weight);
			barrel_MET_et->Fill(*MET_et,weight);
			barrel_Dphiwajj->Fill(*Dphiwajj,weight);
			barrel_Mla->Fill(*Mla,weight);

		}

                //endcap hist fill
                if (fabs(*photonsceta) > 1.566 && fabs(*photonsceta) < 2.5){
                        endcap_nVtx->Fill(*nVtx,weight);
                        endcap_ptlep1->Fill(*ptlep1,weight);
                        endcap_etalep1->Fill(*etalep1,weight);
                        endcap_mtVlepJECnew->Fill(*mtVlepJECnew,weight);
                        endcap_ptVlepJEC->Fill(*ptVlepJEC,weight);
                        endcap_photonet->Fill(*photonet,weight);
                        endcap_photoneta->Fill(*photoneta,weight);
                        endcap_photonsceta->Fill(*photonsceta,weight);
                        endcap_jet1pt->Fill(*jet1pt,weight);
                        endcap_jet1eta->Fill(*jet1eta,weight);
                        endcap_jet2pt->Fill(*jet2pt,weight);
                        endcap_jet2eta->Fill(*jet2eta,weight);
                        endcap_Mjj->Fill(*Mjj,weight);
                        endcap_zepp->Fill(*zepp,weight);
                        endcap_deltaeta->Fill(*deltaeta,weight);
                        endcap_MET_et->Fill(*MET_et,weight);
                        endcap_Dphiwajj->Fill(*Dphiwajj,weight);
                        endcap_Mla->Fill(*Mla,weight);

                }

	}

        output->cd();

        barrel_nVtx->Write();
        barrel_ptlep1->Write();
	barrel_etalep1->Write();
      	barrel_mtVlepJECnew->Write();
       	barrel_ptVlepJEC->Write();
       	barrel_photonet->Write();
        barrel_photoneta->Write();
       	barrel_photonsceta->Write();
       	barrel_jet1pt->Write();
        barrel_jet1eta->Write();
       	barrel_jet2pt->Write();
        barrel_jet2eta->Write();
        barrel_Mjj->Write();
      	barrel_zepp->Write();
       	barrel_deltaeta->Write();
        barrel_MET_et->Write();
       	barrel_Dphiwajj->Write();
        barrel_Mla->Write();

	cb1->cd();
        barrel_nVtx->Draw();
	cb2->cd();
        barrel_ptlep1->Draw();
        cb3->cd();
        barrel_etalep1->Draw();
        cb4->cd();
        barrel_mtVlepJECnew->Draw();
        cb5->cd();
        barrel_ptVlepJEC->Draw();
        cb6->cd();
        barrel_photonet->Draw();
        cb7->cd();
        barrel_photoneta->Draw();
        cb8->cd();
        barrel_photonsceta->Draw();
        cb9->cd();
        barrel_jet1pt->Draw();
        cb10->cd();
        barrel_jet1eta->Draw();
        cb11->cd();
        barrel_jet2pt->Draw();
        cb12->cd();
        barrel_jet2eta->Draw();
        cb13->cd();
        barrel_Mjj->Draw();
        cb14->cd();
        barrel_zepp->Draw();
        cb15->cd();
        barrel_deltaeta->Draw();
        cb16->cd();
        barrel_MET_et->Draw();
        cb17->cd();
        barrel_Dphiwajj->Draw();
        cb18->cd();
        barrel_Mla->Draw();

        endcap_nVtx->Write();
        endcap_ptlep1->Write();
        endcap_etalep1->Write();
        endcap_mtVlepJECnew->Write();
        endcap_ptVlepJEC->Write();
        endcap_photonet->Write();
        endcap_photoneta->Write();
        endcap_photonsceta->Write();
        endcap_jet1pt->Write();
        endcap_jet1eta->Write();
        endcap_jet2pt->Write();
        endcap_jet2eta->Write();
        endcap_Mjj->Write();
        endcap_zepp->Write();
        endcap_deltaeta->Write();
        endcap_MET_et->Write();
        endcap_Dphiwajj->Write();
        endcap_Mla->Write();

        ce1->cd();
        endcap_nVtx->Draw();
        ce2->cd();
        endcap_ptlep1->Draw();
        ce3->cd();
        endcap_etalep1->Draw();
        ce4->cd();
        endcap_mtVlepJECnew->Draw();
        ce5->cd();
        endcap_ptVlepJEC->Draw();
        ce6->cd();
        endcap_photonet->Draw();
        ce7->cd();
        endcap_photoneta->Draw();
        ce8->cd();
        endcap_photonsceta->Draw();
        ce9->cd();
        endcap_jet1pt->Draw();
        ce10->cd();
        endcap_jet1eta->Draw();
        ce11->cd();
        endcap_jet2pt->Draw();
        ce12->cd();
        endcap_jet2eta->Draw();
        ce13->cd();
        endcap_Mjj->Draw();
        ce14->cd();
        endcap_zepp->Draw();
        ce15->cd();
        endcap_deltaeta->Draw();
        ce16->cd();
        endcap_MET_et->Draw();
        ce17->cd();
        endcap_Dphiwajj->Draw();
        ce18->cd();
        endcap_Mla->Draw();
        //output->Write();
        //output->Close();

        cb1->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_nVtx.png");
        cb2->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_ptlep1.png");
        cb3->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_etalep1.png");
        cb4->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_mtVlepJECnew.png");
        cb5->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_ptVlepJEC.png");
        cb6->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_photonet.png");
        cb7->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_photoneta.png");
        cb8->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_photonsceta.png");
        cb9->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_jet1pt.png");
        cb10->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_jet1eta.png");
        cb11->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_jet2pt.png");
        cb12->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_jet2eta.png");
        cb13->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_Mjj.png");
        cb14->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_zepp.png");
        cb15->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_deltaeta.png");
        cb16->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_MET_et.png");
        cb17->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_Dphiwajj.png");
        cb18->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_barrel_Mla.png");

        ce1->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_nVtx.png");
        ce2->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_ptlep1.png");
        ce3->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_etalep1.png");
        ce4->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_mtVlepJECnew.png");
        ce5->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_ptVlepJEC.png");
        ce6->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_photonet.png");
        ce7->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_photoneta.png");
        ce8->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_photonsceta.png");
        ce9->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_jet1pt.png");
        ce10->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_jet1eta.png");
        ce11->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_jet2pt.png");
        ce12->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_jet2eta.png");
        ce13->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_Mjj.png");
        ce14->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_zepp.png");
        ce15->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_deltaeta.png");
        ce16->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_MET_et.png");
        ce17->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_Dphiwajj.png");
        ce18->SaveAs("/home/pku/pengj/VBSWGamma/94X2016/myplot/make_hist/his_tore/png/data_fake_endcap_Mla.png");

        output->Write();
        output->Close();
}
