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
void hist_data_fakephoton(){

	double fake_lepton_weight, barrel_fake_photon_weight, endcap_fake_photon_weight, weight;

        TFile * input = TFile::Open("fakephoton_mu.root");
        TFile* output = new TFile("hist_fakephoton_photonet_mu.root","RECREATE");


    	//TCanvas* cb1 = new TCanvas("cb1", "cb1", 700, 700);

        //TCanvas* ce1 = new TCanvas("ce1", "ce1", 700, 700);


	const Int_t NBINS = 7;
	Double_t edges[NBINS + 1] = {25., 30., 40., 50., 70., 100., 135., 400.};
	//barrel hist
        TH1F* barrel_photonet = new TH1F("barrel_photonet", "barrel_photonet", NBINS, edges);
	//end cap
        TH1F* endcap_photonet = new TH1F("endcap_photonet", "endcap_photonet", NBINS, edges);

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
        TTreeReaderValue<Double_t> photonsceta_f         = {fReader, "photonsceta_f"};

        TTreeReaderValue<Double_t> photonscphi         = {fReader, "photonscphi"};

        TTreeReaderValue<Double_t> Mla                 = {fReader, "Mla"};
        TTreeReaderValue<Double_t> Mla_f                 = {fReader, "Mla_f"};

        TTreeReaderValue<Bool_t>   photonhaspixelseed  = {fReader, "photonhaspixelseed"};
	TTreeReaderValue<Bool_t>   photonhaspixelseed_f  = {fReader, "photonhaspixelseed_f"};
        TTreeReaderValue<Bool_t>   photonpasseleveto   = {fReader, "photonpasseleveto"};
        TTreeReaderValue<Double_t> photonet            = {fReader, "photonet"};
        TTreeReaderValue<Double_t> photoneta           = {fReader, "photoneta"};
        TTreeReaderValue<Double_t> photonet_f            = {fReader, "photonet_f"};
        TTreeReaderValue<Double_t> photoneta_f           = {fReader, "photoneta_f"};
        TTreeReaderValue<Double_t> photonphi           = {fReader, "photonphi"};
        TTreeReaderValue<Double_t> photone             = {fReader, "photone"};
        //TTreeReaderValue<Double_t> photonsieie         = {fReader, "photonsieie"};
        //TTreeReaderValue<Double_t> photonphoiso        = {fReader, "photonphoiso"};
        //TTreeReaderValue<Double_t> photonchiso         = {fReader, "photonchiso"};
        //TTreeReaderValue<Double_t> photonnhiso         = {fReader, "photonnhiso"};
        TTreeReaderValue<Double_t> drla                = {fReader, "drla"};
        TTreeReaderValue<Double_t> jet1pt              = {fReader, "jet1pt"};
        TTreeReaderValue<Double_t> jet1eta             = {fReader, "jet1eta"};
        TTreeReaderValue<Double_t> jet1pt_f              = {fReader, "jet1pt_f"};
        TTreeReaderValue<Double_t> jet1eta_f             = {fReader, "jet1eta_f"};
        TTreeReaderValue<Double_t> jet1phi             = {fReader, "jet1phi"};
        TTreeReaderValue<Double_t> jet1e               = {fReader, "jet1e"};
        TTreeReaderValue<Double_t> jet1icsv            = {fReader, "jet1icsv"};
        TTreeReaderValue<Double_t> jet2pt              = {fReader, "jet2pt"};
        TTreeReaderValue<Double_t> jet2eta             = {fReader, "jet2eta"};
        TTreeReaderValue<Double_t> jet2pt_f              = {fReader, "jet2pt_f"};
        TTreeReaderValue<Double_t> jet2eta_f             = {fReader, "jet2eta_f"};
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
        TTreeReaderValue<Double_t> Mjj_f                 = {fReader, "Mjj_f"};
        TTreeReaderValue<Double_t> deltaeta_f            = {fReader, "deltaeta_f"};
        TTreeReaderValue<Double_t> zepp_f                = {fReader, "zepp_f"};
        TTreeReaderValue<Double_t> j1metPhi            = {fReader, "j1metPhi"};
        TTreeReaderValue<Double_t> j2metPhi            = {fReader, "j2metPhi"};
        TTreeReaderValue<Double_t> Dphiwajj            = {fReader, "Dphiwajj"};
        TTreeReaderValue<Double_t> Dphiwajj_f            = {fReader, "Dphiwajj_f"};

        TTreeReaderValue<Double_t> j1metPhi_f            = {fReader, "j1metPhi_f"};
        TTreeReaderValue<Double_t> j2metPhi_f            = {fReader, "j2metPhi_f"};
        TTreeReaderValue<Double_t> drj1a_f               = {fReader, "drj1a_f"};
        TTreeReaderValue<Double_t> drj2a_f               = {fReader, "drj2a_f"};
        TTreeReaderValue<Double_t> drj1l_f               = {fReader, "drj1l_f"};
        TTreeReaderValue<Double_t> drj2l_f               = {fReader, "drj2l_f"};
        TTreeReaderValue<Double_t> jet1icsv_f            = {fReader, "jet1icsv_f"};
        TTreeReaderValue<Double_t> jet2icsv_f            = {fReader, "jet2icsv_f"};
        TTreeReaderValue<Double_t> drla_f                = {fReader, "drla"};

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
                && *photonhaspixelseed==0 && *photonet_f>25. && *photonet_f < 400.
                && *MET_et>30. && *mtVlepJECnew>30
                && *drla_f>0.5 && *drj1l_f>0.5 && *drj2l_f>0.5 && *drj1a_f>0.5 && *drj2a_f>0.5
                && abs(*j1metPhi_f)>0.5 && abs(*j2metPhi_f)>0.5)) continue;
	
		//cout<<i<<" pass"<<endl;
		etal1 = *etalep1;
		if (*ptlep1 >= 50) ptl1 = 45;
		else ptl1 = *ptlep1;
		
		//barell hist fill
		if (fabs(*photonsceta_f) < 1.4442){

			barrel_photonet->Fill(*photonet_f);
		}

                //endcap hist fill
                if (fabs(*photonsceta_f) > 1.566 && fabs(*photonsceta_f) < 2.5){
			endcap_photonet->Fill(*photonet_f);
                }

	}

        output->cd();
       	barrel_photonet->Write();
        endcap_photonet->Write();
        //output->Write();
        output->Close();
}
