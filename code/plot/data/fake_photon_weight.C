using namespace std;
void fake_photon_weight(){

        double fpf_barrel[7] = {0.332,     0.238,  0.165,  0.132,  0.0958, 0.0532, 0.0256};
        double fpf_endcap[7] = {0.0842,    0.0682, 0.0605, 0.0551, 0.036,  0.0203, 0.0424};
        //double fpf_elec_barrel[7] = {0.347,     0.225,  0.172,  0.124,  0.0814, 0.0502, 0.0255};
        //double fpf_elec_endcap[7] = {0.376,     0.297,  0.215,  0.176,  0.136,  0.101,  0.0};

        TFile* output = new TFile("fake_photon_weight.root","RECREATE");

	TFile * input1 = TFile::Open("hist_fakephoton_photonet_mu.root");
	TFile * input2 = TFile::Open("hist_photon_mu.root");
        TH1F* b = (TH1F*)input2->Get("barrel_photonet");
        TH1F* e = (TH1F*)input2->Get("endcap_photonet");

        TH1F* fb = (TH1F*)input1->Get("barrel_photonet");
        TH1F* fe = (TH1F*)input1->Get("barrel_photonet");

        TH1F* muon_barrel_fake_photon_weight = (TH1F*)b->Clone();
	muon_barrel_fake_photon_weight->SetTitle("muon_barrel_fake_photon_weight");
        muon_barrel_fake_photon_weight->SetName("muon_barrel_fake_photon_weight");

        TH1F* muon_endcap_fake_photon_weight = (TH1F*)e->Clone();
        muon_endcap_fake_photon_weight->SetTitle("muon_endcap_fake_photon_weight");
        muon_endcap_fake_photon_weight->SetName("muon_endcap_fake_photon_weight");

	muon_barrel_fake_photon_weight->Divide(fb);
        muon_endcap_fake_photon_weight->Divide(fe);

	int i;
	for(i=1; i<=7; i++){
		muon_barrel_fake_photon_weight->SetBinContent(i,fpf_barrel[i-1]*muon_barrel_fake_photon_weight->GetBinContent(i));
                muon_endcap_fake_photon_weight->SetBinContent(i,fpf_endcap[i-1]*muon_endcap_fake_photon_weight->GetBinContent(i));


	}


	output->cd();
	muon_barrel_fake_photon_weight->Write();
        muon_endcap_fake_photon_weight->Write();

	output->Close();
}
