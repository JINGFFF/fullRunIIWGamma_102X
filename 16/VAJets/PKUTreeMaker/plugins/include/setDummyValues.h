#ifndef _Setdummy_
#define _Setdummy_

void PKUTreeMaker::setDummyValues() {
     hasphoton=0.;

    npT                 = -1e1;
    npIT                = -1e1;
    nBX                 = -1e1;
    nevent              = -1e1;
    run                 = -1e1;
    ls                  = -1e1;
    nVtx                = -1e1;
    triggerWeight       = -1e1;
    pileupWeight        = -1e1;
    lumiWeight          = -1e1;
    theWeight           = -99;
    lep                 = -1e1;
    nlooseeles          = -1e1;
    nloosemus           = -1e1;
    ngoodmus 			= -1e1;
    ngoodeles 			= -1e1;
 
   	ptVlep              = -1e1;
    yVlep               = -1e1;
    phiVlep             = -1e1;
    massVlep            = -1e1;
    mtVlep              = -1e1;
    mtVlepnew           = -1e1;
    ptVlepJEC           = -1e1;
    yVlepJEC            = -1e1;
    phiVlepJEC          = -1e1;
    massVlepJEC         = -1e1;
    mtVlepJEC           = -1e1;
    mtVlepJECnew        = -1e1;
    ptVlepJEC_new       = -1e1;
    yVlepJEC_new        = -1e1;
    phiVlepJEC_new      = -1e1;
    massVlepJEC_new     = -1e1;
    mtVlepJEC_new       = -1e1;
    mtVlepJECnew_new    = -1e1;
    ptVlepJEC_JEC_up    = -1e1;
    yVlepJEC_JEC_up     = -1e1;
    phiVlepJEC_JEC_up  	= -1e1;
    massVlepJEC_JEC_up  = -1e1;
    mtVlepJEC_JEC_up   	= -1e1;
    mtVlepJECnew_JEC_up = -1e1;
    ptVlepJEC_JEC_down  = -1e1;
    yVlepJEC_JEC_down   = -1e1;
    phiVlepJEC_JEC_down = -1e1;
    massVlepJEC_JEC_down= -1e1;
    mtVlepJEC_JEC_down  = -1e1;
    mtVlepJECnew_JEC_down	= -1e1;
    ptVlepJEC_JER_up    = -1e1;
    yVlepJEC_JER_up     = -1e1;
    phiVlepJEC_JER_up   = -1e1;
    massVlepJEC_JER_up  = -1e1;
    mtVlepJEC_JER_up    = -1e1;
    mtVlepJECnew_JER_up = -1e1;
    ptVlepJEC_JER_down  = -1e1;
    yVlepJEC_JER_down 	= -1e1;
    phiVlepJEC_JER_down = -1e1;
    massVlepJEC_JER_down= -1e1;
    mtVlepJEC_JER_down  = -1e1;
    mtVlepJECnew_JER_down	= -1e1;
    Mla                 = -1e1;
    Mla_f               = -1e1;
    Mva                 = -1e1;
    Mva_f               = -1e1;
    ptlep1              = -1e1;
    etalep1             = -1e1;
    philep1             = -1e1;
    energylep1          = -1e1;
    met                 = -1e1;
    metPhi              = -1e1;
    j1metPhi            = -1e1;
    j1metPhi_f          = -1e1;
    j2metPhi            = -1e1;
    j2metPhi_f          = -1e1;
    j1metPhi_new        = -1e1;
    j1metPhi_new_f      = -1e1;
    j1metPhi_JEC_up     = -1e1;
    j1metPhi_JEC_down   = -1e1;
    j1metPhi_JER_up     = -1e1;
    j1metPhi_JER_down   = -1e1;
    j1metPhi_JEC_up_f   = -1e1;
    j1metPhi_JEC_down_f = -1e1;
    j1metPhi_JER_up_f   = -1e1;
    j1metPhi_JER_down_f = -1e1;
    j2metPhi_new        = -1e1;
    j2metPhi_new_f      = -1e1;
    j2metPhi_JEC_up     = -1e1;
    j2metPhi_JEC_down   = -1e1;
    j2metPhi_JER_up     = -1e1;
    j2metPhi_JER_down   = -1e1;
    j2metPhi_JEC_up_f   = -1e1;
    j2metPhi_JEC_down_f = -1e1;
    j2metPhi_JER_up_f   = -1e1;
    j2metPhi_JER_down_f = -1e1;
    METraw_et           = -99;
    METraw_phi          = -99;
    METraw_sumEt        = -99;
    genMET              = -99;
    MET_et_new          = -99;
    MET_phi_new         = -99;
    MET_sumEt_new       = -99;
    MET_et              = -99;
    MET_phi             = -99;
    MET_sumEt           = -99;
    // Marked for debug
    MET_et_JEC_up      	= -99;
    MET_et_JEC_down    	= -99;
    MET_et_JER_up      	= -99;
    MET_et_JER_down    	= -99;
    MET_phi_JEC_up     	= -99;
    MET_phi_JEC_down   	= -99;
    MET_phi_JER_up     	= -99;
    MET_phi_JER_down   	= -99;
    MET_sumEt_JEC_up   	= -99;
    MET_sumEt_JEC_down 	= -99;
    MET_sumEt_JER_up   	= -99;
    MET_sumEt_JER_down 	= -99;
    MET_corrPx         	= -99;
    MET_corrPy         	= -99;

    Dphiwajj   			= -1e1;
    Dphiwajj_f 			= -1e1;
    Dphiwajj_new   		= -1e1;
    Dphiwajj_JEC_up   	= -1e1;
    Dphiwajj_JEC_down   = -1e1;
    Dphiwajj_JER_up   	= -1e1;
    Dphiwajj_JER_down   = -1e1;
    // Marked for debug

    for (int i = 0; i < 6; i++) {
        ak4jet_hf[i]       	= -1e1;
        ak4jet_pf[i]       	= -1e1;
        genphoton_pt[i]    	= -1e1;
        genphoton_eta[i]   	= -1e1;
        genphoton_phi[i]   	= -1e1;
        genmuon_pt[i]      	= -1e1;
        genmuon_eta[i]     	= -1e1;
        genmuon_phi[i]     	= -1e1;
        genelectron_pt[i]  	= -1e1;
        genelectron_eta[i] 	= -1e1;
        genelectron_phi[i] 	= -1e1;
        photon_pt[i]       	= -1e1;
        photon_eta[i]      	= -1e1;
        photon_phi[i]      	= -1e1;
        photon_e[i]        	= -1e1;
        photonsc_eta[i]    	= -1e1;
        photonsc_phi[i]    	= -1e1;
        photon_pev[i]      	= false;
        photon_pevnew[i]   	= false;
        photon_ppsv[i]     	= false;
        photon_iseb[i]     	= false;
        photon_isee[i]     	= false;
        photon_hoe[i]      	= -1e1;
        photon_sieie[i]    	= -1e1;
        photon_sieie2[i]   	= -1e1;
        photon_chiso[i]    	= -1e1;
        photon_nhiso[i]    	= -1e1;
        photon_phoiso[i]   	= -1e1;
        photon_istrue[i]   	= -1;
        photon_isprompt[i] 	= -1;
        photon_drla[i]     	= -1e1;
        photon_mla[i]      	= -1e1;
        photon_mva[i]      	= -1e1;
        ak4jet_pt[i]       	= -1e1;
        ak4jet_eta[i] 		= -1e1;
        ak4jet_phi[i] 		= -1e1;
        ak4jet_e[i]    		= -1e1;
        ak4jet_csv[i]  		= -1e1;
        ak4jet_icsv[i] 		= -1e1;
    }

    photonet             = -1e1;
    photonet_f           = -1e1;
    photoneta            = -1e1;
    photoneta_f          = -1e1;
    photonphi            = -1e1;
    photonphi_f          = -1e1;
    photone              = -1e1;
    photone_f            = -1e1;
    photonsceta          = -1e1;
    photonsceta_f        = -1e1;
    photonscphi          = -1e1;
    photonscphi_f        = -1e1;
    photonsieie          = -1e1;
    photonsieie_f        = -1e1;
    photonphoiso         = -1e1;
    photonphoiso_f       = -1e1;
    photonchiso          = -1e1;
    photonchiso_f        = -1e1;
    photonnhiso          = -1e1;
    photonnhiso_f        = -1e1;
    iphoton              = -1;
    iphoton_f            = -1;
    drla                 = -1e1;
    drla_f               = -1e1;
    passEleVeto          = false;
    passEleVetonew       = false;
    passPixelSeedVeto    = false;
    photonhaspixelseed   = false;
    photonhaspixelseed_f = false;
    photonpasseleveto    = false;
    photonpasseleveto_f  = false;

    ISRPho              = false;
    dR_                 = 999;
    isTrue_             = -1;
    isprompt_           = -1;


   	jet1hf_orig			=-1e1;
   	jet1pf_orig			=-1e1;
    jet2hf_orig			=-1e1;
    jet2pf_orig			=-1e1;
    jet1pt_orig			=-1e1;
    jet1eta_orig		=-1e1;
    jet1phi_orig		=-1e1;
    jet1e_orig			=-1e1;
    jet2pt_orig			=-1e1;
    jet2eta_orig		=-1e1;
    jet2phi_orig		=-1e1;
   	jet2e_orig			=-1e1;
    jet1csv_orig 		=-1e1;
    jet2csv_orig 		=-1e1;
    jet1icsv_orig 		=-1e1;
    jet2icsv_orig 		=-1e1;
    drj1a_orig			=-1e1;
    drj2a_orig			=-1e1;
    drj1l_orig			=-1e1;
    drj2l_orig			=-1e1;

    jet1pt              = -1e1;
    jet1pt_f            = -1e1;
    jet1pt_new          = -1e1;
    jet1pt_new_f        = -1e1;
    jet1pt_JEC_up       = -1e1;
    jet1pt_JEC_down     = -1e1;
    jet1pt_JER_up       = -1e1;
    jet1pt_JER_down     = -1e1;
    jet1pt_JEC_up_f     = -1e1;
    jet1pt_JEC_down_f   = -1e1;
    jet1pt_JER_up_f     = -1e1;
    jet1pt_JER_down_f   = -1e1;
    jet1eta             = -1e1;
    jet1eta_f           = -1e1;
    jet1eta_new         = -1e1;
    jet1eta_new_f       = -1e1;
    jet1eta_JEC_up      = -1e1;
    jet1eta_JEC_down    = -1e1;
    jet1eta_JER_up      = -1e1;
    jet1eta_JER_down    = -1e1;
    jet1eta_JEC_up_f    = -1e1;
    jet1eta_JEC_down_f  = -1e1;
    jet1eta_JER_up_f    = -1e1;
    jet1eta_JER_down_f  = -1e1;
    jet1phi             = -1e1;
    jet1phi_f           = -1e1;
    jet1phi_new         = -1e1;
    jet1phi_new_f       = -1e1;
    jet1phi_JEC_up      = -1e1;
    jet1phi_JEC_down    = -1e1;
    jet1phi_JER_up      = -1e1;
    jet1phi_JER_down    = -1e1;
    jet1phi_JEC_up_f    = -1e1;
    jet1phi_JEC_down_f  = -1e1;
    jet1phi_JER_up_f    = -1e1;
    jet1phi_JER_down_f  = -1e1;
    jet1e               = -1e1;
    jet1e_f             = -1e1;
    jet1e_new           = -1e1;
    jet1e_new_f         = -1e1;
    jet1e_JEC_up        = -1e1;
    jet1e_JEC_down      = -1e1;
    jet1e_JER_up        = -1e1;
    jet1e_JER_down      = -1e1;
    jet1e_JEC_up_f      = -1e1;
    jet1e_JEC_down_f    = -1e1;
    jet1e_JER_up_f      = -1e1;
    jet1e_JER_down_f    = -1e1;
    jet1csv             = -1e1;
    jet1csv_f           = -1e1;
    jet1csv_new         = -1e1;
    jet1csv_new_f       = -1e1;
    jet1csv_JEC_up      = -1e1;
    jet1csv_JEC_down    = -1e1;
    jet1csv_JER_up      = -1e1;
    jet1csv_JER_down    = -1e1;
    jet1csv_JEC_up_f    = -1e1;
    jet1csv_JEC_down_f  = -1e1;
    jet1csv_JER_up_f    = -1e1;
    jet1csv_JER_down_f  = -1e1;
    jet1icsv            = -1e1;
    jet1icsv_f          = -1e1;
    jet1icsv_new        = -1e1;
    jet1icsv_new_f      = -1e1;
    jet1icsv_JEC_up     = -1e1;
    jet1icsv_JEC_down   = -1e1;
    jet1icsv_JER_up     = -1e1;
    jet1icsv_JER_down   = -1e1;
    jet1icsv_JEC_up_f   = -1e1;
    jet1icsv_JEC_down_f = -1e1;
    jet1icsv_JER_up_f   = -1e1;
    jet1icsv_JER_down_f = -1e1;
    jet2pt              = -1e1;
    jet2pt_f            = -1e1;
    jet2pt_new          = -1e1;
    jet2pt_new_f        = -1e1;
    jet2pt_JEC_up       = -1e1;
    jet2pt_JEC_down     = -1e1;
    jet2pt_JER_up       = -1e1;
    jet2pt_JER_down     = -1e1;
    jet2pt_JEC_up_f     = -1e1;
    jet2pt_JEC_down_f   = -1e1;
    jet2pt_JER_up_f     = -1e1;
    jet2pt_JER_down_f   = -1e1;
    jet2eta             = -1e1;
    jet2eta_f           = -1e1;
    jet2eta_new         = -1e1;
    jet2eta_new_f       = -1e1;
    jet2eta_JEC_up      = -1e1;
    jet2eta_JEC_down    = -1e1;
    jet2eta_JER_up      = -1e1;
    jet2eta_JER_down    = -1e1;
    jet2eta_JEC_up_f    = -1e1;
    jet2eta_JEC_down_f  = -1e1;
    jet2eta_JER_up_f    = -1e1;
    jet2eta_JER_down_f  = -1e1;
    jet2phi             = -1e1;
    jet2phi_f           = -1e1;
    jet2phi_new         = -1e1;
    jet2phi_new_f       = -1e1;
    jet2phi_JEC_up      = -1e1;
    jet2phi_JEC_down    = -1e1;
    jet2phi_JER_up      = -1e1;
    jet2phi_JER_down    = -1e1;
    jet2phi_JEC_up_f    = -1e1;
    jet2phi_JEC_down_f  = -1e1;
    jet2phi_JER_up_f    = -1e1;
    jet2phi_JER_down_f  = -1e1;
    jet2e               = -1e1;
    jet2e_f             = -1e1;
    jet2e_new           = -1e1;
    jet2e_new_f         = -1e1;
    jet2e_JEC_up        = -1e1;
    jet2e_JEC_down      = -1e1;
    jet2e_JER_up        = -1e1;
    jet2e_JER_down      = -1e1;
    jet2e_JEC_up_f      = -1e1;
    jet2e_JEC_down_f    = -1e1;
    jet2e_JER_up_f      = -1e1;
    jet2e_JER_down_f    = -1e1;
    jet2csv             = -1e1;
    jet2csv_f           = -1e1;
    jet2csv_new         = -1e1;
    jet2csv_new_f       = -1e1;
    jet2csv_JEC_up      = -1e1;
    jet2csv_JEC_down    = -1e1;
    jet2csv_JER_up      = -1e1;
    jet2csv_JER_down    = -1e1;
    jet2csv_JEC_up_f    = -1e1;
    jet2csv_JEC_down_f  = -1e1;
    jet2csv_JER_up_f    = -1e1;
    jet2csv_JER_down_f  = -1e1;
    jet2icsv            = -1e1;
    jet2icsv_f          = -1e1;
    jet2icsv_new        = -1e1;
    jet2icsv_new_f      = -1e1;
    jet2icsv_JEC_up     = -1e1;
    jet2icsv_JEC_down   = -1e1;
    jet2icsv_JER_up     = -1e1;
    jet2icsv_JER_down   = -1e1;
    jet2icsv_JEC_up_f   = -1e1;
    jet2icsv_JEC_down_f = -1e1;
    jet2icsv_JER_up_f   = -1e1;
    jet2icsv_JER_down_f = -1e1;
    drj1a               = -1e1;
    drj1a_f             = -1e1;
    drj1a_new           = -1e1;
    drj1a_new_f         = -1e1;
    drj1a_JEC_up        = -1e1;
    drj1a_JEC_down      = -1e1;
    drj1a_JER_up        = -1e1;
    drj1a_JER_down      = -1e1;
    drj1a_JEC_up_f      = -1e1;
    drj1a_JEC_down_f    = -1e1;
    drj1a_JER_up_f      = -1e1;
    drj1a_JER_down_f    = -1e1;
    drj2a               = -1e1;
    drj2a_f             = -1e1;
    drj2a_new           = -1e1;
    drj2a_new_f         = -1e1;
    drj2a_JEC_up        = -1e1;
    drj2a_JEC_down      = -1e1;
    drj2a_JER_up        = -1e1;
    drj2a_JER_down      = -1e1;
    drj2a_JEC_up_f      = -1e1;
    drj2a_JEC_down_f    = -1e1;
    drj2a_JER_up_f      = -1e1;
    drj2a_JER_down_f    = -1e1;
    drj1l               = -1e1;
    drj1l_f             = -1e1;
    drj1l_new           = -1e1;
    drj1l_new_f         = -1e1;
    drj1l_JEC_up        = -1e1;
    drj1l_JEC_down      = -1e1;
    drj1l_JER_up        = -1e1;
    drj1l_JER_down      = -1e1;
    drj1l_JEC_up_f      = -1e1;
    drj1l_JEC_down_f    = -1e1;
    drj1l_JER_up_f      = -1e1;
    drj1l_JER_down_f    = -1e1;
    drj2l               = -1e1;
    drj2l_f             = -1e1;
    drj2l_new           = -1e1;
    drj2l_new_f         = -1e1;
    drj2l_JEC_up        = -1e1;
    drj2l_JEC_down      = -1e1;
    drj2l_JER_up        = -1e1;
    drj2l_JER_down      = -1e1;
    drj2l_JEC_up_f      = -1e1;
    drj2l_JEC_down_f    = -1e1;
    drj2l_JER_up_f      = -1e1;
    drj2l_JER_down_f    = -1e1;
    Mjj                 = -1e1;
    Mjj_f               = -1e1;
    Mjj_new             = -1e1;
    Mjj_new_f           = -1e1;
    Mjj_JEC_up          = -1e1;
    Mjj_JEC_down        = -1e1;
    Mjj_JER_up          = -1e1;
    Mjj_JER_down        = -1e1;
    Mjj_JEC_up_f        = -1e1;
    Mjj_JEC_down_f      = -1e1;
    Mjj_JER_up_f        = -1e1;
    Mjj_JER_down_f      = -1e1;
    deltaeta            = -1e1;
    deltaeta_f          = -1e1;
    deltaeta_new        = -1e1;
    deltaeta_new_f      = -1e1;
    deltaeta_JEC_up     = -1e1;
    deltaeta_JEC_down   = -1e1;
    deltaeta_JER_up     = -1e1;
    deltaeta_JER_down   = -1e1;
    deltaeta_JEC_up_f   = -1e1;
    deltaeta_JEC_down_f = -1e1;
    deltaeta_JER_up_f   = -1e1;
    deltaeta_JER_down_f = -1e1;
    zepp                = -1e1;
    zepp_f              = -1e1;
    zepp_new            = -1e1;
    zepp_new_f          = -1e1;
    zepp_JEC_up         = -1e1;
    zepp_JEC_down       = -1e1;
    zepp_JER_up         = -1e1;
    zepp_JER_down       = -1e1;
    zepp_JEC_up_f       = -1e1;
    zepp_JEC_down_f     = -1e1;
    zepp_JER_up_f       = -1e1;
    zepp_JER_down_f     = -1e1;

    ISRPho        		= false;
    dR1_          		= 999;
    ispromptLep_  		= -1;
    lepton_istrue 		= -1;

    _prefiringweight	=-10;
    _prefiringweightup	=-10;
    _prefiringweightdown=-10;

    HLT_Ele1 			= -99;
    HLT_Ele2 			= -99;
    HLT_Mu1  			= -99;
    HLT_Mu2  			= -99;
    HLT_Mu3  			= -99;

    jet1hf   			= -1e1;
    jet1pf   			= -1e1;
    jet2hf   			= -1e1;
    jet2pf   			= -1e1;
    jet1hf_f 			= -1e1;
    jet1pf_f 			= -1e1;
    jet2hf_f 			= -1e1;
    jet2pf_f 			= -1e1;

    passFilter_HBHE_    = false;
    passFilter_HBHEIso_ = false;
    passFilter_globalTightHalo_	= false;
    passFilter_ECALDeadCell_    = false;
    passFilter_GoodVtx_       	= false;
    passFilter_EEBadSc_       	= false;
    passFilter_badMuon_       	= false;
    passFilter_badChargedHadron_= false;
}

/*
void PKUTreeMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {

    elPaths1.clear();
    elPaths2.clear();
    muPaths1.clear();
    muPaths2.clear();
    muPaths3.clear();

    std::cout << "-----begin-----" << std::endl;
    bool changed;
    if (!hltConfig.init(iRun, iSetup, "HLT", changed)) {
        edm::LogError("HltAnalysis") << "Initialization of HLTConfigProvider failed!!";
        return;
    }
    for (size_t i = 0; i < elPaths1_.size(); i++) {
        std::vector<std::string> foundPaths1 = hltConfig.matched(hltConfig.triggerNames(), elPaths1_[i]);
        while (!foundPaths1.empty()) {
            elPaths1.push_back(foundPaths1.back());
            foundPaths1.pop_back();
        }
    }
    for (size_t i = 0; i < muPaths1_.size(); i++) {
        std::vector<std::string> foundPaths1 = hltConfig.matched(hltConfig.triggerNames(), muPaths1_[i]);
        while (!foundPaths1.empty()) {
            muPaths1.push_back(foundPaths1.back());
            foundPaths1.pop_back();
        }
    }
    std::cout << "\n************** HLT-1 Information **************\n";
    for (size_t i = 0; i < elPaths1.size(); i++)
        std::cout << "\n Electron paths-1:    " << i << "  " << elPaths1[i].c_str() << "\t" << std::endl;
    for (size_t i = 0; i < muPaths1.size(); i++)
        std::cout << "\n Muon paths-1:   " << i << "  " << muPaths1[i].c_str() << "\t" << std::endl;
    std::cout << "\n*********************************************\n\n";

    for (size_t i = 0; i < elPaths2_.size(); i++) {
        std::vector<std::string> foundPaths2 = hltConfig.matched(hltConfig.triggerNames(), elPaths2_[i]);
        while (!foundPaths2.empty()) {
            elPaths2.push_back(foundPaths2.back());
            foundPaths2.pop_back();
        }
    }
    for (size_t i = 0; i < muPaths2_.size(); i++) {
        std::vector<std::string> foundPaths2 = hltConfig.matched(hltConfig.triggerNames(), muPaths2_[i]);
        while (!foundPaths2.empty()) {
            muPaths2.push_back(foundPaths2.back());
            foundPaths2.pop_back();
        }
    }
    std::cout << "\n************** HLT-2 Information **************\n";
    for (size_t i = 0; i < elPaths2.size(); i++)
        std::cout << "\n Electron paths-2:    " << i << "  " << elPaths2[i].c_str() << "\t" << std::endl;
    for (size_t i = 0; i < muPaths2.size(); i++)
        std::cout << "\n Muon paths-2:   " << i << "  " << muPaths2[i].c_str() << "\t" << std::endl;
    std::cout << "\n*********************************************\n\n";
    for (size_t i = 0; i < muPaths3_.size(); i++) {
        std::vector<std::string> foundPaths3 = hltConfig.matched(hltConfig.triggerNames(), muPaths3_[i]);
        while (!foundPaths3.empty()) {
            muPaths3.push_back(foundPaths3.back());
            foundPaths3.pop_back();
        }
    }

    std::cout << "\n************** HLT-3 Information **************\n";
    for (size_t i = 0; i < muPaths3.size(); i++)
        std::cout << "\n Muon paths-3:   " << i << "  " << muPaths3[i].c_str() << "\t" << std::endl;
    std::cout << "\n*********************************************\n\n";
}

void PKUTreeMaker::endJob() {
    std::cout << "PKUTreeMaker endJob()..." << std::endl;
}

*/
#endif
