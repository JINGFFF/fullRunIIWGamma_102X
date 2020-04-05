#ifndef _PKUTreeMakerHead_
#define _PKUTreeMakerHead_

using namespace std;
using namespace edm;

struct sortPt {
    bool operator()(TLorentzVector* s1, TLorentzVector* s2) const {
        return s1->Pt() >= s2->Pt();
    }
} mysortPt;
//
// class declaration
//

class PKUTreeMaker : public edm::EDAnalyzer {
public:
    explicit PKUTreeMaker(const edm::ParameterSet&);
    ~PKUTreeMaker();
    //static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

    enum PhotonMatchType { UNMATCHED = 0,
                           MATCHED_FROM_GUDSCB,
                           MATCHED_FROM_PI0,
                           MATCHED_FROM_OTHER_SOURCES };

private:
    //virtual void   beginJob() override;
    virtual void   analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void   endJob() override;
    virtual void   beginRun(const edm::Run&, const edm::EventSetup&) override;
    //virtual void   endRun(const edm::Run&, const edm::EventSetup&) override;
    virtual void   addTypeICorr(edm::Event const& event);
    virtual void   addTypeICorr_user(edm::Event const& event);  //---for MET, Meng
    virtual double getJEC(reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_);
    virtual double getJECOffset(reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_);

    math::XYZTLorentzVector getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType);

    bool hasMatchedPromptElectron(const reco::SuperClusterRef& sc, const edm::Handle<edm::View<pat::Electron>>& eleCol, const edm::Handle<reco::ConversionCollection>& convCol, const math::XYZPoint& beamspot, float lxyMin = 2.0, float probMin = 1e-6, unsigned int nHitsBeforeVtxMax = 0);
    int  matchToTruth(const reco::Photon& pho, const edm::Handle<edm::View<reco::GenParticle>>& genParticles, bool& ISRPho, double& dR, int& isprompt);

    int matchToTrueLep(double lept_eta, double lept_phi, const edm::Handle<edm::View<reco::GenParticle>>& genParticles, double& dR, int& ispromptLep);  //////////////////////////////////

    void findFirstNonPhotonMother(const reco::Candidate* particle, int& ancestorPID, int& ancestorStatus);
    void setDummyValues();

    float                                             EAch(float x);
    float                                             EAnh(float x);
    float                                             EApho(float x);
    std::vector<std::string>                          offsetCorrLabel_;
    FactorizedJetCorrector*                           jecOffset_;
    std::vector<std::string>                          jetCorrLabel_;

    edm::EDGetTokenT<edm::View<pat::Muon>> goodmuonToken_;
    edm::EDGetTokenT<edm::View<pat::Electron>> goodeleToken_;
    edm::Handle<double>                               rho_;
    edm::EDGetTokenT<double>                          rhoToken_;
    edm::EDGetTokenT<pat::METCollection>              metInputToken_;
    std::vector<edm::EDGetTokenT<pat::METCollection>> mettokens;
    edm::EDGetTokenT<pat::METCollection>              metToken_;
    edm::EDGetTokenT<edm::View<pat::Electron>>        electronToken_;
    edm::EDGetTokenT<edm::View<pat::Photon>>          photonToken_;
    edm::EDGetTokenT<reco::BeamSpot>                  beamSpotToken_;
    edm::EDGetTokenT<std::vector<reco::Conversion>>   conversionsToken_;
    edm::EDGetTokenT<edm::View<pat::Electron>>        looseelectronToken_;
    edm::EDGetTokenT<edm::View<pat::Muon>>            loosemuonToken_;

    //L1 prefiring
    edm::EDGetTokenT< double > prefweight_token;
    edm::EDGetTokenT< double > prefweightup_token;
    edm::EDGetTokenT< double > prefweightdown_token;


    // Filter
    edm::EDGetTokenT<edm::TriggerResults> noiseFilterToken_;
    edm::Handle<edm::TriggerResults>      noiseFilterBits_;
    std::string                           HBHENoiseFilter_Selector_;
    std::string                           HBHENoiseIsoFilter_Selector_;
    std::string                           ECALDeadCellNoiseFilter_Selector_;
    std::string                           GoodVtxNoiseFilter_Selector_;
    std::string                           EEBadScNoiseFilter_Selector_;
    std::string                           globalTightHaloFilter_Selector_;
    edm::EDGetTokenT<bool>                badMuon_Selector_;
    edm::EDGetTokenT<bool>                badChargedHadron_Selector_;

    edm::EDGetTokenT<edm::ValueMap<float>> full5x5SigmaIEtaIEtaMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> phoChargedIsolationToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> phoNeutralHadronIsolationToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> phoPhotonIsolationToken_;
    EffectiveAreas                         effAreaChHadrons_;
    EffectiveAreas                         effAreaNeuHadrons_;
    EffectiveAreas                         effAreaPhotons_;

    // ----------member data ---------------------------
    TTree* outTree_;
    double hasphoton;
    double MW_;
    int    nevent, run, ls;
    int    nVtx;
    double triggerWeight, lumiWeight, pileupWeight;
    double theWeight;
    double nump = 0.;
    double numm = 0.;
    double npT, npIT;
    int    nBX;
    double ptVlep, yVlep, phiVlep, massVlep, mtVlep, mtVlepnew;
    double ptVlepJEC, yVlepJEC, phiVlepJEC, massVlepJEC, mtVlepJEC, mtVlepJECnew;
    double ptVlepJEC_new, yVlepJEC_new, phiVlepJEC_new, massVlepJEC_new, mtVlepJEC_new, mtVlepJECnew_new;
    double ptVlepJEC_JEC_up, yVlepJEC_JEC_up, phiVlepJEC_JEC_up, massVlepJEC_JEC_up, mtVlepJEC_JEC_up, mtVlepJECnew_JEC_up;
    double ptVlepJEC_JEC_down, yVlepJEC_JEC_down, phiVlepJEC_JEC_down, massVlepJEC_JEC_down, mtVlepJEC_JEC_down, mtVlepJECnew_JEC_down;
    double ptVlepJEC_JER_up, yVlepJEC_JER_up, phiVlepJEC_JER_up, massVlepJEC_JER_up, mtVlepJEC_JER_up, mtVlepJECnew_JER_up;
    double ptVlepJEC_JER_down, yVlepJEC_JER_down, phiVlepJEC_JER_down, massVlepJEC_JER_down, mtVlepJEC_JER_down, mtVlepJECnew_JER_down;
    double Mla, Mva;
    double Mla_f, Mva_f;
    double ptlep1, etalep1, philep1, energylep1;
    int    lep, nlooseeles, nloosemus, ngoodeles, ngoodmus;
    double _prefiringweight,_prefiringweightup,_prefiringweightdown;
    double met, metPhi, j1metPhi, j2metPhi;
    double j1metPhi_new, j1metPhi_JEC_up, j1metPhi_JEC_down, j1metPhi_JER_up, j1metPhi_JER_down;
    double j2metPhi_new, j2metPhi_JEC_up, j2metPhi_JEC_down, j2metPhi_JER_up, j2metPhi_JER_down;
    double j1metPhi_f, j2metPhi_f;
    double j1metPhi_new_f, j1metPhi_JEC_up_f, j1metPhi_JEC_down_f, j1metPhi_JER_up_f, j1metPhi_JER_down_f;
    double j2metPhi_new_f, j2metPhi_JEC_up_f, j2metPhi_JEC_down_f, j2metPhi_JER_up_f, j2metPhi_JER_down_f;
    //Met JEC
    float  rawPt;
    double METraw_et, METraw_phi, METraw_sumEt;
    double genMET, MET_et, MET_phi, MET_sumEt, MET_corrPx, MET_corrPy;
    double MET_et_new, MET_phi_new, MET_sumEt_new;
    // Marked for debug
    //-------------- Met uncertainty ----------------//
    double MET_et_JEC_up, MET_et_JEC_down, MET_et_JER_up, MET_et_JER_down;
    double MET_phi_JEC_up, MET_phi_JEC_down, MET_phi_JER_up, MET_phi_JER_down;
    double MET_sumEt_JEC_up, MET_sumEt_JEC_down, MET_sumEt_JER_up, MET_sumEt_JER_down;
    //-------------- Met uncertainty-----------------//
    double useless;
    // AK4 Jets
    int    ak4jet_hf[6], ak4jet_pf[6];
    double ak4jet_pt[6], ak4jet_eta[6], ak4jet_phi[6], ak4jet_e[6];
    double ak4jet_pt_old[6], ak4jet_e_old[6];
    double ak4jet_pt_new[6], ak4jet_e_new[6];
    double ak4jet_pt_JEC_up[6], ak4jet_pt_JEC_down[6], ak4jet_e_JEC_up[6], ak4jet_e_JEC_down[6];
    double ak4jet_pt_JER_up[6], ak4jet_pt_JER_down[6], ak4jet_e_JER_up[6], ak4jet_e_JER_down[6];
    double ak4jet_pt_jer[6];
    double ak4jet_csv[6], ak4jet_icsv[6];
    double drjetlep[6], drjetphoton[6];
    double genjet_pt[6], genjet_eta[6], genjet_phi[6], genjet_e[6];

    //Photon
    double genphoton_pt[6], genphoton_eta[6], genphoton_phi[6];
    double genmuon_pt[6], genmuon_eta[6], genmuon_phi[6];
    double genelectron_pt[6], genelectron_eta[6], genelectron_phi[6];
    double photon_pt[6], photon_eta[6], photon_phi[6], photon_e[6], photonsc_eta[6], photonsc_phi[6];
    bool   photon_pev[6], photon_pevnew[6], photon_ppsv[6], photon_iseb[6], photon_isee[6];
    double photon_hoe[6], photon_sieie[6], photon_sieie2[6], photon_chiso[6], photon_nhiso[6], photon_phoiso[6], photon_drla[6], photon_mla[6], photon_mva[6];
    int    photon_istrue[6], photon_isprompt[6];
    double photonet, photoneta, photonphi, photone, photonsceta, photonscphi;
    double photonet_f, photoneta_f, photonphi_f, photone_f, photonsceta_f, photonscphi_f;
    double photonsieie, photonphoiso, photonchiso, photonnhiso;
    double photonsieie_f, photonphoiso_f, photonchiso_f, photonnhiso_f;
    int    iphoton;
    int    iphoton_f;
    double drla, drla_f;
    bool   passEleVeto, passEleVetonew, passPixelSeedVeto, photonhaspixelseed, photonhaspixelseed_f, photonpasseleveto, photonpasseleveto_f;
    //Photon gen match
    int    isTrue_;
    bool   ISRPho;
    int    isprompt_, ispromptLep_, lepton_istrue;  //////////////////////////////////
    double dR_, dR1_;

    int jet1hf_orig,jet1pf_orig,jet2hf_orig,jet2pf_orig;
    double jet1pt_orig, jet1eta_orig, jet1phi_orig, jet1e_orig, jet1csv_orig, jet1icsv_orig;
    double jet2pt_orig, jet2eta_orig, jet2phi_orig, jet2e_orig, jet2csv_orig, jet2icsv_orig;
    double drj1a_orig, drj2a_orig, drj1l_orig, drj2l_orig;


    //Jets
    int    jet1hf, jet1pf, jet2hf, jet2pf, jet1hf_f, jet1pf_f, jet2hf_f, jet2pf_f;
    double Dphiwajj, Dphiwajj_f,Dphiwajj_new,Dphiwajj_JEC_up,Dphiwajj_JEC_down,Dphiwajj_JER_up,Dphiwajj_JER_down;

    double jet1pt, jet1eta, jet1phi, jet1e, jet1csv, jet1icsv;
    double jet1pt_new, jet1pt_JEC_up, jet1pt_JEC_down, jet1pt_JER_up, jet1pt_JER_down;
    double jet1e_new, jet1e_JEC_up, jet1e_JEC_down, jet1e_JER_up, jet1e_JER_down;
    double jet1eta_new, jet1eta_JEC_up, jet1eta_JEC_down, jet1eta_JER_up, jet1eta_JER_down;
    double jet1phi_new, jet1phi_JEC_up, jet1phi_JEC_down, jet1phi_JER_up, jet1phi_JER_down;
    double jet1csv_new, jet1csv_JEC_up, jet1csv_JEC_down, jet1csv_JER_up, jet1csv_JER_down;
    double jet1icsv_new, jet1icsv_JEC_up, jet1icsv_JEC_down, jet1icsv_JER_up, jet1icsv_JER_down;

    double jet1pt_f, jet1eta_f, jet1phi_f, jet1e_f, jet1csv_f, jet1icsv_f;
    double jet1pt_new_f, jet1pt_JEC_up_f, jet1pt_JEC_down_f, jet1pt_JER_up_f, jet1pt_JER_down_f;
    double jet1e_new_f, jet1e_JEC_up_f, jet1e_JEC_down_f, jet1e_JER_up_f, jet1e_JER_down_f;
    double jet1eta_new_f, jet1eta_JEC_up_f, jet1eta_JEC_down_f, jet1eta_JER_up_f, jet1eta_JER_down_f;
    double jet1phi_new_f, jet1phi_JEC_up_f, jet1phi_JEC_down_f, jet1phi_JER_up_f, jet1phi_JER_down_f;
    double jet1csv_new_f, jet1csv_JEC_up_f, jet1csv_JEC_down_f, jet1csv_JER_up_f, jet1csv_JER_down_f;
    double jet1icsv_new_f, jet1icsv_JEC_up_f, jet1icsv_JEC_down_f, jet1icsv_JER_up_f, jet1icsv_JER_down_f;

    double jet2pt, jet2eta, jet2phi, jet2e, jet2csv, jet2icsv;
    double jet2pt_new, jet2pt_JEC_up, jet2pt_JEC_down, jet2pt_JER_up, jet2pt_JER_down;
    double jet2e_new, jet2e_JEC_up, jet2e_JEC_down, jet2e_JER_up, jet2e_JER_down;
    double jet2eta_new, jet2eta_JEC_up, jet2eta_JEC_down, jet2eta_JER_up, jet2eta_JER_down;
    double jet2phi_new, jet2phi_JEC_up, jet2phi_JEC_down, jet2phi_JER_up, jet2phi_JER_down;
    double jet2csv_new, jet2csv_JEC_up, jet2csv_JEC_down, jet2csv_JER_up, jet2csv_JER_down;
    double jet2icsv_new, jet2icsv_JEC_up, jet2icsv_JEC_down, jet2icsv_JER_up, jet2icsv_JER_down;

    double jet2pt_f, jet2eta_f, jet2phi_f, jet2e_f, jet2csv_f, jet2icsv_f;
    double jet2pt_new_f, jet2e_new_f, jet2pt_JEC_up_f, jet2pt_JEC_down_f, jet2pt_JER_up_f, jet2pt_JER_down_f;
    double jet2e_JEC_up_f, jet2e_JEC_down_f, jet2e_JER_up_f, jet2e_JER_down_f;
    double jet2eta_new_f, jet2eta_JEC_up_f, jet2eta_JEC_down_f, jet2eta_JER_up_f, jet2eta_JER_down_f;
    double jet2phi_new_f, jet2phi_JEC_up_f, jet2phi_JEC_down_f, jet2phi_JER_up_f, jet2phi_JER_down_f;
    double jet2csv_new_f, jet2csv_JEC_up_f, jet2csv_JEC_down_f, jet2csv_JER_up_f, jet2csv_JER_down_f;
    double jet2icsv_new_f, jet2icsv_JEC_up_f, jet2icsv_JEC_down_f, jet2icsv_JER_up_f, jet2icsv_JER_down_f;

    double drj1a, drj2a, drj1l, drj2l;
    double drj1a_f, drj2a_f, drj1l_f, drj2l_f;

    double drj1a_new, drj1a_JEC_up, drj1a_JEC_down, drj1a_JER_up, drj1a_JER_down;
    double drj2a_new, drj2a_JEC_up, drj2a_JEC_down, drj2a_JER_up, drj2a_JER_down;
    double drj1l_new, drj1l_JEC_up, drj1l_JEC_down, drj1l_JER_up, drj1l_JER_down;
    double drj2l_new, drj2l_JEC_up, drj2l_JEC_down, drj2l_JER_up, drj2l_JER_down;

    double drj1a_new_f, drj1a_JEC_up_f, drj1a_JEC_down_f, drj1a_JER_up_f, drj1a_JER_down_f;
    double drj2a_new_f, drj2a_JEC_up_f, drj2a_JEC_down_f, drj2a_JER_up_f, drj2a_JER_down_f;
    double drj1l_new_f, drj1l_JEC_up_f, drj1l_JEC_down_f, drj1l_JER_up_f, drj1l_JER_down_f;
    double drj2l_new_f, drj2l_JEC_up_f, drj2l_JEC_down_f, drj2l_JER_up_f, drj2l_JER_down_f;

    double Mjj, deltaeta, zepp;
    double deltaeta_new, deltaeta_JEC_up, deltaeta_JEC_down, deltaeta_JER_up, deltaeta_JER_down;
    double Mjj_new, Mjj_JEC_up, Mjj_JEC_down, Mjj_JER_up, Mjj_JER_down;
    double zepp_new, zepp_JEC_up, zepp_JEC_down, zepp_JER_up, zepp_JER_down;

    double Mjj_f, deltaeta_f, zepp_f;
    double deltaeta_new_f, deltaeta_JEC_up_f, deltaeta_JEC_down_f, deltaeta_JER_up_f, deltaeta_JER_down_f;
    double Mjj_new_f, Mjj_JEC_up_f, Mjj_JEC_down_f, Mjj_JER_up_f, Mjj_JER_down_f;
    double zepp_new_f, zepp_JEC_up_f, zepp_JEC_down_f, zepp_JER_up_f, zepp_JER_down_f;

	bool is_leptonicVs_Empty;
	double energyVlepJEC;

    /// Parameters to steer the treeDumper
    int                      originalNEvents_;
    double                   crossSectionPb_;
    double                   targetLumiInvPb_;
    std::string              PKUChannel_;
    bool                     isGen_, RunOnMC_;
    std::vector<std::string> jecAK4Labels_;
    std::vector<std::string> jecAK4chsLabels_;
    //correction jet
    FactorizedJetCorrector*       jecAK4_;
    std::string                   gravitonSrc_;
    std::map<std::string, double> TypeICorrMap_;
    std::map<std::string, double> TypeICorrMap_user_;
    edm::InputTag                 mets_;

    //High Level Trigger
    HLTConfigProvider                     hltConfig;
    edm::EDGetTokenT<edm::TriggerResults> hltToken_;
    std::vector<std::string>              elPaths1_, elPaths2_;
    std::vector<std::string>              muPaths1_, muPaths2_, muPaths3_;
    std::vector<std::string>              elPaths1, elPaths2;
    std::vector<std::string>              muPaths1, muPaths2, muPaths3;
    int                                   HLT_Ele1, HLT_Ele2;
    int                                   HLT_Mu1, HLT_Mu2, HLT_Mu3;

    // filter
    bool passFilter_HBHE_;
    bool passFilter_HBHEIso_;
    bool passFilter_globalTightHalo_;
    bool passFilter_ECALDeadCell_;
    bool passFilter_GoodVtx_;
    bool passFilter_EEBadSc_;
    bool passFilter_badMuon_;
    bool passFilter_badChargedHadron_;

    edm::EDGetTokenT<GenEventInfoProduct>            GenToken_;
    edm::EDGetTokenT<reco::GenJetCollection>         genJet_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUToken_;
    edm::EDGetTokenT<edm::View<reco::Candidate>>     leptonicVSrc_;
    edm::EDGetTokenT<edm::View<pat::Jet>>            ak4jetsSrc_;
    edm::EDGetTokenT<edm::View<pat::Photon>>         photonSrc_;
    edm::EDGetTokenT<edm::View<reco::GenParticle>>   genSrc_;
    edm::EDGetTokenT<edm::View<reco::Candidate>>     metSrc_;
    edm::EDGetTokenT<reco::VertexCollection>         VertexToken_;
    edm::EDGetTokenT<pat::JetCollection>             t1jetSrc_;
    edm::EDGetTokenT<pat::JetCollection>             t1jetSrc_user_;
    edm::EDGetTokenT<edm::View<pat::Muon>>           t1muSrc_;
};

float PKUTreeMaker::EAch( float x){
        float EA = 0.0112;
        if(x>1.0)   EA = 0.0108;
        if(x>1.479) EA = 0.0106;
        if(x>2.0)   EA = 0.01002;
        if(x>2.2)   EA = 0.0098;
        if(x>2.3)   EA = 0.0089;
        if(x>2.4)   EA = 0.0087;
        return EA;
}

float PKUTreeMaker::EAnh( float x){
        float EA = 0.0668;
        if(x>1.0)   EA = 0.1054;
        if(x>1.479) EA = 0.0786;
        if(x>2.0)   EA = 0.0233;
        if(x>2.2)   EA = 0.0078;
        if(x>2.3)   EA = 0.0028;
        if(x>2.4)   EA = 0.0137;
        return EA;
}

float PKUTreeMaker::EApho( float x){
        float EA = 0.1113;
        if(x>1.0)   EA = 0.0953;
        if(x>1.479) EA = 0.0619;
        if(x>2.0)   EA = 0.0837;
        if(x>2.2)   EA = 0.1070;
        if(x>2.3)   EA = 0.1212;
        if(x>2.4)   EA = 0.1466;
        return EA;
}

// constructors and destructor

PKUTreeMaker::PKUTreeMaker(const edm::ParameterSet& iConfig)  //:
    : effAreaChHadrons_((iConfig.getParameter<edm::FileInPath>("effAreaChHadFile")).fullPath()), 
	effAreaNeuHadrons_((iConfig.getParameter<edm::FileInPath>("effAreaNeuHadFile")).fullPath()), 
	effAreaPhotons_((iConfig.getParameter<edm::FileInPath>("effAreaPhoFile")).fullPath()) 
{
    hltToken_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltToken"));
    elPaths1_ = iConfig.getParameter<std::vector<std::string>>("elPaths1");
    elPaths2_ = iConfig.getParameter<std::vector<std::string>>("elPaths2");
    muPaths1_ = iConfig.getParameter<std::vector<std::string>>("muPaths1");
    muPaths2_ = iConfig.getParameter<std::vector<std::string>>("muPaths2");
    muPaths3_ = iConfig.getParameter<std::vector<std::string>>("muPaths3");
    GenToken_ = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"));
    genJet_   =consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJet"));
    //  LheToken_=consumes<LHEEventProduct> (iConfig.getParameter<edm::InputTag>( "lhe") ) ;
    PUToken_         = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileup"));
    leptonicVSrc_    = consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("leptonicVSrc"));
    ak4jetsSrc_      = consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("ak4jetsSrc"));
    photonSrc_       = consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photonSrc"));
    genSrc_          = consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genSrc"));
    metSrc_          = consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("metSrc"));
    VertexToken_     = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"));
    t1jetSrc_        = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("t1jetSrc"));
    t1jetSrc_user_   = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("t1jetSrc_user"));
    t1muSrc_         = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("t1muSrc"));
    originalNEvents_ = iConfig.getParameter<int>("originalNEvents");
    crossSectionPb_  = iConfig.getParameter<double>("crossSectionPb");
    targetLumiInvPb_ = iConfig.getParameter<double>("targetLumiInvPb");
    PKUChannel_      = iConfig.getParameter<std::string>("PKUChannel");
    isGen_           = iConfig.getParameter<bool>("isGen");
    RunOnMC_         = iConfig.getParameter<bool>("RunOnMC");
    rhoToken_        = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
    jecAK4chsLabels_ = iConfig.getParameter<std::vector<std::string>>("jecAK4chsPayloadNames");
    jecAK4Labels_    = iConfig.getParameter<std::vector<std::string>>("jecAK4PayloadNames");
    metToken_        = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metSrc"));
    mettokens.push_back(metToken_);
    metInputToken_      = mettokens[0];
    electronToken_      = (consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons")));
    looseelectronToken_ = (consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("looseelectronSrc")));
    loosemuonToken_     = (consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("loosemuonSrc")));

    goodmuonToken_      = (consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("goodmuonSrc")));
    goodeleToken_       = (consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("goodeleSrc")));

    beamSpotToken_      = (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot")));
    conversionsToken_   = (consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions")));

    //L1 prefiring
    prefweight_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"));
    prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
    prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"));
    
    jetCorrLabel_ = jecAK4chsLabels_;
    offsetCorrLabel_.push_back(jetCorrLabel_[0]);

    // filter
    noiseFilterToken_                 = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilter"));
    HBHENoiseFilter_Selector_         = iConfig.getParameter<std::string>("noiseFilterSelection_HBHENoiseFilter");
    HBHENoiseIsoFilter_Selector_      = iConfig.getParameter<std::string>("noiseFilterSelection_HBHENoiseIsoFilter");
    globalTightHaloFilter_Selector_   = iConfig.getParameter<std::string>("noiseFilterSelection_globalTightHaloFilter");
    ECALDeadCellNoiseFilter_Selector_ = iConfig.getParameter<std::string>("noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter");
    GoodVtxNoiseFilter_Selector_      = iConfig.getParameter<std::string>("noiseFilterSelection_goodVertices");
    EEBadScNoiseFilter_Selector_      = iConfig.getParameter<std::string>("noiseFilterSelection_eeBadScFilter");
    badMuon_Selector_                 = consumes<bool>(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_badMuon"));
    badChargedHadron_Selector_        = consumes<bool>(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_badChargedHadron"));

    full5x5SigmaIEtaIEtaMapToken_   = (consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap")));
    phoChargedIsolationToken_       = (consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("phoChargedIsolation")));
    phoNeutralHadronIsolationToken_ = (consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolation")));
    phoPhotonIsolationToken_        = (consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("phoPhotonIsolation")));

    MW_ = 80.385;
    //now do what ever initialization is needed
    edm::Service<TFileService> fs;
    outTree_ = fs->make<TTree>("PKUCandidates", "PKU Candidates");
    /// Basic event quantities

    outTree_->Branch("hasphoton", &hasphoton, "hasphoton/D");
    outTree_->Branch("ngoodmus", &ngoodmus, "ngoodmus/I");
    outTree_->Branch("ngoodeles", &ngoodeles, "ngoodeles/I");

    outTree_->Branch("jet1hf_orig",&jet1hf_orig,"jet1hf_orig/I");
    outTree_->Branch("jet1pf_orig",&jet1pf_orig,"jet1pf_orig/I");
    outTree_->Branch("jet2hf_orig",&jet2hf_orig,"jet2hf_orig/I");
    outTree_->Branch("jet2pf_orig",&jet2pf_orig,"jet2pf_orig/I");

    outTree_->Branch("jet1pt_orig",&jet1pt_orig,"jet1pt_orig/D");
    outTree_->Branch("jet1phi_orig",&jet1phi_orig,"jet1phi_orig/D");
    outTree_->Branch("jet1eta_orig",&jet1eta_orig,"jet1eta_orig/D");
    outTree_->Branch("jet1e_orig",&jet1e_orig,"jet1e_orig/D");
    outTree_->Branch("jet1csv_orig",&jet1csv_orig,"jet1csv_orig/D");
    outTree_->Branch("jet1icsv_orig",&jet1icsv_orig,"jet1icsv_orig/D");
    outTree_->Branch("drj1a_orig",&drj1a_orig,"drj1a_orig/D");
    outTree_->Branch("drj1l_orig",&drj1l_orig,"drj1l_orig/D");

    outTree_->Branch("jet2pt_orig",&jet2pt_orig,"jet2pt_orig/D");
    outTree_->Branch("jet2phi_orig",&jet2phi_orig,"jet2phi_orig/D");
    outTree_->Branch("jet2eta_orig",&jet2eta_orig,"jet2eta_orig/D");
    outTree_->Branch("jet2e_orig",&jet2e_orig,"jet2e_orig/D");
    outTree_->Branch("jet2csv_orig",&jet2csv_orig,"jet2csv_orig/D");
    outTree_->Branch("jet2icsv_orig",&jet2icsv_orig,"jet2icsv_orig/D");
    outTree_->Branch("drj2a_orig",&drj2a_orig,"drj2a_orig/D");
    outTree_->Branch("drj2l_orig",&drj2l_orig,"drj2l_orig/D");

    outTree_->Branch("rawPt",&rawPt,"rawPt/F");
    outTree_->Branch("nevent",&nevent,"nevent/I");
    outTree_->Branch("run",&run,"run/I");
    outTree_->Branch("ls",&ls,"ls/I");
    outTree_->Branch("nVtx", &nVtx, "nVtx/I");
    outTree_->Branch("theWeight", &theWeight, "theWeight/D");
    outTree_->Branch("nump", &nump, "nump/D");
    outTree_->Branch("numm", &numm, "numm/D");
    outTree_->Branch("npT", &npT, "npT/D");
    outTree_->Branch("lep", &lep, "lep/I");
    outTree_->Branch("ptVlep", &ptVlep, "ptVlep/D");
    outTree_->Branch("ptVlepJEC", &ptVlepJEC, "ptVlepJEC/D");
    outTree_->Branch("yVlepJEC", &yVlepJEC, "yVlepJEC/D");
    outTree_->Branch("phiVlepJEC", &phiVlepJEC, "phiVlepJEC/D");
    outTree_->Branch("massVlepJEC", &massVlepJEC, "massVlepJEC/D");
    outTree_->Branch("mtVlepJECnew", &mtVlepJECnew, "mtVlepJECnew/D");
    outTree_->Branch("ptVlepJEC_new", &ptVlepJEC_new, "ptVlepJEC_new/D");
    outTree_->Branch("yVlepJEC_new", &yVlepJEC_new, "yVlepJEC_new/D");
    outTree_->Branch("phiVlepJEC_new", &phiVlepJEC_new, "phiVlepJEC_new/D");
    outTree_->Branch("massVlepJEC_new", &massVlepJEC_new, "massVlepJEC_new/D");
    outTree_->Branch("mtVlepJECnew_new", &mtVlepJECnew_new, "mtVlepJECnew_new/D");
    outTree_->Branch("ptVlepJEC_JEC_up", &ptVlepJEC_JEC_up, "ptVlepJEC_JEC_up/D");
    outTree_->Branch("yVlepJEC_JEC_up", &yVlepJEC_JEC_up, "yVlepJEC_JEC_up/D");
    outTree_->Branch("phiVlepJEC_JEC_up", &phiVlepJEC_JEC_up, "phiVlepJEC_JEC_up/D");
    outTree_->Branch("massVlepJEC_JEC_up", &massVlepJEC_JEC_up, "massVlepJEC_JEC_up/D");
    outTree_->Branch("mtVlepJECnew_JEC_up", &mtVlepJECnew_JEC_up, "mtVlepJECnew_JEC_up/D");
    outTree_->Branch("ptVlepJEC_JEC_down", &ptVlepJEC_JEC_down, "ptVlepJEC_JEC_down/D");
    outTree_->Branch("yVlepJEC_JEC_down", &yVlepJEC_JEC_down, "yVlepJEC_JEC_down/D");
    outTree_->Branch("phiVlepJEC_JEC_down", &phiVlepJEC_JEC_down, "phiVlepJEC_JEC_down/D");
    outTree_->Branch("massVlepJEC_JEC_down", &massVlepJEC_JEC_down, "massVlepJEC_JEC_down/D");
    outTree_->Branch("mtVlepJECnew_JEC_down", &mtVlepJECnew_JEC_down, "mtVlepJECnew_JEC_down/D");
    outTree_->Branch("ptVlepJEC_JER_up", &ptVlepJEC_JER_up, "ptVlepJEC_JER_up/D");
    outTree_->Branch("yVlepJEC_JER_up", &yVlepJEC_JER_up, "yVlepJEC_JER_up/D");
    outTree_->Branch("phiVlepJEC_JER_up", &phiVlepJEC_JER_up, "phiVlepJEC_JER_up/D");
    outTree_->Branch("massVlepJEC_JER_up", &massVlepJEC_JER_up, "massVlepJEC_JER_up/D");
    outTree_->Branch("mtVlepJECnew_JER_up", &mtVlepJECnew_JER_up, "mtVlepJECnew_JER_up/D");
    outTree_->Branch("ptVlepJEC_JER_down", &ptVlepJEC_JER_down, "ptVlepJEC_JER_down/D");
    outTree_->Branch("yVlepJEC_JER_down", &yVlepJEC_JER_down, "yVlepJEC_JER_down/D");
    outTree_->Branch("phiVlepJEC_JER_down", &phiVlepJEC_JER_down, "phiVlepJEC_JER_down/D");
    outTree_->Branch("massVlepJEC_JER_down", &massVlepJEC_JER_down, "massVlepJEC_JER_down/D");
    outTree_->Branch("mtVlepJECnew_JER_down", &mtVlepJECnew_JER_down, "mtVlepJECnew_JER_down/D");
    outTree_->Branch("Mla", &Mla, "Mla/D");
    outTree_->Branch("Mla_f", &Mla_f, "Mla_f/D");
    outTree_->Branch("Mva", &Mva, "Mva/D");
    outTree_->Branch("Mva_f", &Mva_f, "Mva_f/D");
    outTree_->Branch("nlooseeles", &nlooseeles, "nlooseeles/I");
    outTree_->Branch("nloosemus", &nloosemus, "nloosemus/I");
    
    outTree_->Branch("genphoton_pt",genphoton_pt,"genphoton_pt[6]/D");
    outTree_->Branch("genphoton_eta",genphoton_eta,"genphoton_eta[6]/D");
    outTree_->Branch("genphoton_phi",genphoton_phi,"genphoton_phi[6]/D");
    outTree_->Branch("genmuon_pt",genmuon_pt,"genmuon_pt[6]/D");
    outTree_->Branch("genmuon_eta",genmuon_eta,"genmuon_eta[6]/D");
    outTree_->Branch("genmuon_phi", genmuon_phi,"genmuon_phi[6]/D");
    outTree_->Branch("genelectron_pt", genelectron_pt,"genelectron_pt[6]/D");
    outTree_->Branch("genelectron_eta", genelectron_eta,"genelectron_eta[6]/D");
    outTree_->Branch("genelectron_phi", genelectron_phi,"genelectron_phi[6]/D");

    /// Photon
    outTree_->Branch("photon_pt", photon_pt, "photon_pt[6]/D");
    outTree_->Branch("photon_eta", photon_eta, "photon_eta[6]/D");
    outTree_->Branch("photon_phi", photon_phi, "photon_phi[6]/D");
    outTree_->Branch("photon_e", photon_e, "photon_e[6]/D");
    outTree_->Branch("photonsc_eta", photonsc_eta, "photonsc_eta[6]/D");
    outTree_->Branch("photonsc_phi", photonsc_phi, "photonsc_phi[6]/D");
    outTree_->Branch("photon_pev", photon_pev, "photon_pev[6]/O");
    outTree_->Branch("photon_pevnew", photon_pevnew, "photon_pevnew[6]/O");
    outTree_->Branch("photon_ppsv", photon_ppsv, "photon_ppsv[6]/O");
    outTree_->Branch("photon_iseb", photon_iseb,"photon_iseb[6]/O");
    outTree_->Branch("photon_isee", photon_isee,"photon_isee[6]/O");
    outTree_->Branch("photon_hoe", photon_hoe, "photon_hoe[6]/D");
    outTree_->Branch("photon_sieie", photon_sieie, "photon_sieie[6]/D");
    outTree_->Branch("photon_sieie2", photon_sieie2, "photon_sieie2[6]/D");
    outTree_->Branch("photon_chiso", photon_chiso, "photon_chiso[6]/D");
    outTree_->Branch("photon_nhiso", photon_nhiso, "photon_nhiso[6]/D");
    outTree_->Branch("photon_phoiso", photon_phoiso, "photon_phoiso[6]/D");
    outTree_->Branch("photon_istrue", photon_istrue,"photon_istrue[6]/I");
    outTree_->Branch("photon_isprompt", photon_isprompt, "photon_isprompt[6]/I");
    outTree_->Branch("photon_drla", photon_drla, "photon_drla[6]/D");
    outTree_->Branch("photon_mla", photon_mla, "photon_mla[6]/D");
    outTree_->Branch("photon_mva", photon_mva,"photon_mva[6]/D");
    outTree_->Branch("passEleVeto", &passEleVeto,"passEleVeto/O");
    outTree_->Branch("passEleVetonew", &passEleVetonew,"passEleVetonew/O");
    outTree_->Branch("passPixelSeedVeto", &passPixelSeedVeto,"passPixelSeedVeto/O");
    outTree_->Branch("photonhaspixelseed", &photonhaspixelseed, "photonhaspixelseed/O");
    outTree_->Branch("photonhaspixelseed_f", &photonhaspixelseed_f, "photonhaspixelseed_f/O");
    outTree_->Branch("photonpasseleveto", &photonpasseleveto, "photonpasseleveto/O");
    outTree_->Branch("photonpasseleveto_f", &photonpasseleveto_f, "photonpasseleveto_F/O");
    outTree_->Branch("photonet", &photonet, "photonet/D");
    outTree_->Branch("photonet_f", &photonet_f, "photonet_f/D");
    outTree_->Branch("photoneta", &photoneta, "photoneta/D");
    outTree_->Branch("photoneta_f", &photoneta_f, "photoneta_f/D");
    outTree_->Branch("photonphi", &photonphi, "photonphi/D");
    outTree_->Branch("photonphi_f", &photonphi_f, "photonphi_f/D");
    outTree_->Branch("photone", &photone, "photone/D");
    outTree_->Branch("photone_f", &photone_f, "photone_f/D");
    outTree_->Branch("photonsceta", &photonsceta, "photonsceta/D");
    outTree_->Branch("photonsceta_f", &photonsceta_f, "photonsceta_f/D");
    outTree_->Branch("photonscphi", &photonscphi, "photonscphi/D");
    outTree_->Branch("photonscphi_f", &photonscphi_f, "photonscphi_f/D");
    
  	outTree_->Branch("photonsieie",&photonsieie,"photonsieie/D");
  	outTree_->Branch("photonsieie_f",&photonsieie_f,"photonsieie_f/D");
	outTree_->Branch("photonphoiso",&photonphoiso,"photonphoiso/D");
  	outTree_->Branch("photonphoiso_f",&photonphoiso_f,"photonphoiso_f/D");
  	outTree_->Branch("photonchiso",&photonchiso,"photonchiso/D");
  	outTree_->Branch("photonchiso_f",&photonchiso_f,"photonchiso_f/D");
  	outTree_->Branch("photonnhiso",&photonnhiso,"photonnhiso/D");
  	outTree_->Branch("photonnhiso_f",&photonnhiso_f,"photonnhiso_f/D");

    outTree_->Branch("iphoton", &iphoton, "iphoton/I");
    outTree_->Branch("iphoton_f", &iphoton_f, "iphoton_f/I");
    outTree_->Branch("drla", &drla, "drla/D");
    outTree_->Branch("drla_f", &drla_f, "drla_f/D");
    //photon gen match
    //    outTree_->Branch("dR"    , &dR_, "dR/D");
    //    outTree_->Branch("ISRPho"        , &ISRPho       ,"ISRPho/O"       );
    outTree_->Branch("isTrue", &isTrue_, "isTrue/I");
    outTree_->Branch("isprompt", &isprompt_, "isprompt/I");
    outTree_->Branch("ispromptLep", &ispromptLep_, "ispromptLep/I");
    //jets
    outTree_->Branch("ak4jet_hf", ak4jet_hf, "ak4jet_hf[6]/I");
    outTree_->Branch("ak4jet_pf", ak4jet_pf, "ak4jet_pf[6]/I");
    outTree_->Branch("jet1hf", &jet1hf, "jet1hf/I");
    outTree_->Branch("jet1pf", &jet1pf, "jet1pf/I");
    outTree_->Branch("jet2hf", &jet2hf, "jet2hf/I");
    outTree_->Branch("jet2pf", &jet2pf, "jet2pf/I");
    outTree_->Branch("jet1hf_f", &jet1hf_f, "jet1hf_f/I");
    outTree_->Branch("jet1pf_f", &jet1pf_f, "jet1pf_f/I");
    outTree_->Branch("jet2hf_f", &jet2hf_f, "jet2hf_f/I");
    outTree_->Branch("jet2pf_f", &jet2pf_f, "jet2pf_f/I");

    outTree_->Branch("ak4jet_pt", ak4jet_pt, "ak4jet_pt[6]/D");
    outTree_->Branch("ak4jet_eta", ak4jet_eta, "ak4jet_eta[6]/D");
    outTree_->Branch("ak4jet_phi", ak4jet_phi, "ak4jet_phi[6]/D");
    outTree_->Branch("ak4jet_e", ak4jet_e, "ak4jet_e[6]/D");
    outTree_->Branch("ak4jet_csv", ak4jet_csv, "ak4jet_csv[6]/D");
    outTree_->Branch("ak4jet_icsv", ak4jet_icsv, "ak4jet_icsv[6]/D");

    outTree_->Branch("ak4jet_pt_old", ak4jet_pt_old, "ak4jet_pt_old[6]/D");
    outTree_->Branch("ak4jet_pt_new", ak4jet_pt_new, "ak4jet_pt_new[6]/D");
    outTree_->Branch("ak4jet_pt_JEC_up", ak4jet_pt_JEC_up, "ak4jet_pt_JEC_up[6]/D");
    outTree_->Branch("ak4jet_pt_JEC_down", ak4jet_pt_JEC_down, "ak4jet_pt_JEC_down[6]/D");
    outTree_->Branch("ak4jet_pt_JER_up", ak4jet_pt_JER_up, "ak4jet_pt_JER_up[6]/D");
    outTree_->Branch("ak4jet_pt_JER_down", ak4jet_pt_JER_down, "ak4jet_pt_JER_down[6]/D");
    outTree_->Branch("ak4jet_eta", ak4jet_eta, "ak4jet_eta[6]/D");
    outTree_->Branch("ak4jet_phi", ak4jet_phi, "ak4jet_phi[6]/D");
    outTree_->Branch("ak4jet_e_old", ak4jet_e_old, "ak4jet_e_old[6]/D");
    outTree_->Branch("ak4jet_e_new", ak4jet_e_new, "ak4jet_e_new[6]/D");
    outTree_->Branch("ak4jet_e_JEC_up", ak4jet_e_JEC_up, "ak4jet_e_JEC_up[6]/D");
    outTree_->Branch("ak4jet_e_JEC_down", ak4jet_e_JEC_down, "ak4jet_e_JEC_down[6]/D");
    outTree_->Branch("ak4jet_e_JER_up", ak4jet_e_JER_up, "ak4jet_e_JER_up[6]/D");
    outTree_->Branch("ak4jet_e_JER_down", ak4jet_e_JER_down, "ak4jet_e_JER_down[6]/D");
    outTree_->Branch("ak4jet_csv", ak4jet_csv, "ak4jet_csv[6]/D");
    outTree_->Branch("ak4jet_icsv", ak4jet_icsv, "ak4jet_icsv[6]/D");
    outTree_->Branch("jet1pt", &jet1pt, "jet1pt/D");
    outTree_->Branch("jet1pt_new", &jet1pt_new, "jet1pt_new/D");
    outTree_->Branch("jet1pt_JEC_up", &jet1pt_JEC_up, "jet1pt_JEC_up/D");
    outTree_->Branch("jet1pt_JER_up", &jet1pt_JER_up, "jet1pt_JER_up/D");
    outTree_->Branch("jet1pt_JEC_down", &jet1pt_JEC_down, "jet1pt_JEC_down/D");
    outTree_->Branch("jet1pt_JER_down", &jet1pt_JER_down, "jet1pt_JER_down/D");
    outTree_->Branch("jet1pt_f", &jet1pt_f, "jet1pt_f/D");
    outTree_->Branch("jet1pt_new_f", &jet1pt_new_f, "jet1pt_new_f/D");
    //outTree_->Branch("jet1pt_JEC_up_f", &jet1pt_JEC_up_f, "jet1pt_JEC_up_f/D");
    //outTree_->Branch("jet1pt_JER_up_f", &jet1pt_JER_up_f, "jet1pt_JER_up_f/D");
    //outTree_->Branch("jet1pt_JEC_down_f", &jet1pt_JEC_down_f, "jet1pt_JEC_down_f/D");
    //outTree_->Branch("jet1pt_JER_down_f", &jet1pt_JER_down_f, "jet1pt_JER_down_f/D");
    outTree_->Branch("jet1eta", &jet1eta, "jet1eta/D");
    outTree_->Branch("jet1eta_new", &jet1eta_new, "jet1eta_new/D");
    outTree_->Branch("jet1eta_JEC_up", &jet1eta_JEC_up, "jet1eta_JEC_up/D");
    outTree_->Branch("jet1eta_JEC_down", &jet1eta_JEC_down, "jet1eta_JEC_down/D");
    outTree_->Branch("jet1eta_JER_up", &jet1eta_JER_up, "jet1eta_JER_up/D");
    outTree_->Branch("jet1eta_JER_down", &jet1eta_JER_down, "jet1eta_JER_down/D");
    outTree_->Branch("jet1eta_f", &jet1eta_f, "jet1eta_f/D");
    outTree_->Branch("jet1eta_new_f", &jet1eta_new_f, "jet1eta_new_f/D");
    //outTree_->Branch("jet1eta_JEC_up_f", &jet1eta_JEC_up_f, "jet1eta_JEC_up_f/D");
    //outTree_->Branch("jet1eta_JEC_down_f", &jet1eta_JEC_down_f, "jet1eta_JEC_down_f/D");
    //outTree_->Branch("jet1eta_JER_up_f", &jet1eta_JER_up_f, "jet1eta_JER_up_f/D");
    //outTree_->Branch("jet1eta_JER_down_f", &jet1eta_JER_down_f, "jet1eta_JER_down_f/D");
    outTree_->Branch("jet1phi", &jet1phi, "jet1phi/D");
    outTree_->Branch("jet1phi_new", &jet1phi_new, "jet1phi_new/D");
    outTree_->Branch("jet1phi_JEC_up", &jet1phi_JEC_up, "jet1phi_JEC_up/D");
    outTree_->Branch("jet1phi_JEC_down", &jet1phi_JEC_down, "jet1phi_JEC_down/D");
    outTree_->Branch("jet1phi_JER_up", &jet1phi_JER_up, "jet1phi_JER_up/D");
    outTree_->Branch("jet1phi_JER_down", &jet1phi_JER_down, "jet1phi_JER_down/D");
    outTree_->Branch("jet1phi_f", &jet1phi_f, "jet1phi_f/D");
    outTree_->Branch("jet1phi_new_f", &jet1phi_new_f, "jet1phi_new_f/D");
    //outTree_->Branch("jet1phi_JEC_up_f", &jet1phi_JEC_up_f, "jet1phi_JEC_up_f/D");
    //outTree_->Branch("jet1phi_JEC_down_f", &jet1phi_JEC_down_f, "jet1phi_JEC_down_f/D");
    //outTree_->Branch("jet1phi_JER_up_f", &jet1phi_JER_up_f, "jet1phi_JER_up_f/D");
    //outTree_->Branch("jet1phi_JER_down_f", &jet1phi_JER_down_f, "jet1phi_JER_down_f/D");
    outTree_->Branch("jet1e", &jet1e, "jet1e/D");
    outTree_->Branch("jet1e_new", &jet1e_new, "jet1e_new/D");
    outTree_->Branch("jet1e_JEC_up", &jet1e_JEC_up, "jet1e_JEC_up/D");
    outTree_->Branch("jet1e_JER_up", &jet1e_JER_up, "jet1e_JER_up/D");
    outTree_->Branch("jet1e_JEC_down", &jet1e_JEC_down, "jet1e_JEC_down/D");
    outTree_->Branch("jet1e_JER_down", &jet1e_JER_down, "jet1e_JER_down/D");
    outTree_->Branch("jet1e_f", &jet1e_f, "jet1e_f/D");
    outTree_->Branch("jet1e_new_f", &jet1e_new_f, "jet1e_new_f/D");
    //outTree_->Branch("jet1e_JEC_up_f", &jet1e_JEC_up_f, "jet1e_JEC_up_f/D");
    //outTree_->Branch("jet1e_JER_up_f", &jet1e_JER_up_f, "jet1e_JER_up_f/D");
    //outTree_->Branch("jet1e_JEC_down_f", &jet1e_JEC_down_f, "jet1e_JEC_down_f/D");
    //outTree_->Branch("jet1e_JER_down_f", &jet1e_JER_down_f, "jet1e_JER_down_f/D");
    outTree_->Branch("jet1csv", &jet1csv, "jet1csv/D");
    outTree_->Branch("jet1csv_new", &jet1csv_new, "jet1csv_new/D");
    outTree_->Branch("jet1csv_JEC_up", &jet1csv_JEC_up, "jet1csv_JEC_up/D");
    outTree_->Branch("jet1csv_JER_up", &jet1csv_JER_up, "jet1csv_JER_up/D");
    outTree_->Branch("jet1csv_JEC_down", &jet1csv_JEC_down, "jet1csv_JEC_down/D");
    outTree_->Branch("jet1csv_JER_down", &jet1csv_JER_down, "jet1csv_JER_down/D");
    outTree_->Branch("jet1csv_f", &jet1csv_f, "jet1csv_f/D");
    outTree_->Branch("jet1csv_new_f", &jet1csv_new_f, "jet1csv_new_f/D");
    //outTree_->Branch("jet1csv_JEC_up_f", &jet1csv_JEC_up_f, "jet1csv_JEC_up_f/D");
    //outTree_->Branch("jet1csv_JER_up_f", &jet1csv_JER_up_f, "jet1csv_JER_up_f/D");
    //outTree_->Branch("jet1csv_JEC_down_f", &jet1csv_JEC_down_f, "jet1csv_JEC_down_f/D");
    //outTree_->Branch("jet1csv_JER_down_f", &jet1csv_JER_down_f, "jet1csv_JER_down_f/D");
    outTree_->Branch("jet1icsv", &jet1icsv, "jet1icsv/D");
    outTree_->Branch("jet1icsv_new", &jet1icsv_new, "jet1icsv_new/D");
    outTree_->Branch("jet1icsv_JEC_up", &jet1icsv_JEC_up, "jet1icsv_JEC_up/D");
    outTree_->Branch("jet1icsv_JER_up", &jet1icsv_JER_up, "jet1icsv_JER_up/D");
    outTree_->Branch("jet1icsv_JEC_down", &jet1icsv_JEC_down, "jet1icsv_JEC_down/D");
    outTree_->Branch("jet1icsv_JER_down", &jet1icsv_JER_down, "jet1icsv_JER_down/D");
    outTree_->Branch("jet1icsv_f", &jet1icsv_f, "jet1icsv_f/D");
    outTree_->Branch("jet1icsv_new_f", &jet1icsv_new_f, "jet1icsv_new_f/D");
    //outTree_->Branch("jet1icsv_JEC_up_f", &jet1icsv_JEC_up_f, "jet1icsv_JEC_up_f/D");
    //outTree_->Branch("jet1icsv_JER_up_f", &jet1icsv_JER_up_f, "jet1icsv_JER_up_f/D");
    //outTree_->Branch("jet1icsv_JEC_down_f", &jet1icsv_JEC_down_f, "jet1icsv_JEC_down_f/D");
    //outTree_->Branch("jet1icsv_JER_down_f", &jet1icsv_JER_down_f, "jet1icsv_JER_down_f/D");
    outTree_->Branch("jet2pt", &jet2pt, "jet2pt/D");
    outTree_->Branch("jet2pt_new", &jet2pt_new, "jet2pt_new/D");
    outTree_->Branch("jet2pt_JEC_up", &jet2pt_JEC_up, "jet2pt_JEC_up/D");
    outTree_->Branch("jet2pt_JER_up", &jet2pt_JER_up, "jet2pt_JER_up/D");
    outTree_->Branch("jet2pt_JEC_down", &jet2pt_JEC_down, "jet2pt_JEC_down/D");
    outTree_->Branch("jet2pt_JER_down", &jet2pt_JER_down, "jet2pt_JER_down/D");
    outTree_->Branch("jet2pt_f", &jet2pt_f, "jet2pt_f/D");
    outTree_->Branch("jet2pt_new_f", &jet2pt_new_f, "jet2pt_new_f/D");
    //outTree_->Branch("jet2pt_JEC_up_f", &jet2pt_JEC_up_f, "jet2pt_JEC_up_f/D");
    //outTree_->Branch("jet2pt_JER_up_f", &jet2pt_JER_up_f, "jet2pt_JER_up_f/D");
    //outTree_->Branch("jet2pt_JEC_down_f", &jet2pt_JEC_down_f, "jet2pt_JEC_down_f/D");
    //outTree_->Branch("jet2pt_JER_down_f", &jet2pt_JER_down_f, "jet2pt_JER_down_f/D");
    outTree_->Branch("jet2eta", &jet2eta, "jet2eta/D");
    outTree_->Branch("jet2eta_new", &jet2eta_new, "jet2eta_new/D");
    outTree_->Branch("jet2eta_JEC_up", &jet2eta_JEC_up, "jet2eta_JEC_up/D");
    outTree_->Branch("jet2eta_JEC_down", &jet2eta_JEC_down, "jet2eta_JEC_down/D");
    outTree_->Branch("jet2eta_JER_up", &jet2eta_JER_up, "jet2eta_JER_up/D");
    outTree_->Branch("jet2eta_JER_down", &jet2eta_JER_down, "jet2eta_JER_down/D");
    outTree_->Branch("jet2phi", &jet2phi, "jet2phi/D");
    outTree_->Branch("jet2phi_new", &jet2phi_new, "jet2phi_new/D");
    outTree_->Branch("jet2phi_JEC_up", &jet2phi_JEC_up, "jet2phi_JEC_up/D");
    outTree_->Branch("jet2phi_JEC_down", &jet2phi_JEC_down, "jet2phi_JEC_down/D");
    outTree_->Branch("jet2phi_JER_up", &jet2phi_JER_up, "jet2phi_JER_up/D");
    outTree_->Branch("jet2phi_JER_down", &jet2phi_JER_down, "jet2phi_JER_down/D");
    outTree_->Branch("jet2phi_f", &jet2phi_f, "jet2phi_f/D");
    outTree_->Branch("jet2phi_new_f", &jet2phi_new_f, "jet2phi_new_f/D");
    //outTree_->Branch("jet2phi_JEC_up_f", &jet2phi_JEC_up_f, "jet2phi_JEC_up_f/D");
    //outTree_->Branch("jet2phi_JEC_down_f", &jet2phi_JEC_down_f, "jet2phi_JEC_down_f/D");
    //outTree_->Branch("jet2phi_JER_up_f", &jet2phi_JER_up_f, "jet2phi_JER_up_f/D");
    //outTree_->Branch("jet2phi_JER_down_f", &jet2phi_JER_down_f, "jet2phi_JER_down_f/D");
	outTree_->Branch("jet2eta_f", &jet2eta_f, "jet2eta_f/D");
    outTree_->Branch("jet2eta_new_f", &jet2eta_new_f, "jet2eta_new_f/D");
    //outTree_->Branch("jet2eta_JEC_up_f", &jet2eta_JEC_up_f, "jet2eta_JEC_up_f/D");
    //outTree_->Branch("jet2eta_JEC_down_f", &jet2eta_JEC_down_f, "jet2eta_JEC_down_f/D");
    //outTree_->Branch("jet2eta_JER_up_f", &jet2eta_JER_up_f, "jet2eta_JER_up_f/D");
    //outTree_->Branch("jet2eta_JER_down_f", &jet2eta_JER_down_f, "jet2eta_JER_down_f/D");
    outTree_->Branch("jet2e", &jet2e, "jet2e/D");
    outTree_->Branch("jet2e_new", &jet2e_new, "jet2e_new/D");
    outTree_->Branch("jet2e_JEC_up", &jet2e_JEC_up, "jet2e_JEC_up/D");
    outTree_->Branch("jet2e_JER_up", &jet2e_JER_up, "jet2e_JER_up/D");
    outTree_->Branch("jet2e_JEC_down", &jet2e_JEC_down, "jet2e_JEC_down/D");
    outTree_->Branch("jet2e_JER_down", &jet2e_JER_down, "jet2e_JER_down/D");
    outTree_->Branch("jet2e_f", &jet2e_f, "jet2e_f/D");
    outTree_->Branch("jet2e_new_f", &jet2e_new_f, "jet2e_new_f/D");
    //outTree_->Branch("jet2e_JEC_up_f", &jet2e_JEC_up_f, "jet2e_JEC_up_f/D");
    //outTree_->Branch("jet2e_JER_up_f", &jet2e_JER_up_f, "jet2e_JER_up_f/D");
    //outTree_->Branch("jet2e_JEC_down_f", &jet2e_JEC_down_f, "jet2e_JEC_down_f/D");
    //outTree_->Branch("jet2e_JER_down_f", &jet2e_JER_down_f, "jet2e_JER_down_f/D");
    outTree_->Branch("jet2csv", &jet2csv, "jet2csv/D");
    outTree_->Branch("jet2csv_new", &jet2csv_new, "jet2csv_new/D");
    outTree_->Branch("jet2csv_JEC_up", &jet2csv_JEC_up, "jet2csv_JEC_up/D");
    outTree_->Branch("jet2csv_JER_up", &jet2csv_JER_up, "jet2csv_JER_up/D");
    outTree_->Branch("jet2csv_JEC_down", &jet2csv_JEC_down, "jet2csv_JEC_down/D");
    outTree_->Branch("jet2csv_JER_down", &jet2csv_JER_down, "jet2csv_JER_down/D");
    outTree_->Branch("jet2csv_f", &jet2csv_f, "jet2csv_f/D");
    //outTree_->Branch("jet2csv_new_f", &jet2csv_new_f, "jet2csv_new_f/D");
    //outTree_->Branch("jet2csv_JEC_up_f", &jet2csv_JEC_up_f, "jet2csv_JEC_up_f/D");
    //outTree_->Branch("jet2csv_JER_up_f", &jet2csv_JER_up_f, "jet2csv_JER_up_f/D");
    //outTree_->Branch("jet2csv_JEC_down_f", &jet2csv_JEC_down_f, "jet2csv_JEC_down_f/D");
    //outTree_->Branch("jet2csv_JER_down_f", &jet2csv_JER_down_f, "jet2csv_JER_down_f/D");
    outTree_->Branch("jet2icsv", &jet2icsv, "jet2icsv/D");
    outTree_->Branch("jet2icsv_new", &jet2icsv_new, "jet2icsv_new/D");
    outTree_->Branch("jet2icsv_JEC_up", &jet2icsv_JEC_up, "jet2icsv_JEC_up/D");
    outTree_->Branch("jet2icsv_JER_up", &jet2icsv_JER_up, "jet2icsv_JER_up/D");
    outTree_->Branch("jet2icsv_JEC_down", &jet2icsv_JEC_down, "jet2icsv_JEC_down/D");
    outTree_->Branch("jet2icsv_JER_down", &jet2icsv_JER_down, "jet2icsv_JER_down/D");
    outTree_->Branch("jet2icsv_f", &jet2icsv_f, "jet2icsv_f/D");
    outTree_->Branch("jet2icsv_new_f", &jet2icsv_new_f, "jet2icsv_new_f/D");
    //outTree_->Branch("jet2icsv_JEC_up_f", &jet2icsv_JEC_up_f, "jet2icsv_JEC_up_f/D");
    //outTree_->Branch("jet2icsv_JER_up_f", &jet2icsv_JER_up_f, "jet2icsv_JER_up_f/D");
    //outTree_->Branch("jet2icsv_JEC_down_f", &jet2icsv_JEC_down_f, "jet2icsv_JEC_down_f/D");
    //outTree_->Branch("jet2icsv_JER_down_f", &jet2icsv_JER_down_f, "jet2icsv_JER_down_f/D");
    outTree_->Branch("drj1a", &drj1a, "drj1a/D");
    outTree_->Branch("drj1a_new", &drj1a_new, "drj1a_new/D");
    outTree_->Branch("drj1a_JEC_up", &drj1a_JEC_up, "drj1a_JEC_up/D");
    outTree_->Branch("drj1a_JEC_down", &drj1a_JEC_down, "drj1a_JEC_down/D");
    outTree_->Branch("drj1a_JER_up", &drj1a_JER_up, "drj1a_JER_up/D");
    outTree_->Branch("drj1a_JER_down", &drj1a_JER_down, "drj1a_JER_down/D");
    outTree_->Branch("drj1a_f", &drj1a_f, "drj1a_f/D");
    outTree_->Branch("drj1a_new_f", &drj1a_new_f, "drj1a_new_f/D");
    //outTree_->Branch("drj1a_JEC_up_f", &drj1a_JEC_up_f, "drj1a_JEC_up_f/D");
    //outTree_->Branch("drj1a_JEC_down_f", &drj1a_JEC_down_f, "drj1a_JEC_down_f/D");
    //outTree_->Branch("drj1a_JER_up_f", &drj1a_JER_up_f, "drj1a_JER_up_f/D");
    //outTree_->Branch("drj1a_JER_down_f", &drj1a_JER_down_f, "drj1a_JER_down_f/D");
    outTree_->Branch("drj2a", &drj2a, "drj2a/D");
    outTree_->Branch("drj2a_new", &drj2a_new, "drj2a_new/D");
    outTree_->Branch("drj2a_JEC_up", &drj2a_JEC_up, "drj2a_JEC_up/D");
    outTree_->Branch("drj2a_JEC_down", &drj2a_JEC_down, "drj2a_JEC_down/D");
    outTree_->Branch("drj2a_JER_up", &drj2a_JER_up, "drj2a_JER_up/D");
    outTree_->Branch("drj2a_JER_down", &drj2a_JER_down, "drj2a_JER_down/D");
    outTree_->Branch("drj2a_f", &drj2a_f, "drj2a_f/D");
    outTree_->Branch("drj2a_new_f", &drj2a_new_f, "drj2a_new_f/D");
    //outTree_->Branch("drj2a_JEC_up_f", &drj2a_JEC_up_f, "drj2a_JEC_up_f/D");
    //outTree_->Branch("drj2a_JEC_down_f", &drj2a_JEC_down_f, "drj2a_JEC_down_f/D");
    //outTree_->Branch("drj2a_JER_up_f", &drj2a_JER_up_f, "drj2a_JER_up_f/D");
    //outTree_->Branch("drj2a_JER_down_f", &drj2a_JER_down_f, "drj2a_JER_down_f/D");
    outTree_->Branch("drj1l", &drj1l, "drj1l/D");
    outTree_->Branch("drj1l_new", &drj1l_new, "drj1l_new/D");
    outTree_->Branch("drj1l_JEC_up", &drj1l_JEC_up, "drj1l_JEC_up/D");
    outTree_->Branch("drj1l_JEC_down", &drj1l_JEC_down, "drj1l_JEC_down/D");
    outTree_->Branch("drj1l_JER_up", &drj1l_JER_up, "drj1l_JER_up/D");
    outTree_->Branch("drj1l_JER_down", &drj1l_JER_down, "drj1l_JER_down/D");
    outTree_->Branch("drj1l_f", &drj1l_f, "drj1l_f/D");
    outTree_->Branch("drj1l_new_f", &drj1l_new_f, "drj1l_new_f/D");
    //outTree_->Branch("drj1l_JEC_up_f", &drj1l_JEC_up_f, "drj1l_JEC_up_f/D");
    //outTree_->Branch("drj1l_JEC_down_f", &drj1l_JEC_down_f, "drj1l_JEC_down_f/D");
    //outTree_->Branch("drj1l_JER_up_f", &drj1l_JER_up_f, "drj1l_JER_up_f/D");
    //outTree_->Branch("drj1l_JER_down_f", &drj1l_JER_down_f, "drj1l_JER_down_f/D");
    outTree_->Branch("drj2l", &drj2l, "drj2l/D");
    outTree_->Branch("drj2l_new", &drj2l_new, "drj2l_new/D");
    outTree_->Branch("drj2l_JEC_up", &drj2l_JEC_up, "drj2l_JEC_up/D");
    outTree_->Branch("drj2l_JEC_down", &drj2l_JEC_down, "drj2l_JEC_down/D");
    outTree_->Branch("drj2l_JER_up", &drj2l_JER_up, "drj2l_JER_up/D");
    outTree_->Branch("drj2l_JER_down", &drj2l_JER_down, "drj2l_JER_down/D");
    outTree_->Branch("drj2l_f", &drj2l_f, "drj2l_f/D");
    outTree_->Branch("drj2l_new_f", &drj2l_new_f, "drj2l_new_f/D");
    //outTree_->Branch("drj2l_JEC_up_f", &drj2l_JEC_up_f, "drj2l_JEC_up_f/D");
    //outTree_->Branch("drj2l_JEC_down_f", &drj2l_JEC_down_f, "drj2l_JEC_down_f/D");
    //outTree_->Branch("drj2l_JER_up_f", &drj2l_JER_up_f, "drj2l_JER_up_f/D");
    //outTree_->Branch("drj2l_JER_down_f", &drj2l_JER_down_f, "drj2l_JER_down_f/D");

    outTree_->Branch("Mjj", &Mjj, "Mjj/D");
    outTree_->Branch("Mjj_new", &Mjj_new, "Mjj_new/D");
    outTree_->Branch("Mjj_JEC_up", &Mjj_JEC_up, "Mjj_JEC_up/D");
    outTree_->Branch("Mjj_JEC_down", &Mjj_JEC_down, "Mjj_JEC_down/D");
    outTree_->Branch("Mjj_JER_up", &Mjj_JER_up, "Mjj_JER_up/D");
    outTree_->Branch("Mjj_JER_down", &Mjj_JER_down, "Mjj_JER_down/D");
    outTree_->Branch("Mjj_f", &Mjj_f, "Mjj_f/D");
    outTree_->Branch("Mjj_new_f", &Mjj_new_f, "Mjj_new_f/D");
    //outTree_->Branch("Mjj_JEC_up_f", &Mjj_JEC_up_f, "Mjj_JEC_up_f/D");
    //outTree_->Branch("Mjj_JEC_down_f", &Mjj_JEC_down_f, "Mjj_JEC_down_f/D");
    //outTree_->Branch("Mjj_JER_up_f", &Mjj_JER_up_f, "Mjj_JER_up_f/D");
    //outTree_->Branch("Mjj_JER_down_f", &Mjj_JER_down_f, "Mjj_JER_down_f/D");
    outTree_->Branch("deltaeta", &deltaeta, "deltaeta/D");
    outTree_->Branch("deltaeta_new", &deltaeta_new, "deltaeta_new/D");
    outTree_->Branch("deltaeta_JEC_up", &deltaeta_JEC_up, "deltaeta_JEC_up/D");
    outTree_->Branch("deltaeta_JEC_down", &deltaeta_JEC_down, "deltaeta_JEC_down/D");
    outTree_->Branch("deltaeta_JER_up", &deltaeta_JER_up, "deltaeta_JER_up/D");
    outTree_->Branch("deltaeta_JER_down", &deltaeta_JER_down, "deltaeta_JER_down/D");
    outTree_->Branch("deltaeta_f", &deltaeta_f, "deltaeta_f/D");
    outTree_->Branch("deltaeta_new_f", &deltaeta_new_f, "deltaeta_new_f/D");
    //outTree_->Branch("deltaeta_JEC_up_f", &deltaeta_JEC_up_f, "deltaeta_JEC_up_f/D");
    //outTree_->Branch("deltaeta_JEC_down_f", &deltaeta_JEC_down_f, "deltaeta_JEC_down_f/D");
    //outTree_->Branch("deltaeta_JER_up_f", &deltaeta_JER_up_f, "deltaeta_JER_up_f/D");
    //outTree_->Branch("deltaeta_JER_down_f", &deltaeta_JER_down_f, "deltaeta_JER_down_f/D");
    outTree_->Branch("zepp", &zepp, "zepp/D");
    outTree_->Branch("zepp_new", &zepp_new, "zepp_new/D");
    outTree_->Branch("zepp_JEC_up", &zepp_JEC_up, "zepp_JEC_up/D");
    outTree_->Branch("zepp_JEC_down", &zepp_JEC_down, "zepp_JEC_down/D");
    outTree_->Branch("zepp_JER_up", &zepp_JER_up, "zepp_JER_up/D");
    outTree_->Branch("zepp_JER_down", &zepp_JER_down, "zepp_JER_down/D");
    outTree_->Branch("zepp_f", &zepp_f, "zepp_f/D");
    outTree_->Branch("zepp_new_f", &zepp_new_f, "zepp_new_f/D");
    //outTree_->Branch("zepp_JEC_up_f", &zepp_JEC_up_f, "zepp_JEC_up_f/D");
    //outTree_->Branch("zepp_JEC_down_f", &zepp_JEC_down_f, "zepp_JEC_down_f/D");
    //outTree_->Branch("zepp_JER_up_f", &zepp_JER_up_f, "zepp_JER_up_f/D");
    //outTree_->Branch("zepp_JER_down_f", &zepp_JER_down_f, "zepp_JER_down_f/D");
    // Generic kinematic quantities
    outTree_->Branch("ptlep1", &ptlep1, "ptlep1/D");
    outTree_->Branch("etalep1", &etalep1, "etalep1/D");
    outTree_->Branch("philep1", &philep1, "philep1/D");
    outTree_->Branch("energylep1", &energylep1, "energylep1/D");
    outTree_->Branch("j1metPhi", &j1metPhi, "j1metPhi/D");
    outTree_->Branch("j1metPhi_new", &j1metPhi_new, "j1metPhi_new/D");
    outTree_->Branch("j1metPhi_JEC_up", &j1metPhi_JEC_up, "j1metPhi_JEC_up/D");
    outTree_->Branch("j1metPhi_JEC_down", &j1metPhi_JEC_down, "j1metPhi_JEC_down/D");
    outTree_->Branch("j1metPhi_JER_up", &j1metPhi_JER_up, "j1metPhi_JER_up/D");
    outTree_->Branch("j1metPhi_JER_down", &j1metPhi_JER_down, "j1metPhi_JER_down/D");
    outTree_->Branch("j1metPhi_f", &j1metPhi_f, "j1metPhi_f/D");
    outTree_->Branch("j1metPhi_new_f", &j1metPhi_new_f, "j1metPhi_new_f/D");
    //outTree_->Branch("j1metPhi_JEC_up_f", &j1metPhi_JEC_up_f, "j1metPhi_JEC_up_f/D");
    //outTree_->Branch("j1metPhi_JEC_down_f", &j1metPhi_JEC_down_f, "j1metPhi_JEC_down_f/D");
    ///outTree_->Branch("j1metPhi_JER_up_f", &j1metPhi_JER_up_f, "j1metPhi_JER_up_f/D");
    //outTree_->Branch("j1metPhi_JER_down_f", &j1metPhi_JER_down_f, "j1metPhi_JER_down_f/D");
    outTree_->Branch("j2metPhi", &j2metPhi, "j2metPhi/D");
    outTree_->Branch("j2metPhi_new", &j2metPhi_new, "j2metPhi_new/D");
    outTree_->Branch("j2metPhi_JEC_up", &j2metPhi_JEC_up, "j2metPhi_JEC_up/D");
    outTree_->Branch("j2metPhi_JEC_down", &j2metPhi_JEC_down, "j2metPhi_JEC_down/D");
    outTree_->Branch("j2metPhi_JER_up", &j2metPhi_JER_up, "j2metPhi_JER_up/D");
    outTree_->Branch("j2metPhi_JER_down", &j2metPhi_JER_down, "j2metPhi_JER_down/D");
    outTree_->Branch("j2metPhi_f", &j2metPhi_f, "j2metPhi_f/D");
    outTree_->Branch("j2metPhi_new_f", &j2metPhi_new_f, "j2metPhi_new_f/D");
    //outTree_->Branch("j2metPhi_JEC_up_f", &j2metPhi_JEC_up_f, "j2metPhi_JEC_up_f/D");
    //outTree_->Branch("j2metPhi_JEC_down_f", &j2metPhi_JEC_down_f, "j2metPhi_JEC_down_f/D");
    //outTree_->Branch("j2metPhi_JER_up_f", &j2metPhi_JER_up_f, "j2metPhi_JER_up_f/D");
    //outTree_->Branch("j2metPhi_JER_down_f", &j2metPhi_JER_down_f, "j2metPhi_JER_down_f/D");
    outTree_->Branch("Dphiwajj", &Dphiwajj, "Dphiwajj/D");
    outTree_->Branch("Dphiwajj_f", &Dphiwajj_f, "Dphiwajj_f/D");
    outTree_->Branch("Dphiwajj_new", &Dphiwajj_new, "Dphiwajj_new/D");
    outTree_->Branch("Dphiwajj_JEC_up", &Dphiwajj_JEC_up, "Dphiwajj_JEC_up/D");
    outTree_->Branch("Dphiwajj_JEC_down", &Dphiwajj_JEC_down, "Dphiwajj_JEC_down/D");
    outTree_->Branch("Dphiwajj_JER_up", &Dphiwajj_JER_up, "Dphiwajj_JER_up/D");
    outTree_->Branch("Dphiwajj_JER_down", &Dphiwajj_JER_down, "Dphiwajj_JER_down/D");
    // MET
    outTree_->Branch("METraw_et", &METraw_et, "METraw_et/D");
    outTree_->Branch("METraw_phi", &METraw_phi, "METraw_phi/D");
    outTree_->Branch("METraw_sumEt", &METraw_sumEt, "METraw_sumEt/D");
    outTree_->Branch("genMET", &genMET, "genMET/D");
    outTree_->Branch("MET_et", &MET_et, "MET_et/D");
    outTree_->Branch("MET_et_new", &MET_et_new, "MET_et_new/D");
    // Marked for debug
    outTree_->Branch("MET_et_JEC_up", &MET_et_JEC_up, "MET_et_JEC_up/D");
    outTree_->Branch("MET_et_JEC_down", &MET_et_JEC_down, "MET_et_JEC_down/D");
    outTree_->Branch("MET_et_JER_up", &MET_et_JER_up, "MET_et_JER_up/D");
    outTree_->Branch("MET_et_JER_down", &MET_et_JER_down, "MET_et_JER_down/D");
    // Marked for debug
    outTree_->Branch("MET_phi", &MET_phi, "MET_phi/D");
    // Marked for debug
    outTree_->Branch("MET_phi_new", &MET_phi_new, "MET_phi_new/D");
    outTree_->Branch("MET_phi_JEC_up", &MET_phi_JEC_up, "MET_phi_JEC_up/D");
    outTree_->Branch("MET_phi_JEC_down", &MET_phi_JEC_down, "MET_phi_JEC_down/D");
    outTree_->Branch("MET_phi_JER_up", &MET_phi_JER_up, "MET_phi_JER_up/D");
    outTree_->Branch("MET_phi_JER_down", &MET_phi_JER_down, "MET_phi_JER_down/D");
    // Marked for debug
    //HLT bits
    outTree_->Branch("HLT_Ele1", &HLT_Ele1, "HLT_Ele1/I");
    outTree_->Branch("HLT_Ele2", &HLT_Ele2, "HLT_Ele2/I");
    outTree_->Branch("HLT_Mu1", &HLT_Mu1, "HLT_Mu1/I");
    outTree_->Branch("HLT_Mu2", &HLT_Mu2, "HLT_Mu2/I");
    outTree_->Branch("HLT_Mu3", &HLT_Mu3, "HLT_Mu3/I");
    // filter
    //  outTree_->Branch("triggerWeight"   ,&triggerWeight  ,"triggerWeight/D"  );
    outTree_->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    outTree_->Branch("pileupWeight", &pileupWeight, "pileupWeight/D");

    //L1 prefiring
  	outTree_->Branch("prefWeight"   ,&_prefiringweight,"prefWeight/D"  );
   	outTree_->Branch("prefWeightUp" ,&_prefiringweightup,"prefWeightUp/D"  );
 	outTree_->Branch("prefWeightDown",&_prefiringweightdown,"prefWeightDown/D"  );

}

PKUTreeMaker::~PKUTreeMaker()
{
}
#endif
