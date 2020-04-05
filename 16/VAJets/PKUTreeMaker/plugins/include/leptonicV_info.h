#ifndef _leptonicV_
#define _leptonicV_

void PKUTreeMaker::leptonicV_info(edm::Event const & iEvent) {

    edm::Handle<edm::View<reco::Candidate> > leptonicVs;
    iEvent.getByToken(leptonicVSrc_, leptonicVs);
	is_leptonicVs_Empty = leptonicVs->empty();	

    const reco::Candidate& leptonicV = leptonicVs->at(0);
    const reco::Candidate& metCand   = metHandle->at(0);
    const reco::Candidate& lepton    = (*leptonicV.daughter(0));
   
	edm::Handle<edm::View<pat::Muon>> goodmus;
	iEvent.getByToken(goodmuonToken_, goodmus);	
	
	edm::Handle<edm::View<pat::Electron>> goodeles;
	iEvent.getByToken(goodeleToken_, goodeles);

	edm::Handle<edm::View<pat::Muon>> loosemus;
	iEvent.getByToken(loosemuonToken_, loosemus);
	
	edm::Handle<edm::View<pat::Electron>> looseeles;
	iEvent.getByToken(looseelectronToken_, looseeles);
	
	edm::Handle<edm::View<reco::Candidate>> metHandle;
	iEvent.getByToken(metSrc_, metHandle);	

    triggerWeight       = 1.0;
    pileupWeight        = 1.0;
    double targetEvents = targetLumiInvPb_ * crossSectionPb_;
    lumiWeight          = targetEvents / originalNEvents_;
    lep                 = abs(leptonicV.daughter(0)->pdgId());

    ptVlep              = leptonicV.pt();
    yVlep               = leptonicV.eta();
    phiVlep             = leptonicV.phi();
    massVlep            = leptonicV.mass();
    mtVlep              = leptonicV.mt();
    ptlep1              = leptonicV.daughter(1)->pt();
    etalep1             = leptonicV.daughter(1)->eta();
    philep1             = leptonicV.daughter(1)->phi();
    energylep1          = leptonicV.daughter(1)->energy();
    if (leptonicV.daughter(0)->isElectron() || leptonicV.daughter(0)->isMuon()) {
        ptlep1     = leptonicV.daughter(0)->pt();
        etalep1    = leptonicV.daughter(0)->eta();
        philep1    = leptonicV.daughter(0)->phi();
        energylep1 = leptonicV.daughter(0)->energy();
    }

    met        = metCand.pt();
    metPhi     = metCand.phi();
    mtVlepnew  = sqrt(2 * ptlep1 * met * (1.0 - cos(philep1 - metPhi)));
    nlooseeles = looseeles->size();
    nloosemus  = loosemus->size();
	ngoodmus   = goodmus->size();
    ngoodeles   = goodeles->size();

    TLorentzVector glepton;
    glepton.SetPtEtaPhiE(ptlep1, etalep1, philep1, energylep1);
    math::XYZTLorentzVector     neutrinoP4 = getNeutrinoP4(MET_et, MET_phi, glepton, 1);
    math::XYZTLorentzVector     neutrinoP4_new = getNeutrinoP4(MET_et_new, MET_phi_new, glepton, 1);
    math::XYZTLorentzVector     neutrinoP4_JEC_up = getNeutrinoP4(MET_et_JEC_up, MET_phi_JEC_up, glepton, 1);
    math::XYZTLorentzVector     neutrinoP4_JEC_down = getNeutrinoP4(MET_et_JEC_down, MET_phi_JEC_down, glepton, 1);
    math::XYZTLorentzVector     neutrinoP4_JER_up = getNeutrinoP4(MET_et_JER_up, MET_phi_JER_up, glepton, 1);
    math::XYZTLorentzVector     neutrinoP4_JER_down = getNeutrinoP4(MET_et_JER_down, MET_phi_JER_down, glepton, 1);
    reco::CandidateBaseRef      METBaseRef = metHandle->refAt(0);  
    reco::ShallowCloneCandidate neutrino(METBaseRef, 0, neutrinoP4);
    reco::ShallowCloneCandidate neutrino_new(METBaseRef, 0, neutrinoP4_new);
    reco::ShallowCloneCandidate neutrino_JEC_up(METBaseRef, 0, neutrinoP4_JEC_up);
    reco::ShallowCloneCandidate neutrino_JEC_down(METBaseRef, 0, neutrinoP4_JEC_down);
    reco::ShallowCloneCandidate neutrino_JER_up(METBaseRef, 0, neutrinoP4_JER_up);
    reco::ShallowCloneCandidate neutrino_JER_down(METBaseRef, 0, neutrinoP4_JER_down);
    reco::CompositeCandidate    WLeptonic;
    reco::CompositeCandidate    WLeptonic_new;
    reco::CompositeCandidate    WLeptonic_JEC_up;
    reco::CompositeCandidate    WLeptonic_JEC_down;
    reco::CompositeCandidate    WLeptonic_JER_up;
    reco::CompositeCandidate    WLeptonic_JER_down;
    WLeptonic.addDaughter(lepton);
    WLeptonic.addDaughter(neutrino);
    WLeptonic_new.addDaughter(lepton);
    WLeptonic_new.addDaughter(neutrino_new);
    WLeptonic_JEC_up.addDaughter(lepton);
    WLeptonic_JEC_up.addDaughter(neutrino_JEC_up);
    WLeptonic_JEC_down.addDaughter(lepton);
    WLeptonic_JEC_down.addDaughter(neutrino_JEC_down);
    WLeptonic_JER_up.addDaughter(lepton);
    WLeptonic_JER_up.addDaughter(neutrino_JER_up);
    WLeptonic_JER_down.addDaughter(lepton);
    WLeptonic_JER_down.addDaughter(neutrino_JER_down);

    AddFourMomenta addP4;
    addP4.set(WLeptonic);
    AddFourMomenta addP4_new;
    addP4_new.set(WLeptonic_new);
    AddFourMomenta addP4_JEC_up;
    addP4_JEC_up.set(WLeptonic_JEC_up);
    AddFourMomenta addP4_JEC_down;
    addP4_JEC_down.set(WLeptonic_JEC_down);
    AddFourMomenta addP4_JER_up;
    addP4_JER_up.set(WLeptonic_JER_up);
    AddFourMomenta addP4_JER_down;
    addP4_JER_down.set(WLeptonic_JER_down);

    ptVlepJEC    = WLeptonic.pt();
    yVlepJEC     = WLeptonic.eta();
    phiVlepJEC   = WLeptonic.phi();
	energyVlepJEC= WLeptonic.energy();
    massVlepJEC  = WLeptonic.mass();
    mtVlepJEC    = WLeptonic.mt();
    mtVlepJECnew = sqrt(2 * ptlep1 * MET_et * (1.0 - cos(philep1 - MET_phi)));

    ptVlepJEC_new    = WLeptonic_new.pt();
    yVlepJEC_new     = WLeptonic_new.eta();
    phiVlepJEC_new   = WLeptonic_new.phi();
    massVlepJEC_new  = WLeptonic_new.mass();
    mtVlepJEC_new    = WLeptonic_new.mt();
    mtVlepJECnew_new = sqrt(2 * ptlep1 * MET_et_new * (1.0 - cos(philep1 - MET_phi_new)));

    ptVlepJEC_JEC_up    = WLeptonic_JEC_up.pt();
    yVlepJEC_JEC_up     = WLeptonic_JEC_up.eta();
    phiVlepJEC_JEC_up   = WLeptonic_JEC_up.phi();
    massVlepJEC_JEC_up  = WLeptonic_JEC_up.mass();
    mtVlepJEC_JEC_up    = WLeptonic_JEC_up.mt();
    mtVlepJECnew_JEC_up = sqrt(2 * ptlep1 * MET_et_JEC_up * (1.0 - cos(philep1 - MET_phi_JEC_up)));

    ptVlepJEC_JEC_down    = WLeptonic_JEC_down.pt();
    yVlepJEC_JEC_down     = WLeptonic_JEC_down.eta();
    phiVlepJEC_JEC_down   = WLeptonic_JEC_down.phi();
    massVlepJEC_JEC_down  = WLeptonic_JEC_down.mass();
    mtVlepJEC_JEC_down    = WLeptonic_JEC_down.mt();
    mtVlepJECnew_JEC_down = sqrt(2 * ptlep1 * MET_et_JEC_down * (1.0 - cos(philep1 - MET_phi_JEC_down)));

    ptVlepJEC_JER_up    = WLeptonic_JER_up.pt();
    yVlepJEC_JER_up     = WLeptonic_JER_up.eta();
    phiVlepJEC_JER_up   = WLeptonic_JER_up.phi();
    massVlepJEC_JER_up  = WLeptonic_JER_up.mass();
    mtVlepJEC_JER_up    = WLeptonic_JER_up.mt();
    mtVlepJECnew_JER_up = sqrt(2 * ptlep1 * MET_et_JER_up * (1.0 - cos(philep1 - MET_phi_JER_up)));

    ptVlepJEC_JER_down    = WLeptonic_JER_down.pt();
    yVlepJEC_JER_down     = WLeptonic_JER_down.eta();
    phiVlepJEC_JER_down   = WLeptonic_JER_down.phi();
    massVlepJEC_JER_down  = WLeptonic_JER_down.mass();
    mtVlepJEC_JER_down    = WLeptonic_JER_down.mt();
    mtVlepJECnew_JER_down = sqrt(2 * ptlep1 * MET_et_JER_down * (1.0 - cos(philep1 - MET_phi_JER_down)));



}

#endif
