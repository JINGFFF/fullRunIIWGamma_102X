#ifndef _Photon_
#define _Photon_

void PKUTreeMaker::photon_info(edm::Event const & iEvent) {
    edm::Handle<edm::View<reco::GenParticle>> genParticles;
    iEvent.getByToken(genSrc_, genParticles);

    iEvent.getByToken(rhoToken_, rho_);
    double fastJetRho = *(rho_.product());
    useless           = fastJetRho;

    edm::Handle<edm::View<pat::Photon>> photons;
    iEvent.getByToken(photonSrc_, photons);
	
    if (photons->empty()) {
        hasphoton = 0.;
    }
    else {
        hasphoton =1.;
    }

    TLorentzVector glepton;
    glepton.SetPtEtaPhiE(ptlep1, etalep1, philep1, energylep1);

    rhoVal_ = -99.;
    rhoVal_ = *rho_;
    edm::Handle<edm::ValueMap<float>> full5x5SigmaIEtaIEtaMap;
    iEvent.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);
    edm::Handle<edm::ValueMap<float>> phoChargedIsolationMap;
    iEvent.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
    edm::Handle<edm::ValueMap<float>> phoNeutralHadronIsolationMap;
    iEvent.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
    edm::Handle<edm::ValueMap<float>> phoPhotonIsolationMap;
    iEvent.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);

	int cachecount1 = 0;
	int cachecount2 = 0;	

    for (size_t ip = 0; ip < photons->size(); ip++) {
        const auto pho = photons->ptrAt(ip);

        double phosc_eta = pho->superCluster()->eta();
        double phosc_phi = pho->superCluster()->phi();


         double pho_ieie = (*photons)[ip].full5x5_sigmaIetaIeta();
         double chIso1 = (*photons)[ip].userFloat("phoChargedIsolation");
         double nhIso1 = (*photons)[ip].userFloat("phoNeutralHadronIsolation");
         double phIso1 = (*photons)[ip].userFloat("phoPhotonIsolation");

		double chiso=std::max(0.0, chIso1 - rhoVal_*EAch(fabs((*photons)[ip].eta())));
		double nhiso=std::max(0.0, nhIso1 - rhoVal_*EAnh(fabs((*photons)[ip].eta())));
		double phoiso=std::max(0.0, phIso1 - rhoVal_*EApho(fabs((*photons)[ip].eta())));

        int                                   ismedium_photon   = 0;
        int                                   ismedium_photon_f = 0;
        edm::Handle<edm::View<pat::Electron>> electrons;
        iEvent.getByToken(electronToken_, electrons);
        edm::Handle<reco::BeamSpot> beamSpot;
        iEvent.getByToken(beamSpotToken_, beamSpot);
        edm::Handle<std::vector<reco::Conversion>> conversions;
        iEvent.getByToken(conversionsToken_, conversions);
        passEleVeto       = (!hasMatchedPromptElectron((*photons)[ip].superCluster(), electrons, conversions, beamSpot->position()));
        passEleVetonew    = (*photons)[ip].passElectronVeto();
        passPixelSeedVeto = (*photons)[ip].hasPixelSeed();

        if (ip < 6) {
            photon_pt[ip]     = (*photons)[ip].pt();
            photon_eta[ip]    = (*photons)[ip].eta();
            photon_phi[ip]    = (*photons)[ip].phi();
            photon_e[ip]      = (*photons)[ip].energy();
            photonsc_eta[ip]  = phosc_eta;
            photonsc_phi[ip]  = phosc_phi;
            photon_pev[ip]    = passEleVeto;
            photon_pevnew[ip] = passEleVetonew;
            photon_ppsv[ip]   = passPixelSeedVeto;
            photon_iseb[ip]   = (*photons)[ip].isEB();
            photon_isee[ip]   = (*photons)[ip].isEE();
            photon_hoe[ip]    = (*photons)[ip].hadTowOverEm();
            photon_sieie[ip]  = pho_ieie;  //(*photons)[ip].sigmaIetaIeta();
            photon_sieie2[ip] = (*photons)[ip].sigmaIetaIeta();
            photon_chiso[ip]  = chiso;
            photon_nhiso[ip]  = nhiso;
            photon_phoiso[ip] = phoiso;
            if (RunOnMC_ && photon_pt[ip] > 10) {
                const auto pho    = photons->ptrAt(ip);
                photon_istrue[ip] = matchToTruth(*pho, genParticles, ISRPho, dR_, photon_isprompt[ip]);
            }
            photon_drla[ip] = deltaR(photonsc_eta[ip], photonsc_phi[ip], etalep1, philep1);
            TLorentzVector tp4;
            tp4.SetPtEtaPhiE(photon_pt[ip], photon_eta[ip], photon_phi[ip], photon_e[ip]);
            photon_mla[ip] = (tp4 + glepton).M();
            TLorentzVector fwp4;
            fwp4.SetPtEtaPhiE(ptVlepJEC, yVlepJEC, phiVlepJEC, energyVlepJEC);
            //fwp4.SetPtEtaPhiE(WLeptonic.pt(), WLeptonic.eta(), WLeptonic.phi(), WLeptonic.energy());
            photon_mva[ip] = (tp4 + fwp4).M();
        }

		// photon ID
        if (fabs(phosc_eta) < 1.4442 && (*photons)[ip].hadTowOverEm() < 0.02197 && photon_sieie[ip] < 0.01015 && chiso < 1.141 && nhiso < (1.189 + (0.01512 * (*photons)[ip].pt() + 0.00002259 * (*photons)[ip].pt() * (*photons)[ip].pt())) && phoiso < (2.08 + 0.004017 * (*photons)[ip].pt())) {ismedium_photon = 1;}
        if (fabs(phosc_eta) > 1.566 && fabs(phosc_eta) < 2.5 && (*photons)[ip].hadTowOverEm() < 0.0326 && photon_sieie[ip] < 0.0272 && chiso < 1.051 && nhiso < (2.718 + (0.0117 * (*photons)[ip].pt() + 0.000023 * (*photons)[ip].pt() * (*photons)[ip].pt())) && phoiso < (3.867 + 0.0037 * (*photons)[ip].pt())) {ismedium_photon = 1;}

        if (ismedium_photon == 1 && deltaR(phosc_eta, phosc_phi, etalep1, philep1) > 0.5) {
            if (cachecount1 == 0) {
                photonet = (*photons)[ip].pt();
                iphoton  = ip;
            }
            cachecount1++;
            if ((*photons)[ip].pt() > photonet) {
                photonet = (*photons)[ip].pt();
                iphoton  = ip;
            }
        }
	
	// fake photon
        //Inverting loose ID
        //            if(passEleVetonew && (*photons)[ip].isEB() && (*photons)[ip].hadTowOverEm()<5*0.0597 && chiso<10 && nhiso<std::min(0.2*(*photons)[ip].pt(), 5.*(10.910 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt())))  && phoiso<std::min(0.2*(*photons)[ip].pt(), 5.*(3.630+0.0047*(*photons)[ip].pt())) && !(photon_sieie[ip]<0.01031 && chiso<4)) {ismedium_photon_f=1;}  // && nhiso<(10.910 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt())) && phoiso<(3.630+0.0047*(*photons)[ip].pt())
        //            if(passEleVetonew && (*photons)[ip].isEE() && (*photons)[ip].hadTowOverEm()<5*0.0481 && chiso<10 && nhiso<std::min(0.2*(*photons)[ip].pt(), 5.*(5.931 + (0.0163*(*photons)[ip].pt()+0.000014*(*photons)[ip].pt()*(*photons)[ip].pt())))  && phoiso<std::min(0.2*(*photons)[ip].pt(), 5.*(6.641+0.0034*(*photons)[ip].pt())) && !(photon_sieie[ip]<0.03013 && chiso<4)) {ismedium_photon_f=1;}  // && nhiso<(5.931 + (0.0163*(*photons)[ip].pt()+0.000014*(*photons)[ip].pt()*(*photons)[ip].pt())) && phoiso<(6.641+0.0034*(*photons)[ip].pt())

        //            if(passEleVetonew && phosc_eta<1.4442 && (*photons)[ip].hadTowOverEm()<0.0597 && chiso<10 && nhiso<std::min(0.2*(*photons)[ip].pt(), (10.910 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt())))  && phoiso<std::min(0.2*(*photons)[ip].pt(), (3.630+0.0047*(*photons)[ip].pt())) && !(photon_sieie[ip]<0.01031 && chiso<4)) {ismedium_photon_f=1;}
        //            if(passEleVetonew && phosc_eta>1.566 && phosc_eta<2.5 && (*photons)[ip].hadTowOverEm()<0.0481 && chiso<10 && nhiso<std::min(0.2*(*photons)[ip].pt(), (5.931 + (0.0163*(*photons)[ip].pt()+0.000014*(*photons)[ip].pt()*(*photons)[ip].pt())))  && phoiso<std::min(0.2*(*photons)[ip].pt(), (6.641+0.0034*(*photons)[ip].pt())) && !(photon_sieie[ip]<0.03013 && chiso<4)) {ismedium_photon_f=1;}

        //            if(passEleVetonew && phosc_eta<1.4442 && (*photons)[ip].hadTowOverEm()<0.0396 && chiso<10 && nhiso<(2.725 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt()))  && phoiso<(2.571+0.0047*(*photons)[ip].pt()) && !(photon_sieie[ip]<0.01022 && chiso<4)) {ismedium_photon_f=1;}

        if (fabs(phosc_eta) < 1.4442 && !((*photons)[ip].hadTowOverEm() < 0.02197 && photon_sieie[ip] < 0.01015 && chiso < 1.141 && nhiso < (1.189 + (0.01512 * (*photons)[ip].pt() + 0.00002259 * (*photons)[ip].pt() * (*photons)[ip].pt())) && phoiso < (2.08 + 0.004017 * (*photons)[ip].pt()))) {
            ismedium_photon_f = 1;
        }
        if (fabs(phosc_eta) > 1.566 && fabs(phosc_eta) < 2.5 && !((*photons)[ip].hadTowOverEm() < 0.0326 && photon_sieie[ip] < 0.0272 && chiso < 1.051 && nhiso < (2.718 + (0.0117 * (*photons)[ip].pt() + 0.000023 * (*photons)[ip].pt() * (*photons)[ip].pt())) && phoiso < (3.867 + 0.0037 * (*photons)[ip].pt()))) {
            ismedium_photon_f = 1;
        }

        if (ismedium_photon_f == 1 && deltaR(phosc_eta, phosc_phi, etalep1, philep1) > 0.5) {
            if (cachecount2 == 0) {
                photonet_f = (*photons)[ip].pt();
                iphoton_f  = ip;
            }
            cachecount2++;
            if ((*photons)[ip].pt() > photonet_f) {
                photonet_f = (*photons)[ip].pt();
                iphoton_f  = ip;
            }
        }
	}

	// gen photon match
    if (RunOnMC_ && iphoton > -1) {
        const auto pho1 = photons->ptrAt(iphoton);
        isTrue_         = matchToTruth(*pho1, genParticles, ISRPho, dR_, isprompt_);
    }

    if (RunOnMC_ && iphoton > -1) {
        const auto pho1 = photons->ptrAt(iphoton);
        isTrue_         = matchToTruth(*pho1, genParticles, ISRPho, dR_, isprompt_);
    }

    if (iphoton > -1 && iphoton < 6) {
        photonet     = photon_pt[iphoton];   //(*photons)[iphoton].pt();
        photoneta    = photon_eta[iphoton];  //(*photons)[iphoton].eta();
        photonphi    = photon_phi[iphoton];  //(*photons)[iphoton].phi();
        photone      = photon_e[iphoton];    //(*photons)[iphoton].energy();
        photonsceta  = photonsc_eta[iphoton];
        photonscphi  = photonsc_phi[iphoton];   
        photonsieie  = photon_sieie[iphoton];   //(*photons)[iphoton].sigmaIetaIeta();
        photonphoiso = photon_phoiso[iphoton];  //std::max((*photons)[iphoton].photonIso()-rhoVal_*EApho(fabs((*photons)[iphoton].eta())),0.0);
        photonchiso  = photon_chiso[iphoton];   //std::max((*photons)[iphoton].chargedHadronIso()-rhoVal_*EAch(fabs((*photons)[iphoton].eta())),0.0);
        photonnhiso  = photon_nhiso[iphoton];   //std::max((*photons)[iphoton].neutralHadronIso()-rhoVal_*EAnh(fabs((*photons)[iphoton].eta())),0.0);
        drla         = deltaR(photonsc_eta[iphoton], photonsc_phi[iphoton], etalep1, philep1);
        TLorentzVector photonp4;
        photonp4.SetPtEtaPhiE(photonet, photoneta, photonphi, photone);
        Mla = (photonp4 + glepton).M();
        TLorentzVector wp4;
       	//wp4.SetPtEtaPhiE(WLeptonic.pt(), WLeptonic.eta(), WLeptonic.phi(), WLeptonic.energy());
		wp4.SetPtEtaPhiE(ptVlepJEC, yVlepJEC, phiVlepJEC, energyVlepJEC);
        Mva                = (photonp4 + wp4).M();
        photonhaspixelseed = photon_ppsv[iphoton];
        photonpasseleveto  = photon_pevnew[iphoton];
    }

    if (iphoton_f > -1 && iphoton_f < 6) {
        photonet_f     = photon_pt[iphoton_f];   //(*photons)[iphoton_f].pt();
        photoneta_f    = photon_eta[iphoton_f];  //(*photons)[iphoton_f].eta();
        photonphi_f    = photon_phi[iphoton_f];  //(*photons)[iphoton_f].phi();
        photone_f      = photon_e[iphoton_f];    //(*photons)[iphoton_f].energy();
        photonsceta_f  = photonsc_eta[iphoton_f];
        photonscphi_f  = photonsc_phi[iphoton_f];
        photonsieie_f  = photon_sieie[iphoton_f];   //(*photons)[iphoton_f].sigmaIetaIeta();
        photonphoiso_f = photon_phoiso[iphoton_f];  //std::max((*photons)[iphoton_f].photonIso()-rhoVal_*EApho(fabs((*photons)[iphoton_f].eta())),0.0);
        photonchiso_f  = photon_chiso[iphoton_f];   //std::max((*photons)[iphoton_f].chargedHadronIso()-rhoVal_*EAch(fabs((*photons)[iphoton_f].eta())),0.0);
        photonnhiso_f  = photon_nhiso[iphoton_f];   //std::max((*photons)[iphoton_f].neutralHadronIso()-rhoVal_*EAnh(fabs((*photons)[iphoton_f].eta())),0.0);
        drla_f         = deltaR(photonsc_eta[iphoton_f], photonsc_phi[iphoton_f], etalep1, philep1);
        TLorentzVector photonp4_f;
        photonp4_f.SetPtEtaPhiE(photonet_f, photoneta_f, photonphi_f, photone_f);
        Mla_f = (photonp4_f + glepton).M();
        TLorentzVector wp4_f;
        //wp4_f.SetPtEtaPhiE(WLeptonic.pt(), WLeptonic.eta(), WLeptonic.phi(), WLeptonic.energy());
        wp4_f.SetPtEtaPhiE(ptVlepJEC, yVlepJEC, phiVlepJEC, energyVlepJEC);
        Mva_f                = (photonp4_f + wp4_f).M();
        photonhaspixelseed_f = photon_ppsv[iphoton_f];
        photonpasseleveto_f  = photon_pevnew[iphoton_f];
    }


}



#endif
