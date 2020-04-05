#ifndef _Gen_
#define _Gen_

void PKUTreeMaker::gen_photon_lepton_info(edm::Event const & iEvent) {
	//define genParticle
	edm::Handle<edm::View<reco::GenParticle>> genParticles;  
	iEvent.getByToken(genSrc_, genParticles);

	int ipp = 0, imm = 0, iee = 0;
        
	for (size_t i = 0; i < genParticles->size(); i++) {
		const reco::Candidate* particle = &(*genParticles)[i];
            
		if (abs(particle->pdgId()) == 22 && particle->status() == 1 && (*genParticles)[i].isPromptFinalState() > 0 && ipp < 6) {
			genphoton_pt[ipp]  = particle->pt();
			genphoton_eta[ipp] = particle->eta();
			genphoton_phi[ipp] = particle->phi();
			ipp++;  
            
		}
            
		if (abs(particle->pdgId()) == 13 && particle->status() == 1 && (*genParticles)[i].isPromptFinalState() > 0 && imm < 6) {
			genmuon_pt[imm]  = particle->pt();
			genmuon_eta[imm] = particle->eta();
			genmuon_phi[imm] = particle->phi();
			imm++;
           
		}
            
		if (abs(particle->pdgId()) == 11 && particle->status() == 1 && (*genParticles)[i].isPromptFinalState() > 0 && iee < 6) {
			genelectron_pt[iee]  = particle->pt();
			genelectron_eta[iee] = particle->eta();
			genelectron_phi[iee] = particle->phi();
			iee++;
            
		}
	}	

}

void PKUTreeMaker::gen_jet_info(edm::Event const & iEvent) {
        
	int ijj=0;
	edm::Handle<reco::GenJetCollection> genJets;
	iEvent.getByToken(genJet_,genJets);
	reco::GenJetCollection::const_iterator i_jet;
        
	for( i_jet=genJets->begin(); i_jet != genJets->end();i_jet++){
		genjet_e[ijj] = i_jet->energy();
		genjet_pt[ijj]= i_jet->pt();
		genjet_eta[ijj]= i_jet->eta();
		genjet_phi[ijj]=i_jet->phi();
		ijj++;
 
	}
}

#endif
