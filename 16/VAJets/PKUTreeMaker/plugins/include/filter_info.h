#ifndef _Filter_
#define _Filter_

void PKUTreeMaker::filter_info(edm::Event const & iEvent) {
	iEvent.getByToken(noiseFilterToken_, noiseFilterBits_);
	const edm::TriggerNames& names = iEvent.triggerNames(*noiseFilterBits_);
    
	for (unsigned int i = 0, n = noiseFilterBits_->size(); i < n; ++i) {
		if (names.triggerName(i) == HBHENoiseFilter_Selector_)
			passFilter_HBHE_ = noiseFilterBits_->accept(i);
        
		if (names.triggerName(i) == HBHENoiseIsoFilter_Selector_)
			passFilter_HBHEIso_ = noiseFilterBits_->accept(i);
        
		if (names.triggerName(i) == ECALDeadCellNoiseFilter_Selector_)
			passFilter_globalTightHalo_ = noiseFilterBits_->accept(i);
        
		if (names.triggerName(i) == ECALDeadCellNoiseFilter_Selector_)
			passFilter_ECALDeadCell_ = noiseFilterBits_->accept(i);
        
		if (names.triggerName(i) == GoodVtxNoiseFilter_Selector_)
			passFilter_GoodVtx_ = noiseFilterBits_->accept(i);
        
		if (names.triggerName(i) == EEBadScNoiseFilter_Selector_)
			passFilter_EEBadSc_ = noiseFilterBits_->accept(i);
    }
	
	edm::Handle<bool> badMuonResultHandle;
	edm::Handle<bool> badChargedHadronResultHandle;
    
	iEvent.getByToken(badMuon_Selector_, badMuonResultHandle);
	iEvent.getByToken(badChargedHadron_Selector_, badChargedHadronResultHandle);
	
	passFilter_badMuon_          = *badMuonResultHandle;
	passFilter_badChargedHadron_ = *badChargedHadronResultHandle;	

}

#endif
