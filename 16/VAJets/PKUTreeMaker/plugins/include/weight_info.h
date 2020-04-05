#ifndef _Weight_
#define _Weight_

void PKUTreeMaker::weight_info(edm::Event const & iEvent) {

	edm::Handle<GenEventInfoProduct> genEvtInfo;
	iEvent.getByToken(GenToken_, genEvtInfo);
        
	theWeight = genEvtInfo->weight();
        
	if (theWeight > 0)
		nump = nump + 1;
        
	if (theWeight < 0)
		numm = numm + 1;
        
	edm::Handle<std::vector<PileupSummaryInfo>> PupInfo;
	iEvent.getByToken(PUToken_, PupInfo);
        
	std::vector<PileupSummaryInfo>::const_iterator PVI;
        
	for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
		nBX = PVI->getBunchCrossing();
            
		if (nBX == 0) {
			npT  = PVI->getTrueNumInteractions();
			npIT = PVI->getPU_NumInteractions();
		}
	}
        
	edm::Handle< double > theprefweight;
	iEvent.getByToken(prefweight_token, theprefweight ) ;
	_prefiringweight =(*theprefweight);

	edm::Handle< double > theprefweightup;
	iEvent.getByToken(prefweightup_token, theprefweightup ) ;
	_prefiringweightup =(*theprefweightup);

 	edm::Handle< double > theprefweightdown;
	iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
	_prefiringweightdown =(*theprefweightdown);
	

}

#endif
