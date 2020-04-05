#ifndef _HLT_
#define _HLT_

void PKUTreeMaker::hlt_info(edm::Event const & iEvent) {
	
    Handle<TriggerResults> trigRes;
    iEvent.getByToken(hltToken_, trigRes);

    int xtemp1 = 0;
    for (size_t i = 0; i < elPaths1.size(); i++) {
        xtemp1 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths1[i]));
        if (HLT_Ele1 < xtemp1)
            HLT_Ele1 = xtemp1;
    }

    int xtemp2 = 0;
    for (size_t i = 0; i < elPaths2.size(); i++) {
        xtemp2 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths2[i]));
        if (HLT_Ele2 < xtemp2)
            HLT_Ele2 = xtemp2;
    }

    int mtemp1 = 0;
    for (size_t i = 0; i < muPaths1.size(); i++) {
        mtemp1 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths1[i]));
        if (HLT_Mu1 < mtemp1)
            HLT_Mu1 = mtemp1;
    }

    int mtemp2 = 0;
    for (size_t i = 0; i < muPaths2.size(); i++) {
        mtemp2 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths2[i]));
        if (HLT_Mu2 < mtemp2)
            HLT_Mu2 = mtemp2;
    }

    int mtemp3 = 0;
    for (size_t i = 0; i < muPaths3.size(); i++) {
        mtemp3 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths3[i]));
        if (HLT_Mu3 < mtemp3)
            HLT_Mu3 = mtemp3;
    }

}

#endif
