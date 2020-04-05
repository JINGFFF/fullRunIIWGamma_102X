#ifndef _MET_
#define _MET_

void PKUTreeMaker::met_info(edm::Event const & iEvent) {
    edm::Handle<pat::METCollection> METs_;
    bool defaultMET = iEvent.getByToken(metInputToken_ , METs_ );

    if(RunOnMC_){
        const pat::MET &xmet = METs_->front();
        genMET=xmet.genMET()->pt();
    }

    if(defaultMET){
        addTypeICorr(iEvent);
        addTypeICorr_user(iEvent);
        for (const pat::MET &met : *METs_) {

            const float rawPt = met.uncorPt();
            const float rawPhi = met.uncorPhi();
            const float rawSumEt = met.uncorSumEt();
            TVector2 rawMET_;
            rawMET_.SetMagPhi (rawPt, rawPhi );
            Double_t rawPx = rawMET_.Px();
            Double_t rawPy = rawMET_.Py();
            Double_t rawEt = std::hypot(rawPx,rawPy);
            METraw_et = rawEt;
            METraw_phi = rawPhi;
            METraw_sumEt = rawSumEt;

            double pxcorr = rawPx+TypeICorrMap_["corrEx"];
            double pycorr = rawPy+TypeICorrMap_["corrEy"];
            double et     = std::hypot(pxcorr,pycorr);
            double sumEtcorr = rawSumEt+TypeICorrMap_["corrSumEt"];

            // Marked for debug
            //------------------central value, correction from JetuserData---------------------
            double pxcorr_new= rawPx+TypeICorrMap_user_["corrEx_JEC"]+TypeICorrMap_user_["corrEx_JER"];
            double pycorr_new= rawPy+TypeICorrMap_user_["corrEy_JEC"]+TypeICorrMap_user_["corrEy_JER"];
            double et_new     = std::hypot(pxcorr_new,pycorr_new);
            double sumEtcorr_new = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC"]+TypeICorrMap_user_["corrSumEt_JER"];

            //----for JEC uncertainty study
            double pxcorr_JEC_up = rawPx+TypeICorrMap_user_["corrEx_JEC_up"]+TypeICorrMap_user_["corrEx_JER"];
            double pycorr_JEC_up = rawPy+TypeICorrMap_user_["corrEy_JEC_up"]+TypeICorrMap_user_["corrEy_JER"];
            double et_JEC_up     = std::hypot(pxcorr_JEC_up, pycorr_JEC_up);
            double sumEtcorr_JEC_up = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC_up"]+TypeICorrMap_user_["corrSumEt_JER"];
            double pxcorr_JEC_down = rawPx+TypeICorrMap_user_["corrEx_JEC_down"]+TypeICorrMap_user_["corrEx_JER"];
            double pycorr_JEC_down = rawPy+TypeICorrMap_user_["corrEy_JEC_down"]+TypeICorrMap_user_["corrEy_JER"];
            double et_JEC_down     = std::hypot(pxcorr_JEC_down, pycorr_JEC_down);
            double sumEtcorr_JEC_down = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC_down"]+TypeICorrMap_user_["corrSumEt_JER"];

            //----for JER uncertainty study
            double pxcorr_JER_up = rawPx+TypeICorrMap_user_["corrEx_JEC"]+TypeICorrMap_user_["corrEx_JER_up"];
            double pycorr_JER_up = rawPy+TypeICorrMap_user_["corrEy_JEC"]+TypeICorrMap_user_["corrEy_JER_up"];
            double et_JER_up     = std::hypot(pxcorr_JER_up, pycorr_JER_up);
            double sumEtcorr_JER_up = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC"]+TypeICorrMap_user_["corrSumEt_JER_up"];
            double pxcorr_JER_down = rawPx+TypeICorrMap_user_["corrEx_JEC"]+TypeICorrMap_user_["corrEx_JER_down"];
            double pycorr_JER_down = rawPy+TypeICorrMap_user_["corrEy_JEC"]+TypeICorrMap_user_["corrEy_JER_down"];
            double et_JER_down     = std::hypot(pxcorr_JER_down,pycorr_JER_down);
            double sumEtcorr_JER_down = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC"]+TypeICorrMap_user_["corrSumEt_JER_down"];

            //------------------ correction from JetuserData---------------------
            // Marked for debug
            TLorentzVector corrmet;

            corrmet.SetPxPyPzE(pxcorr,pycorr,0.,et);
            MET_et = et;
            MET_phi = corrmet.Phi();
            MET_sumEt = sumEtcorr;
            useless = sumEtcorr;
            useless = rawEt;
            MET_corrPx = TypeICorrMap_["corrEx"];
            MET_corrPy = TypeICorrMap_["corrEy"];

            // Marked for debug
            MET_et_new= et_new;
            MET_et_JEC_up = et_JEC_up;
            MET_et_JEC_down = et_JEC_down;
            MET_et_JER_up = et_JER_up;
            MET_et_JER_down = et_JER_down;

            corrmet.SetPxPyPzE(pxcorr_new,pycorr_new,0.,et_new);
            MET_phi_new = corrmet.Phi();
            corrmet.SetPxPyPzE(pxcorr_JEC_up,pycorr_JEC_up,0.,et_JEC_up);
            MET_phi_JEC_up = corrmet.Phi();
            corrmet.SetPxPyPzE(pxcorr_JEC_down,pycorr_JEC_down,0.,et_JEC_down);
            MET_phi_JEC_down = corrmet.Phi();
            corrmet.SetPxPyPzE(pxcorr_JER_up,pycorr_JER_up,0.,et_JER_up);
            MET_phi_JER_up = corrmet.Phi();
            corrmet.SetPxPyPzE(pxcorr_JER_down,pycorr_JER_down,0.,et_JER_down);
            MET_phi_JER_down = corrmet.Phi();

            MET_sumEt_new = sumEtcorr_new;
            MET_sumEt_JEC_up = sumEtcorr_JEC_up;
            MET_sumEt_JEC_down = sumEtcorr_JEC_down;
            MET_sumEt_JER_up = sumEtcorr_JER_up;
            MET_sumEt_JER_down = sumEtcorr_JER_down;
            // Marked for debug
		}
	}	

}

#endif
