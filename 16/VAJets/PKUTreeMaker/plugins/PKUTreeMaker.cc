// system include files

#include "TMath.h"
#include <iostream>
#include <memory>
// user include files
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include <algorithm>
#define Pi 3.141593
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "Math/VectorUtil.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TMath.h"
#include <TFormula.h>

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "include/PKUTreeMaker.h"
#include "include/setDummyValues.h"
#include "include/function.h"
using namespace std;

       
//------------------------------------
//------------------------------------
void PKUTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    setDummyValues();  //Initalize variables with dummy values
    nevent = iEvent.eventAuxiliary().event();
    run    = iEvent.eventAuxiliary().run();
    ls     = iEvent.eventAuxiliary().luminosityBlock();

    //events weight
    if (RunOnMC_) {
        //        std::cout<<lheEvtInfo->hepeup().NUP<<std::endl;
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
            if (nBX == 0) {  // "0" is the in-time crossing, negative values are the early crossings, positive are late
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
        //std::cout<<muPaths1.size()<<" "<<muPaths1[i]<<std::endl;
        //std::cout<<(hltConfig.triggerIndex(muPaths1[i]))<<std::endl;
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


   edm::Handle<edm::View<reco::Candidate> > leptonicVs;
   iEvent.getByToken(leptonicVSrc_, leptonicVs);

    if (leptonicVs->empty()) {
        outTree_->Fill();
        return;
    }  //outTree_->Fill();

    iEvent.getByToken(rhoToken_, rho_);
    double fastJetRho = *(rho_.product());
    useless           = fastJetRho;

    edm::Handle<edm::View<pat::Jet>> ak4jets;
    iEvent.getByToken(ak4jetsSrc_, ak4jets);

    edm::Handle<edm::View<pat::Photon>> photons;
    iEvent.getByToken(photonSrc_, photons);

    //if (photons->empty()) {
    //    outTree_->Fill();
    //    return;
    //}  //outTree_->Fill();

   if (photons->empty()) {   hasphoton = 0.;  }
   else {hasphoton =1.;}//outTree_->Fill();

    edm::Handle<edm::View<reco::GenParticle>> genParticles;  //define genParticle
    iEvent.getByToken(genSrc_, genParticles);
    //   iEvent.getByLabel(InputTag("packedGenParticles"), genParticles);

    if (RunOnMC_) {
        int ipp = 0, imm = 0, iee = 0;
        for (size_t i = 0; i < genParticles->size(); i++) {  // std::cout<<"i = "<<i<<std::endl;
            const reco::Candidate* particle = &(*genParticles)[i];
            if (abs(particle->pdgId()) == 22 && particle->status() == 1 && (*genParticles)[i].isPromptFinalState() > 0 && ipp < 6) {
                genphoton_pt[ipp]  = particle->pt();
                genphoton_eta[ipp] = particle->eta();
                genphoton_phi[ipp] = particle->phi();
                ipp++;  //	std::cout<<"ipp = "<<ipp<<std::endl;
            }           //std::cout<<"ipp = "<<ipp<<std::endl;
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

	if(RunOnMC_){
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

    //filter
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

    const reco::Candidate& leptonicV = leptonicVs->at(0);
    const reco::Candidate& metCand   = metHandle->at(0);
    const reco::Candidate& lepton    = (*leptonicV.daughter(0));

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(VertexToken_, vertices);
    if (vertices->empty()) {
        outTree_->Fill();
        return;
    }  // skip the event if no PV foundoutTree_->Fill();
    nVtx                                                   = vertices->size();
    reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx) {
        // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
        // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
        if (  // !vtx->isFake() &&
            !(vtx->chi2() == 0 && vtx->ndof() == 0)
            && vtx->ndof() >= 4. && vtx->position().Rho() <= 2.0
            && fabs(vtx->position().Z()) <= 24.0) {
            firstGoodVertex = vtx;
            break;
        }
    }
    if (firstGoodVertex == vertices->end()) {
        outTree_->Fill();
        return;
    }  // skip event if there are no good PVsoutTree_->Fill();

    // ************************* MET ********************** //
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
			//         const float  rawPt    = met.shiftedPt(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
			//         const float  rawPhi   = met.shiftedPhi(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
			//         const float  rawSumEt = met.shiftedSumEt(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);


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
    //------------------------------------
    /// For the time being, set these to 1
    triggerWeight       = 1.0;
    pileupWeight        = 1.0;
    double targetEvents = targetLumiInvPb_ * crossSectionPb_;
    lumiWeight          = targetEvents / originalNEvents_;
    lep                 = abs(leptonicV.daughter(0)->pdgId());  //std::max(abs(leptonicV.daughter(0)->pdgId()), abs(leptonicV.daughter(1)->pdgId()));
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
    reco::CandidateBaseRef      METBaseRef = metHandle->refAt(0);  //?????
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
    if (RunOnMC_ && ptlep1 > 10) {
        //const auto lept = lepton;
        lepton_istrue = matchToTrueLep(etalep1, philep1, genParticles, dR1_, ispromptLep_);
    }

    // ************************* Photon Jets Information****************** //
    // *************************************************************//
    double rhoVal_;
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

    photonet        = -10.;
    photonet_f      = -10.;
    iphoton         = -1;
    iphoton_f       = -1;
    int cachecount1 = 0;  //added by Qianming Huang !!!
    int cachecount2 = 0;  //added by Qianming Huang !!!
    for (size_t ip = 0; ip < photons->size(); ip++) {
        const auto pho = photons->ptrAt(ip);

        double phosc_eta = pho->superCluster()->eta();
        double phosc_phi = pho->superCluster()->phi();


         double pho_ieie = (*photons)[ip].full5x5_sigmaIetaIeta();
         double chIso1 = (*photons)[ip].userFloat("phoChargedIsolation");
         double nhIso1 = (*photons)[ip].userFloat("phoNeutralHadronIsolation");
         double phIso1 = (*photons)[ip].userFloat("phoPhotonIsolation");
        
            double chiso=std::max(0.0, chIso1 - rhoVal_*EAch(fabs((*photons)[ip].eta()))); //effAreaChHadrons_.getEffectiveArea(fabs(phosc_eta)));
//            double chiso=std::max((*photons)[ip].chargedHadronIso()-rhoVal_*EAch(fabs((*photons)[ip].eta())),0.0);
            double nhiso=std::max(0.0, nhIso1 - rhoVal_*EAnh(fabs((*photons)[ip].eta()))); //effAreaNeuHadrons_.getEffectiveArea(fabs(phosc_eta)));
//            double nhiso=std::max((*photons)[ip].neutralHadronIso()-rhoVal_*EAnh(fabs((*photons)[ip].eta())),0.0);
            double phoiso=std::max(0.0, phIso1 - rhoVal_*EApho(fabs((*photons)[ip].eta()))); //effAreaPhotons_.getEffectiveArea(fabs(phosc_eta)));
//            double phoiso=std::max((*photons)[ip].photonIso()-rhoVal_*EApho(fabs((*photons)[ip].eta())),0.0);


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
            fwp4.SetPtEtaPhiE(WLeptonic.pt(), WLeptonic.eta(), WLeptonic.phi(), WLeptonic.energy());
            photon_mva[ip] = (tp4 + fwp4).M();
        }

        //if (fabs(phosc_eta) < 1.4442 && (*photons)[ip].hadTowOverEm() < 0.0396 && photon_sieie[ip] < 0.01022 && chiso < 0.441 && nhiso < (2.725 + (0.0148 * (*photons)[ip].pt() + 0.000017 * (*photons)[ip].pt() * (*photons)[ip].pt())) && phoiso < (2.571 + 0.0047 * (*photons)[ip].pt())) {ismedium_photon = 1;}
        //if (fabs(phosc_eta) > 1.566 && fabs(phosc_eta) < 2.5 && (*photons)[ip].hadTowOverEm() < 0.0219 && photon_sieie[ip] < 0.03001 && chiso < 0.442 && nhiso < (1.715 + (0.0163 * (*photons)[ip].pt() + 0.000014 * (*photons)[ip].pt() * (*photons)[ip].pt())) && phoiso < (3.863 + 0.0034 * (*photons)[ip].pt())) {ismedium_photon = 1;}

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
        //////////////////////////////////for fake photon study, store photon without sieie cut
        //Inverting loose ID
        //            if(passEleVetonew && (*photons)[ip].isEB() && (*photons)[ip].hadTowOverEm()<5*0.0597 && chiso<10 && nhiso<std::min(0.2*(*photons)[ip].pt(), 5.*(10.910 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt())))  && phoiso<std::min(0.2*(*photons)[ip].pt(), 5.*(3.630+0.0047*(*photons)[ip].pt())) && !(photon_sieie[ip]<0.01031 && chiso<4)) {ismedium_photon_f=1;}  // && nhiso<(10.910 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt())) && phoiso<(3.630+0.0047*(*photons)[ip].pt())
        //            if(passEleVetonew && (*photons)[ip].isEE() && (*photons)[ip].hadTowOverEm()<5*0.0481 && chiso<10 && nhiso<std::min(0.2*(*photons)[ip].pt(), 5.*(5.931 + (0.0163*(*photons)[ip].pt()+0.000014*(*photons)[ip].pt()*(*photons)[ip].pt())))  && phoiso<std::min(0.2*(*photons)[ip].pt(), 5.*(6.641+0.0034*(*photons)[ip].pt())) && !(photon_sieie[ip]<0.03013 && chiso<4)) {ismedium_photon_f=1;}  // && nhiso<(5.931 + (0.0163*(*photons)[ip].pt()+0.000014*(*photons)[ip].pt()*(*photons)[ip].pt())) && phoiso<(6.641+0.0034*(*photons)[ip].pt())

        //            if(passEleVetonew && phosc_eta<1.4442 && (*photons)[ip].hadTowOverEm()<0.0597 && chiso<10 && nhiso<std::min(0.2*(*photons)[ip].pt(), (10.910 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt())))  && phoiso<std::min(0.2*(*photons)[ip].pt(), (3.630+0.0047*(*photons)[ip].pt())) && !(photon_sieie[ip]<0.01031 && chiso<4)) {ismedium_photon_f=1;}
        //            if(passEleVetonew && phosc_eta>1.566 && phosc_eta<2.5 && (*photons)[ip].hadTowOverEm()<0.0481 && chiso<10 && nhiso<std::min(0.2*(*photons)[ip].pt(), (5.931 + (0.0163*(*photons)[ip].pt()+0.000014*(*photons)[ip].pt()*(*photons)[ip].pt())))  && phoiso<std::min(0.2*(*photons)[ip].pt(), (6.641+0.0034*(*photons)[ip].pt())) && !(photon_sieie[ip]<0.03013 && chiso<4)) {ismedium_photon_f=1;}

        //            if(passEleVetonew && phosc_eta<1.4442 && (*photons)[ip].hadTowOverEm()<0.0396 && chiso<10 && nhiso<(2.725 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt()))  && phoiso<(2.571+0.0047*(*photons)[ip].pt()) && !(photon_sieie[ip]<0.01022 && chiso<4)) {ismedium_photon_f=1;}
        //            if(passEleVetonew && phosc_eta>1.566 && phosc_eta<2.5 && (*photons)[ip].hadTowOverEm()<0.0219 && chiso<10 && nhiso<(1.715 + (0.0163*(*photons)[ip].pt()+0.000014*(*photons)[ip].pt()*(*photons)[ip].pt()))  && phoiso<(3.863+0.0034*(*photons)[ip].pt()) && !(photon_sieie[ip]<0.03001 && chiso<4)) {ismedium_photon_f=1;}

        if (fabs(phosc_eta) < 1.4442 && !((*photons)[ip].hadTowOverEm() < 0.02197 && photon_sieie[ip] < 0.01015 && chiso < 1.141 && nhiso < (1.189 + (0.01512 * (*photons)[ip].pt() + 0.00002259 * (*photons)[ip].pt() * (*photons)[ip].pt())) && phoiso < (2.08 + 0.004017 * (*photons)[ip].pt()))) {
            ismedium_photon_f = 1;
        }
        if (fabs(phosc_eta) > 1.566 && fabs(phosc_eta) < 2.5 && !((*photons)[ip].hadTowOverEm() < 0.0326 && photon_sieie[ip] < 0.0272 && chiso < 1.051 && nhiso < (2.718 + (0.0117 * (*photons)[ip].pt() + 0.000023 * (*photons)[ip].pt() * (*photons)[ip].pt())) && phoiso < (3.867 + 0.0037 * (*photons)[ip].pt()))) {
            ismedium_photon_f = 1;
        }

        //            if(phosc_eta<1.4442 && (*photons)[ip].hadTowOverEm()<0.0597 && chiso<15 && (nhiso<std::min(0.2*(*photons)[ip].pt(), 5.*(10.910 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt())))  || phoiso<std::min(0.2*(*photons)[ip].pt(), 5.*(3.630+0.0047*(*photons)[ip].pt())) || !(photon_sieie[ip]<0.01031 && chiso<4))) {ismedium_photon_f=1;} // && nhiso<(10.910 + (0.0148*(*photons)[ip].pt()+0.000017*(*photons)[ip].pt()*(*photons)[ip].pt())) && phoiso<(3.630+0.0047*(*photons)[ip].pt())
        //            if(phosc_eta>1.566 && phosc_eta<2.5  && (*photons)[ip].hadTowOverEm()<0.0481 && chiso<15 && (nhiso<std::min(0.2*(*photons)[ip].pt(), 5.*(5.931 + (0.0163*(*photons)[ip].pt()+0.000014*(*photons)[ip].pt()*(*photons)[ip].pt())))  || phoiso<std::min(0.2*(*photons)[ip].pt(), 5.*(6.641+0.0034*(*photons)[ip].pt())) || !(photon_sieie[ip]<0.03013 && chiso<4))) {ismedium_photon_f=1;}  // && nhiso<(5.931 + (0.0163*(*photons)[ip].pt()+0.000014*(*photons)[ip].pt()*(*photons)[ip].pt())) && phoiso<(6.641+0.0034*(*photons)[ip].pt())

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
        ////////////////////////////////////////////////////////////////////////////////////////
    }

    //Gen photon matching
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
        wp4.SetPtEtaPhiE(WLeptonic.pt(), WLeptonic.eta(), WLeptonic.phi(), WLeptonic.energy());
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
        wp4_f.SetPtEtaPhiE(WLeptonic.pt(), WLeptonic.eta(), WLeptonic.phi(), WLeptonic.energy());
        Mva_f                = (photonp4_f + wp4_f).M();
        photonhaspixelseed_f = photon_ppsv[iphoton_f];
        photonpasseleveto_f  = photon_pevnew[iphoton_f];
    }

    // ************************* AK4 Jets Information****************** //
    // ***********************************************************//
    Int_t jetindexphoton12[2]            = {-1, -1};
    Int_t jetindexphoton12_f[2]          = {-1, -1};
    Int_t jetindexphoton12_new[2]        = {-1, -1};
    Int_t jetindexphoton12_new_f[2]      = {-1, -1};
    Int_t jetindexphoton12_JEC_up[2]     = {-1, -1};
    Int_t jetindexphoton12_JEC_up_f[2]   = {-1, -1};
    Int_t jetindexphoton12_JEC_down[2]   = {-1, -1};
    Int_t jetindexphoton12_JEC_down_f[2] = {-1, -1};
    Int_t jetindexphoton12_JER_up[2]     = {-1, -1};
    Int_t jetindexphoton12_JER_up_f[2]   = {-1, -1};
    Int_t jetindexphoton12_JER_down[2]   = {-1, -1};
    Int_t jetindexphoton12_JER_down_f[2] = {-1, -1};

    std::vector<JetCorrectorParameters> vPar;
    for (std::vector<std::string>::const_iterator payloadBegin = jecAK4Labels_.begin(), payloadEnd = jecAK4Labels_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    jecAK4_ = new FactorizedJetCorrector(vPar);
    vPar.clear();

    int                          nujets      = 0;
    double                       tmpjetptcut = 20.0;
    std::vector<TLorentzVector*> jets;


    //################Jet Correction##########################
    //two leading jets without JER
    for (size_t ik = 0; ik < ak4jets->size(); ik++) {
        reco::Candidate::LorentzVector uncorrJet = (*ak4jets)[ik].correctedP4(0);
        jecAK4_->setJetEta(uncorrJet.eta());
        jecAK4_->setJetPt(uncorrJet.pt());
        jecAK4_->setJetE(uncorrJet.energy());
        jecAK4_->setRho(rhoVal_);
        jecAK4_->setNPV(vertices->size());
        jecAK4_->setJetA((*ak4jets)[ik].jetArea());
        double corr = jecAK4_->getCorrection();

        if (corr * uncorrJet.pt() > tmpjetptcut) {
            TLorentzVector* dummy = new TLorentzVector(0, 0, 0, 0);
            dummy->SetPtEtaPhiE(corr * uncorrJet.pt(), uncorrJet.eta(), uncorrJet.phi(), corr * uncorrJet.energy());
            jets.push_back(dummy);
            ++nujets;
        }
        if (ik < 6) {
            ak4jet_hf[ik] = (*ak4jets)[ik].hadronFlavour();
            ak4jet_pf[ik] = (*ak4jets)[ik].partonFlavour();
            ak4jet_pt[ik]   = corr * uncorrJet.pt();
            ak4jet_eta[ik] = (*ak4jets)[ik].eta();
            ak4jet_phi[ik] = (*ak4jets)[ik].phi();
            ak4jet_e[ik]    = corr * uncorrJet.energy();
            ak4jet_csv[ik]  = (*ak4jets)[ik].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
            ak4jet_icsv[ik] = (*ak4jets)[ik].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        }
    }
    // two leading jets with JER


    std::vector<TLorentzVector*> jets_new;
    for (size_t ik = 0; ik < ak4jets->size(); ik++) {
        reco::Candidate::LorentzVector uncorrJet = (*ak4jets)[ik].correctedP4(0);
        if ((*ak4jets)[ik].userFloat("SmearedPt") > tmpjetptcut) {
            TLorentzVector* dummy = new TLorentzVector(0, 0, 0, 0);
            dummy->SetPtEtaPhiE((*ak4jets)[ik].userFloat("SmearedPt"), uncorrJet.eta(), uncorrJet.phi(), (*ak4jets)[ik].userFloat("SmearedE"));
            jets_new.push_back(dummy);
        }
        if (ik < 6) {
            ak4jet_pt_new[ik] = (*ak4jets)[ik].userFloat("SmearedPt");
            ak4jet_e_new[ik]  = (*ak4jets)[ik].userFloat("SmearedE");
        }
    }
    //////----------------------------------------------
    //------------- jet pt energy JEC up uncertaity
    //////----------------------------------------------

    std::vector<TLorentzVector*> jets_JEC_up;
    for (size_t ik = 0; ik < ak4jets->size(); ik++) {
        reco::Candidate::LorentzVector uncorrJet = (*ak4jets)[ik].correctedP4(0);
        if ((*ak4jets)[ik].userFloat("SmearedPt_JEC_up") > tmpjetptcut) {
            TLorentzVector* dummy = new TLorentzVector(0, 0, 0, 0);
            dummy->SetPtEtaPhiE((*ak4jets)[ik].userFloat("SmearedPt_JEC_up"), uncorrJet.eta(), uncorrJet.phi(), (*ak4jets)[ik].userFloat("SmearedE_JEC_up"));
            jets_JEC_up.push_back(dummy);
        }
        if (ik < 6) {
            ak4jet_pt_JEC_up[ik] = (*ak4jets)[ik].userFloat("SmearedPt_JEC_up");
            ak4jet_e_JEC_up[ik]  = (*ak4jets)[ik].userFloat("SmearedE_JEC_up");
        }
    }
    //////----------------------------------------------
    ////------------- jet pt energy JEC down uncertaity
    ////////----------------------------------------------


    std::vector<TLorentzVector*> jets_JEC_down;
    for (size_t ik = 0; ik < ak4jets->size(); ik++) {
        reco::Candidate::LorentzVector uncorrJet = (*ak4jets)[ik].correctedP4(0);
        if ((*ak4jets)[ik].userFloat("SmearedPt_JEC_down") > tmpjetptcut) {
            TLorentzVector* dummy = new TLorentzVector(0, 0, 0, 0);
            dummy->SetPtEtaPhiE((*ak4jets)[ik].userFloat("SmearedPt_JEC_down"), uncorrJet.eta(), uncorrJet.phi(), (*ak4jets)[ik].userFloat("SmearedE_JEC_down"));
            jets_JEC_down.push_back(dummy);
        }
        if (ik < 6) {
            ak4jet_pt_JEC_down[ik] = (*ak4jets)[ik].userFloat("SmearedPt_JEC_down");
            ak4jet_e_JEC_down[ik]  = (*ak4jets)[ik].userFloat("SmearedE_JEC_down");
        }
    }
    //////----------------------------------------------
    ////------------- jet pt energy JER up uncertaity
    ////////----------------------------------------------


    std::vector<TLorentzVector*> jets_JER_up;
    for (size_t ik = 0; ik < ak4jets->size(); ik++) {
        reco::Candidate::LorentzVector uncorrJet = (*ak4jets)[ik].correctedP4(0);
        if ((*ak4jets)[ik].userFloat("SmearedPt_JER_up") > tmpjetptcut) {
            TLorentzVector* dummy = new TLorentzVector(0, 0, 0, 0);
            dummy->SetPtEtaPhiE((*ak4jets)[ik].userFloat("SmearedPt_JER_up"), uncorrJet.eta(), uncorrJet.phi(), (*ak4jets)[ik].userFloat("SmearedE_JER_up"));
            jets_JER_up.push_back(dummy);
        }
        if (ik < 6) {
            ak4jet_pt_JER_up[ik] = (*ak4jets)[ik].userFloat("SmearedPt_JER_up");
            ak4jet_e_JER_up[ik]  = (*ak4jets)[ik].userFloat("SmearedE_JER_up");
        }
    }
    //////----------------------------------------------
    ////------------- jet pt energy JER down uncertaity
    ////////----------------------------------------------

    std::vector<TLorentzVector*> jets_JER_down;
    for (size_t ik = 0; ik < ak4jets->size(); ik++) {
        reco::Candidate::LorentzVector uncorrJet = (*ak4jets)[ik].correctedP4(0);
        if ((*ak4jets)[ik].userFloat("SmearedPt_JER_down") > tmpjetptcut) {
            TLorentzVector* dummy = new TLorentzVector(0, 0, 0, 0);
            dummy->SetPtEtaPhiE((*ak4jets)[ik].userFloat("SmearedPt_JER_down"), uncorrJet.eta(), uncorrJet.phi(), (*ak4jets)[ik].userFloat("SmearedE_JER_down"));
            jets_JER_down.push_back(dummy);
        }
        if (ik < 6) {
            ak4jet_pt_JER_down[ik] = (*ak4jets)[ik].userFloat("SmearedPt_JER_down");
            ak4jet_e_JER_down[ik]  = (*ak4jets)[ik].userFloat("SmearedE_JER_down");
        }
    }

    sort(jets.begin(), jets.end(), mysortPt);
    sort(jets_new.begin(), jets_new.end(), mysortPt);
    sort(jets_JEC_up.begin(), jets_JEC_up.end(), mysortPt);
    sort(jets_JEC_down.begin(), jets_JEC_down.end(), mysortPt);
    sort(jets_JER_up.begin(), jets_JER_up.end(), mysortPt);
    sort(jets_JER_down.begin(), jets_JER_down.end(), mysortPt);

    //two leading jets
    for (size_t i = 0; i < jets.size(); i++) {
        if (iphoton > -1) {
            double drtmp1 = deltaR(jets.at(i)->Eta(), jets.at(i)->Phi(), photonsceta, photonscphi);
            if (drtmp1 > 0.5 && jetindexphoton12[0] == -1 && jetindexphoton12[1] == -1) {
                jetindexphoton12[0] = i;
                continue;  // the first num
            }
            if (drtmp1 > 0.5 && jetindexphoton12[0] != -1 && jetindexphoton12[1] == -1) {
                jetindexphoton12[1] = i;
                continue;  // the second num
            }
        }
    }

    for (size_t i = 0; i < jets.size(); i++) {
        if (iphoton_f > -1) {
            double drtmp1_f = deltaR(jets.at(i)->Eta(), jets.at(i)->Phi(), photonsceta_f, photonscphi_f);
            if (drtmp1_f > 0.5 && jetindexphoton12_f[0] == -1 && jetindexphoton12_f[1] == -1) {
                jetindexphoton12_f[0] = i;
                continue;  // the first num
            }
            if (drtmp1_f > 0.5 && jetindexphoton12_f[0] != -1 && jetindexphoton12_f[1] == -1) {
                jetindexphoton12_f[1] = i;
                continue;  // the second num
            }
        }
    }

    //two leading jets, new
    for (size_t i = 0; i < jets_new.size(); i++) {
        if (iphoton > -1) {
            double drtmp1 = deltaR(jets_new.at(i)->Eta(), jets_new.at(i)->Phi(), photoneta, photonphi);
            if (drtmp1 > 0.5 && jetindexphoton12_new[0] == -1 && jetindexphoton12_new[1] == -1) {
                jetindexphoton12_new[0] = i;
                continue;
            }
            if (drtmp1 > 0.5 && jetindexphoton12_new[0] != -1 && jetindexphoton12_new[1] == -1) {
                jetindexphoton12_new[1] = i;
                continue;
            }
        }
    }

    for (size_t i = 0; i < jets_new.size(); i++) {
        if (iphoton_f > -1) {
            double drtmp1_f = deltaR(jets_new.at(i)->Eta(), jets_new.at(i)->Phi(), photoneta_f, photonphi_f);
            if (drtmp1_f > 0.5 && jetindexphoton12_new_f[0] == -1 && jetindexphoton12_new_f[1] == -1) {
                jetindexphoton12_new_f[0] = i;
                continue;
            }
            if (drtmp1_f > 0.5 && jetindexphoton12_new_f[0] != -1 && jetindexphoton12_new_f[1] == -1) {
                jetindexphoton12_new_f[1] = i;
                continue;
            }
        }
    }


    //two leading jets, JEC up
    for (size_t i = 0; i < jets_JEC_up.size(); i++) {
        if (iphoton > -1) {
            double drtmp1 = deltaR(jets_JEC_up.at(i)->Eta(), jets_JEC_up.at(i)->Phi(), photoneta, photonphi);
            if (drtmp1 > 0.5 && jetindexphoton12_JEC_up[0] == -1 && jetindexphoton12_JEC_up[1] == -1) {
                jetindexphoton12_JEC_up[0] = i;
                continue;
            }
            if (drtmp1 > 0.5 && jetindexphoton12_JEC_up[0] != -1 && jetindexphoton12_JEC_up[1] == -1) {
                jetindexphoton12_JEC_up[1] = i;
                continue;
            }
        }
    }

    for (size_t i = 0; i < jets_JEC_up.size(); i++) {
        if (iphoton_f > -1) {
            double drtmp1_f = deltaR(jets_JEC_up.at(i)->Eta(), jets_JEC_up.at(i)->Phi(), photoneta_f, photonphi_f);
            if (drtmp1_f > 0.5 && jetindexphoton12_JEC_up_f[0] == -1 && jetindexphoton12_JEC_up_f[1] == -1) {
                jetindexphoton12_JEC_up_f[0] = i;
                continue;
            }
            if (drtmp1_f > 0.5 && jetindexphoton12_JEC_up_f[0] != -1 && jetindexphoton12_JEC_up_f[1] == -1) {
                jetindexphoton12_JEC_up_f[1] = i;
                continue;
            }
        }
    }
    //two leading jets, JEC down
    for (size_t i = 0; i < jets_JEC_down.size(); i++) {
        if (iphoton > -1) {
            double drtmp1 = deltaR(jets_JEC_down.at(i)->Eta(), jets_JEC_down.at(i)->Phi(), photoneta, photonphi);
            if (drtmp1 > 0.5 && jetindexphoton12_JEC_down[0] == -1 && jetindexphoton12_JEC_down[1] == -1) {
                jetindexphoton12_JEC_down[0] = i;
                continue;
            }
            if (drtmp1 > 0.5 && jetindexphoton12_JEC_down[0] != -1 && jetindexphoton12_JEC_down[1] == -1) {
                jetindexphoton12_JEC_down[1] = i;
                continue;
            }
        }
    }

    for (size_t i = 0; i < jets_JEC_down.size(); i++) {
        if (iphoton_f > -1) {
            double drtmp1_f = deltaR(jets_JEC_down.at(i)->Eta(), jets_JEC_down.at(i)->Phi(), photoneta_f, photonphi_f);
            if (drtmp1_f > 0.5 && jetindexphoton12_JEC_down_f[0] == -1 && jetindexphoton12_JEC_down_f[1] == -1) {
                jetindexphoton12_JEC_down_f[0] = i;
                continue;
            }
            if (drtmp1_f > 0.5 && jetindexphoton12_JEC_down_f[0] != -1 && jetindexphoton12_JEC_down_f[1] == -1) {
                jetindexphoton12_JEC_down_f[1] = i;
                continue;
            }
        }
    }

    //two leading jets, JER up
    for (size_t i = 0; i < jets_JER_up.size(); i++) {
        if (iphoton > -1) {
            double drtmp1 = deltaR(jets_JER_up.at(i)->Eta(), jets_JER_up.at(i)->Phi(), photoneta, photonphi);
            if (drtmp1 > 0.5 && jetindexphoton12_JER_up[0] == -1 && jetindexphoton12_JER_up[1] == -1) {
                jetindexphoton12_JER_up[0] = i;
                continue;
            }
            if (drtmp1 > 0.5 && jetindexphoton12_JER_up[0] != -1 && jetindexphoton12_JER_up[1] == -1) {
                jetindexphoton12_JER_up[1] = i;
                continue;
            }
        }
    }

    for (size_t i = 0; i < jets_JER_up.size(); i++) {
        if (iphoton_f > -1) {
            double drtmp1_f = deltaR(jets_JER_up.at(i)->Eta(), jets_JER_up.at(i)->Phi(), photoneta_f, photonphi_f);
            if (drtmp1_f > 0.5 && jetindexphoton12_JER_up_f[0] == -1 && jetindexphoton12_JER_up_f[1] == -1) {
                jetindexphoton12_JER_up_f[0] = i;
                continue;
            }
            if (drtmp1_f > 0.5 && jetindexphoton12_JER_up_f[0] != -1 && jetindexphoton12_JER_up_f[1] == -1) {
                jetindexphoton12_JER_up_f[1] = i;
                continue;
            }
        }
    }
    //two leading jets, JER down
    for (size_t i = 0; i < jets_JER_down.size(); i++) {
        if (iphoton > -1) {
            double drtmp1 = deltaR(jets_JER_down.at(i)->Eta(), jets_JER_down.at(i)->Phi(), photoneta, photonphi);
            if (drtmp1 > 0.5 && jetindexphoton12_JER_down[0] == -1 && jetindexphoton12_JER_down[1] == -1) {
                jetindexphoton12_JER_down[0] = i;
                continue;
            }
            if (drtmp1 > 0.5 && jetindexphoton12_JER_down[0] != -1 && jetindexphoton12_JER_down[1] == -1) {
                jetindexphoton12_JER_down[1] = i;
                continue;
            }
        }
    }
    for (size_t i = 0; i < jets_JER_down.size(); i++) {
        if (iphoton_f > -1) {
            double drtmp1_f = deltaR(jets_JER_down.at(i)->Eta(), jets_JER_down.at(i)->Phi(), photoneta_f, photonphi_f);
            if (drtmp1_f > 0.5 && jetindexphoton12_JER_down_f[0] == -1 && jetindexphoton12_JER_down_f[1] == -1) {
                jetindexphoton12_JER_down_f[0] = i;
                continue;
            }
            if (drtmp1_f > 0.5 && jetindexphoton12_JER_down_f[0] != -1 && jetindexphoton12_JER_down_f[1] == -1) {
                jetindexphoton12_JER_down_f[1] = i;
                continue;
            }
        }
    }


if(ak4jets->size()>=1){
        jet1hf_orig=(*ak4jets)[0].hadronFlavour();
        jet1pf_orig=(*ak4jets)[0].partonFlavour();
            jet1pt_orig=ak4jet_pt[0];
            jet1eta_orig=ak4jet_eta[0];
            jet1phi_orig=ak4jet_phi[0];
            jet1e_orig=ak4jet_e[0];
            jet1csv_orig =(*ak4jets)[0].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
            jet1icsv_orig =(*ak4jets)[0].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
            drj1a_orig=deltaR(jet1eta_orig,jet1phi_orig,photonsceta,photonscphi);
            drj1l_orig=deltaR(jet1eta_orig,jet1phi_orig,etalep1,philep1);

        if(ak4jets->size()>=2){

        jet2hf_orig=(*ak4jets)[1].hadronFlavour();
        jet2pf_orig=(*ak4jets)[1].partonFlavour();
            jet2pt_orig=ak4jet_pt[1];
            jet2eta_orig=ak4jet_eta[1];
            jet2phi_orig=ak4jet_phi[1];
            jet2e_orig=ak4jet_e[1];
            jet2csv_orig =(*ak4jets)[1].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
            jet2icsv_orig =(*ak4jets)[1].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
            drj2a_orig=deltaR(jet2eta_orig,jet2phi_orig,photonsceta,photonscphi);
            drj2l_orig=deltaR(jet2eta_orig,jet2phi_orig,etalep1,philep1);


        }


}


    // variable concerning jet, old
    if (jetindexphoton12[0] > -1 && jetindexphoton12[1] > -1) {

        jet1hf   = (*ak4jets)[jetindexphoton12[0]].hadronFlavour();
        jet1pf   = (*ak4jets)[jetindexphoton12[0]].partonFlavour();
        jet2hf   = (*ak4jets)[jetindexphoton12[1]].hadronFlavour();
        jet2pf   = (*ak4jets)[jetindexphoton12[1]].partonFlavour();
        jet1pt   = jets[jetindexphoton12[0]]->Pt();
        jet1eta  = jets[jetindexphoton12[0]]->Eta();
        jet1phi  = jets[jetindexphoton12[0]]->Phi();
        jet1e    = jets[jetindexphoton12[0]]->E();
        jet2pt   = jets[jetindexphoton12[1]]->Pt();
        jet2eta  = jets[jetindexphoton12[1]]->Eta();
        jet2phi  = jets[jetindexphoton12[1]]->Phi();
        jet2e    = jets[jetindexphoton12[1]]->E();
        jet1csv  = (*ak4jets)[jetindexphoton12[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet2csv  = (*ak4jets)[jetindexphoton12[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet1icsv = (*ak4jets)[jetindexphoton12[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet2icsv = (*ak4jets)[jetindexphoton12[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        drj1a    = deltaR(jet1eta, jet1phi, photonsceta, photonscphi);
        drj2a    = deltaR(jet2eta, jet2phi, photonsceta, photonscphi);
        drj1l    = deltaR(jet1eta, jet1phi, etalep1, philep1);
        drj2l    = deltaR(jet2eta, jet2phi, etalep1, philep1);
        TLorentzVector j1p4;
        j1p4.SetPtEtaPhiE(jet1pt, jet1eta, jet1phi, jet1e);
        TLorentzVector j2p4;
        j2p4.SetPtEtaPhiE(jet2pt, jet2eta, jet2phi, jet2e);
        TLorentzVector photonp42;
        photonp42.SetPtEtaPhiE(photonet, photoneta, photonphi, photone);
        TLorentzVector vp4;
        //            vp4.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
        vp4.SetPtEtaPhiE(WLeptonic.pt(), WLeptonic.eta(), WLeptonic.phi(), WLeptonic.energy());
        j1metPhi = fabs(jet1phi - MET_phi);
        if (j1metPhi > Pi) {
            j1metPhi = 2.0 * Pi - j1metPhi;
        }
        j2metPhi = fabs(jet2phi - MET_phi);
        if (j2metPhi > Pi) {
            j2metPhi = 2.0 * Pi - j2metPhi;
        }
        Mjj      = (j1p4 + j2p4).M();
        deltaeta = fabs(jet1eta - jet2eta);
        zepp     = fabs((vp4 + photonp42).Rapidity() - (j1p4.Rapidity() + j2p4.Rapidity()) / 2.0);
        Dphiwajj = fabs((vp4 + photonp42).Phi() - (j1p4 + j2p4).Phi());
        if (Dphiwajj > Pi) {
            Dphiwajj = 2.0 * Pi - Dphiwajj;
        }
    }
    // variable concerning jet,new
    if (jetindexphoton12_new[0] > -1 && jetindexphoton12_new[1] > -1) {
        jet1pt_new   = jets_new[jetindexphoton12_new[0]]->Pt();
        jet1eta_new  = jets_new[jetindexphoton12_new[0]]->Eta();
        jet1phi_new  = jets_new[jetindexphoton12_new[0]]->Phi();
        jet1e_new    = jets_new[jetindexphoton12_new[0]]->E();
        jet2pt_new   = jets_new[jetindexphoton12_new[1]]->Pt();
        jet2eta_new  = jets_new[jetindexphoton12_new[1]]->Eta();
        jet2phi_new  = jets_new[jetindexphoton12_new[1]]->Phi();
        jet2e_new    = jets_new[jetindexphoton12_new[1]]->E();
        jet1csv_new  = (*ak4jets)[jetindexphoton12_new[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet2csv_new  = (*ak4jets)[jetindexphoton12_new[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet1icsv_new = (*ak4jets)[jetindexphoton12_new[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet2icsv_new = (*ak4jets)[jetindexphoton12_new[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        drj1a_new    = deltaR(jet1eta_new, jet1phi_new, photoneta, photonphi);
        drj2a_new    = deltaR(jet2eta_new, jet2phi_new, photoneta, photonphi);
        drj1l_new    = deltaR(jet1eta_new, jet1phi_new, etalep1, philep1);
        drj2l_new    = deltaR(jet2eta_new, jet2phi_new, etalep1, philep1);
        TLorentzVector j1p4;
        j1p4.SetPtEtaPhiE(jet1pt_new, jet1eta_new, jet1phi_new, jet1e_new);
        TLorentzVector j2p4;
        j2p4.SetPtEtaPhiE(jet2pt_new, jet2eta_new, jet2phi_new, jet2e_new);
        TLorentzVector photonp42;
        photonp42.SetPtEtaPhiE(photonet, photoneta, photonphi, photone);
        TLorentzVector vp4;
        vp4.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
        j1metPhi_new = fabs(jet1phi_new - MET_phi_new);
        if (j1metPhi_new > Pi) {
            j1metPhi_new = 2.0 * Pi - j1metPhi_new;
        }
        j2metPhi_new = fabs(jet2phi_new - MET_phi_new);
        if (j2metPhi_new > Pi) {
            j2metPhi_new = 2.0 * Pi - j2metPhi_new;
        }
        Mjj_new      = (j1p4 + j2p4).M();
        zepp_new     = fabs((vp4 + photonp42).Rapidity() - (j1p4.Rapidity() + j2p4.Rapidity()) / 2.0);
        deltaeta_new = fabs(jet1eta_new - jet2eta_new);
        Dphiwajj_new=fabs((vp4+photonp42).Phi()-(j1p4+j2p4).Phi());
        if(Dphiwajj_new>Pi){Dphiwajj_new=2.0*Pi-Dphiwajj_new;}
    }
    // variable concerning jet, JEC up
    if (jetindexphoton12_JEC_up[0] > -1 && jetindexphoton12_JEC_up[1] > -1) {
        jet1pt_JEC_up   = jets_JEC_up[jetindexphoton12_JEC_up[0]]->Pt();
        jet1eta_JEC_up  = jets_JEC_up[jetindexphoton12_JEC_up[0]]->Eta();
        jet1phi_JEC_up  = jets_JEC_up[jetindexphoton12_JEC_up[0]]->Phi();
        jet1e_JEC_up    = jets_JEC_up[jetindexphoton12_JEC_up[0]]->E();
        jet2pt_JEC_up   = jets_JEC_up[jetindexphoton12_JEC_up[1]]->Pt();
        jet2eta_JEC_up  = jets_JEC_up[jetindexphoton12_JEC_up[1]]->Eta();
        jet2phi_JEC_up  = jets_JEC_up[jetindexphoton12_JEC_up[1]]->Phi();
        jet2e_JEC_up    = jets_JEC_up[jetindexphoton12_JEC_up[1]]->E();
        jet1csv_JEC_up  = (*ak4jets)[jetindexphoton12_JEC_up[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet2csv_JEC_up  = (*ak4jets)[jetindexphoton12_JEC_up[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet1icsv_JEC_up = (*ak4jets)[jetindexphoton12_JEC_up[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet2icsv_JEC_up = (*ak4jets)[jetindexphoton12_JEC_up[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        drj1a_JEC_up    = deltaR(jet1eta_JEC_up, jet1phi_JEC_up, photoneta, photonphi);
        drj2a_JEC_up    = deltaR(jet2eta_JEC_up, jet2phi_JEC_up, photoneta, photonphi);
        drj1l_JEC_up    = deltaR(jet1eta_JEC_up, jet1phi_JEC_up, etalep1, philep1);
        drj2l_JEC_up    = deltaR(jet2eta_JEC_up, jet2phi_JEC_up, etalep1, philep1);
        TLorentzVector j1p4;
        j1p4.SetPtEtaPhiE(jet1pt_JEC_up, jet1eta_JEC_up, jet1phi_JEC_up, jet1e_JEC_up);
        TLorentzVector j2p4;
        j2p4.SetPtEtaPhiE(jet2pt_JEC_up, jet2eta_JEC_up, jet2phi_JEC_up, jet2e_JEC_up);
        TLorentzVector photonp42;
        photonp42.SetPtEtaPhiE(photonet, photoneta, photonphi, photone);
        TLorentzVector vp4;
        vp4.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
        j1metPhi_JEC_up = fabs(jet1phi_JEC_up - MET_phi_JEC_up);
        if (j1metPhi_JEC_up > Pi) {
            j1metPhi_JEC_up = 2.0 * Pi - j1metPhi_JEC_up;
        }
        j2metPhi_JEC_up = fabs(jet2phi_JEC_up - MET_phi_JEC_up);
        if (j2metPhi_JEC_up > Pi) {
            j2metPhi_JEC_up = 2.0 * Pi - j2metPhi_JEC_up;
        }
        Mjj_JEC_up      = (j1p4 + j2p4).M();
        zepp_JEC_up     = fabs((vp4 + photonp42).Rapidity() - (j1p4.Rapidity() + j2p4.Rapidity()) / 2.0);
        deltaeta_JEC_up = fabs(jet1eta_JEC_up - jet2eta_JEC_up);
        Dphiwajj_JEC_up=fabs((vp4+photonp42).Phi()-(j1p4+j2p4).Phi());
        if(Dphiwajj_JEC_up>Pi){Dphiwajj_JEC_up=2.0*Pi-Dphiwajj_JEC_up;}
    }
    // variable concerning jet, JEC down
    if (jetindexphoton12_JEC_down[0] > -1 && jetindexphoton12_JEC_down[1] > -1) {
        jet1pt_JEC_down   = jets_JEC_down[jetindexphoton12_JEC_down[0]]->Pt();
        jet1eta_JEC_down  = jets_JEC_down[jetindexphoton12_JEC_down[0]]->Eta();
        jet1phi_JEC_down  = jets_JEC_down[jetindexphoton12_JEC_down[0]]->Phi();
        jet1e_JEC_down    = jets_JEC_down[jetindexphoton12_JEC_down[0]]->E();
        jet2pt_JEC_down   = jets_JEC_down[jetindexphoton12_JEC_down[1]]->Pt();
        jet2eta_JEC_down  = jets_JEC_down[jetindexphoton12_JEC_down[1]]->Eta();
        jet2phi_JEC_down  = jets_JEC_down[jetindexphoton12_JEC_down[1]]->Phi();
        jet2e_JEC_down    = jets_JEC_down[jetindexphoton12_JEC_down[1]]->E();
        jet1csv_JEC_down  = (*ak4jets)[jetindexphoton12_JEC_down[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet2csv_JEC_down  = (*ak4jets)[jetindexphoton12_JEC_down[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet1icsv_JEC_down = (*ak4jets)[jetindexphoton12_JEC_down[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet2icsv_JEC_down = (*ak4jets)[jetindexphoton12_JEC_down[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        drj1a_JEC_down    = deltaR(jet1eta_JEC_down, jet1phi_JEC_down, photoneta, photonphi);
        drj2a_JEC_down    = deltaR(jet2eta_JEC_down, jet2phi_JEC_down, photoneta, photonphi);
        drj1l_JEC_down    = deltaR(jet1eta_JEC_down, jet1phi_JEC_down, etalep1, philep1);
        drj2l_JEC_down    = deltaR(jet2eta_JEC_down, jet2phi_JEC_down, etalep1, philep1);
        //drj1l2_JEC_down   = deltaR(jet1eta_JEC_down, jet1phi_JEC_down, etalep2, philep2);
        //drj2l2_JEC_down   = deltaR(jet2eta_JEC_down, jet2phi_JEC_down, etalep2, philep2);
        TLorentzVector j1p4;
        j1p4.SetPtEtaPhiE(jet1pt_JEC_down, jet1eta_JEC_down, jet1phi_JEC_down, jet1e_JEC_down);
        TLorentzVector j2p4;
        j2p4.SetPtEtaPhiE(jet2pt_JEC_down, jet2eta_JEC_down, jet2phi_JEC_down, jet2e_JEC_down);
        TLorentzVector photonp42;
        photonp42.SetPtEtaPhiE(photonet, photoneta, photonphi, photone);
        TLorentzVector vp4;
        vp4.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
        j1metPhi_JEC_down = fabs(jet1phi_JEC_down - MET_phi_JEC_down);
        if (j1metPhi_JEC_down > Pi) {
            j1metPhi_JEC_down = 2.0 * Pi - j1metPhi_JEC_down;
        }
        j2metPhi_JEC_down = fabs(jet2phi_JEC_down - MET_phi_JEC_down);
        if (j2metPhi_JEC_down > Pi) {
            j2metPhi_JEC_down = 2.0 * Pi - j2metPhi_JEC_down;
        }
        Mjj_JEC_down      = (j1p4 + j2p4).M();
        zepp_JEC_down     = fabs((vp4 + photonp42).Rapidity() - (j1p4.Rapidity() + j2p4.Rapidity()) / 2.0);
        deltaeta_JEC_down = fabs(jet1eta_JEC_down - jet2eta_JEC_down);
        Dphiwajj_JEC_down=fabs((vp4+photonp42).Phi()-(j1p4+j2p4).Phi());
        if(Dphiwajj_JEC_down>Pi){Dphiwajj_JEC_down=2.0*Pi-Dphiwajj_JEC_down;}
    }
    // variable concerning jet, JER up
    if (jetindexphoton12_JER_up[0] > -1 && jetindexphoton12_JER_up[1] > -1) {
        jet1pt_JER_up   = jets_JER_up[jetindexphoton12_JER_up[0]]->Pt();
        jet1eta_JER_up  = jets_JER_up[jetindexphoton12_JER_up[0]]->Eta();
        jet1phi_JER_up  = jets_JER_up[jetindexphoton12_JER_up[0]]->Phi();
        jet1e_JER_up    = jets_JER_up[jetindexphoton12_JER_up[0]]->E();
        jet2pt_JER_up   = jets_JER_up[jetindexphoton12_JER_up[1]]->Pt();
        jet2eta_JER_up  = jets_JER_up[jetindexphoton12_JER_up[1]]->Eta();
        jet2phi_JER_up  = jets_JER_up[jetindexphoton12_JER_up[1]]->Phi();
        jet2e_JER_up    = jets_JER_up[jetindexphoton12_JER_up[1]]->E();
        jet1csv_JER_up  = (*ak4jets)[jetindexphoton12_JER_up[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet2csv_JER_up  = (*ak4jets)[jetindexphoton12_JER_up[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet1icsv_JER_up = (*ak4jets)[jetindexphoton12_JER_up[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet2icsv_JER_up = (*ak4jets)[jetindexphoton12_JER_up[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        drj1a_JER_up    = deltaR(jet1eta_JER_up, jet1phi_JER_up, photoneta, photonphi);
        drj2a_JER_up    = deltaR(jet2eta_JER_up, jet2phi_JER_up, photoneta, photonphi);
        drj1l_JER_up    = deltaR(jet1eta_JER_up, jet1phi_JER_up, etalep1, philep1);
        drj2l_JER_up    = deltaR(jet2eta_JER_up, jet2phi_JER_up, etalep1, philep1);
        //drj1l2_JER_up   = deltaR(jet1eta_JER_up, jet1phi_JER_up, etalep2, philep2);
        //drj2l2_JER_up   = deltaR(jet2eta_JER_up, jet2phi_JER_up, etalep2, philep2);
        TLorentzVector j1p4;
        j1p4.SetPtEtaPhiE(jet1pt_JER_up, jet1eta_JER_up, jet1phi_JER_up, jet1e_JER_up);
        TLorentzVector j2p4;
        j2p4.SetPtEtaPhiE(jet2pt_JER_up, jet2eta_JER_up, jet2phi_JER_up, jet2e_JER_up);
        TLorentzVector photonp42;
        photonp42.SetPtEtaPhiE(photonet, photoneta, photonphi, photone);
        TLorentzVector vp4;
        vp4.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
        j1metPhi_JER_up = fabs(jet1phi_JER_up - MET_phi_JER_up);
        if (j1metPhi_JER_up > Pi) {
            j1metPhi_JER_up = 2.0 * Pi - j1metPhi_JER_up;
        }
        j2metPhi_JER_up = fabs(jet2phi_JER_up - MET_phi_JER_up);
        if (j2metPhi_JER_up > Pi) {
            j2metPhi_JER_up = 2.0 * Pi - j2metPhi_JER_up;
        }
        Mjj_JER_up      = (j1p4 + j2p4).M();
        zepp_JER_up     = fabs((vp4 + photonp42).Rapidity() - (j1p4.Rapidity() + j2p4.Rapidity()) / 2.0);
        deltaeta_JER_up = fabs(jet1eta_JER_up - jet2eta_JER_up);
        Dphiwajj_JER_up=fabs((vp4+photonp42).Phi()-(j1p4+j2p4).Phi());
        if(Dphiwajj_JER_up>Pi){Dphiwajj_JER_up=2.0*Pi-Dphiwajj_JER_up;}
    }
    // variable concerning jet, JER down
    if (jetindexphoton12_JER_down[0] > -1 && jetindexphoton12_JER_down[1] > -1) {
        jet1pt_JER_down   = jets_JER_down[jetindexphoton12_JER_down[0]]->Pt();
        jet1eta_JER_down  = jets_JER_down[jetindexphoton12_JER_down[0]]->Eta();
        jet1phi_JER_down  = jets_JER_down[jetindexphoton12_JER_down[0]]->Phi();
        jet1e_JER_down    = jets_JER_down[jetindexphoton12_JER_down[0]]->E();
        jet2pt_JER_down   = jets_JER_down[jetindexphoton12_JER_down[1]]->Pt();
        jet2eta_JER_down  = jets_JER_down[jetindexphoton12_JER_down[1]]->Eta();
        jet2phi_JER_down  = jets_JER_down[jetindexphoton12_JER_down[1]]->Phi();
        jet2e_JER_down    = jets_JER_down[jetindexphoton12_JER_down[1]]->E();
        jet1csv_JER_down  = (*ak4jets)[jetindexphoton12_JER_down[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet2csv_JER_down  = (*ak4jets)[jetindexphoton12_JER_down[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet1icsv_JER_down = (*ak4jets)[jetindexphoton12_JER_down[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet2icsv_JER_down = (*ak4jets)[jetindexphoton12_JER_down[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        drj1a_JER_down    = deltaR(jet1eta_JER_down, jet1phi_JER_down, photoneta, photonphi);
        drj2a_JER_down    = deltaR(jet2eta_JER_down, jet2phi_JER_down, photoneta, photonphi);
        drj1l_JER_down    = deltaR(jet1eta_JER_down, jet1phi_JER_down, etalep1, philep1);
        drj2l_JER_down    = deltaR(jet2eta_JER_down, jet2phi_JER_down, etalep1, philep1);
        //drj1l2_JER_down   = deltaR(jet1eta_JER_down, jet1phi_JER_down, etalep2, philep2);
        //drj2l2_JER_down   = deltaR(jet2eta_JER_down, jet2phi_JER_down, etalep2, philep2);
        TLorentzVector j1p4;
        j1p4.SetPtEtaPhiE(jet1pt_JER_down, jet1eta_JER_down, jet1phi_JER_down, jet1e_JER_down);
        TLorentzVector j2p4;
        j2p4.SetPtEtaPhiE(jet2pt_JER_down, jet2eta_JER_down, jet2phi_JER_down, jet2e_JER_down);
        TLorentzVector photonp42;
        photonp42.SetPtEtaPhiE(photonet, photoneta, photonphi, photone);
        TLorentzVector vp4;
        vp4.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
        j1metPhi_JER_down = fabs(jet1phi_JER_down - MET_phi_JER_down);
        if (j1metPhi_JER_down > Pi) {
            j1metPhi_JER_down = 2.0 * Pi - j1metPhi_JER_down;
        }
        j2metPhi_JER_down = fabs(jet2phi_JER_down - MET_phi_JER_down);
        if (j2metPhi_JER_down > Pi) {
            j2metPhi_JER_down = 2.0 * Pi - j2metPhi_JER_down;
        }
        Mjj_JER_down      = (j1p4 + j2p4).M();
        zepp_JER_down     = fabs((vp4 + photonp42).Rapidity() - (j1p4.Rapidity() + j2p4.Rapidity()) / 2.0);
        deltaeta_JER_down = fabs(jet1eta_JER_down - jet2eta_JER_down);
        Dphiwajj_JER_down=fabs((vp4+photonp42).Phi()-(j1p4+j2p4).Phi());
        if(Dphiwajj_JER_down>Pi){Dphiwajj_JER_down=2.0*Pi-Dphiwajj_JER_down;}
        //	std::cout<<"Mjj_new "<<Mjj_new<<" Mjj_JEC_up "<<Mjj_JEC_up<<" Mjj_JEC_down "<<Mjj_JEC_down<<" Mjj_JER_up "<<Mjj_JER_up<<" Mjj_JER_down "<<Mjj_JER_down<<std::endl;
    }


std::cout<<"begin process old jets _f !!!"<<std::endl;

    if (jetindexphoton12_f[0] > -1 && jetindexphoton12_f[1] > -1) {

        jet1hf_f   = (*ak4jets)[jetindexphoton12_f[0]].hadronFlavour();
        jet1pf_f   = (*ak4jets)[jetindexphoton12_f[0]].partonFlavour();
        jet2hf_f   = (*ak4jets)[jetindexphoton12_f[1]].hadronFlavour();
        jet2pf_f   = (*ak4jets)[jetindexphoton12_f[1]].partonFlavour();
        jet1pt_f   = jets[jetindexphoton12_f[0]]->Pt();
        jet1eta_f  = jets[jetindexphoton12_f[0]]->Eta();
        jet1phi_f  = jets[jetindexphoton12_f[0]]->Phi();
        jet1e_f    = jets[jetindexphoton12_f[0]]->E();
        jet2pt_f   = jets[jetindexphoton12_f[1]]->Pt();
        jet2eta_f  = jets[jetindexphoton12_f[1]]->Eta();
        jet2phi_f  = jets[jetindexphoton12_f[1]]->Phi();
        jet2e_f    = jets[jetindexphoton12_f[1]]->E();
        jet1csv_f  = (*ak4jets)[jetindexphoton12_f[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet2csv_f  = (*ak4jets)[jetindexphoton12_f[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet1icsv_f = (*ak4jets)[jetindexphoton12_f[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet2icsv_f = (*ak4jets)[jetindexphoton12_f[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        drj1a_f    = deltaR(jet1eta_f, jet1phi_f, photonsceta_f, photonscphi_f);
        drj2a_f    = deltaR(jet2eta_f, jet2phi_f, photonsceta_f, photonscphi_f);
        drj1l_f    = deltaR(jet1eta_f, jet1phi_f, etalep1, philep1);
        drj2l_f    = deltaR(jet2eta_f, jet2phi_f, etalep1, philep1);
        TLorentzVector j1p4_f;
        j1p4_f.SetPtEtaPhiE(jet1pt_f, jet1eta_f, jet1phi_f, jet1e_f);
        TLorentzVector j2p4_f;
        j2p4_f.SetPtEtaPhiE(jet2pt_f, jet2eta_f, jet2phi_f, jet2e_f);
        TLorentzVector photonp42_f;
        photonp42_f.SetPtEtaPhiE(photonet_f, photoneta_f, photonphi_f, photone_f);
        TLorentzVector vp4_f;
        //            vp4.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
        vp4_f.SetPtEtaPhiE(WLeptonic.pt(), WLeptonic.eta(), WLeptonic.phi(), WLeptonic.energy());
        j1metPhi_f = fabs(jet1phi_f - MET_phi);
        if (j1metPhi_f > Pi) {
            j1metPhi_f = 2.0 * Pi - j1metPhi_f;
        }
        j2metPhi_f = fabs(jet2phi_f - MET_phi);
        if (j2metPhi_f > Pi) {
            j2metPhi_f = 2.0 * Pi - j2metPhi_f;
        }
        Mjj_f      = (j1p4_f + j2p4_f).M();
        deltaeta_f = fabs(jet1eta_f - jet2eta_f);
        zepp_f     = fabs((vp4_f + photonp42_f).Rapidity() - (j1p4_f.Rapidity() + j2p4_f.Rapidity()) / 2.0);
        Dphiwajj_f = fabs((vp4_f + photonp42_f).Phi() - (j1p4_f + j2p4_f).Phi());
        if (Dphiwajj_f > Pi) {
            Dphiwajj_f = 2.0 * Pi - Dphiwajj_f;
        }
    }
/*
std::cout<<"begin process old jets _f new !!!"<<std::endl;

    if (jetindexphoton12_new_f[0] > -1 && jetindexphoton12_new_f[1] > -1) {

std::cout<<"begin process old jets _f new array!!!"<<std::endl;

        jet1pt_new_f   = jets[jetindexphoton12_new_f[0]]->Pt();

std::cout<<"begin process old jets _f new array2!!!"<<std::endl;

        jet1eta_new_f  = jets[jetindexphoton12_new_f[0]]->Eta();

std::cout<<"begin process old jets _f new array3!!!"<<std::endl;

        jet1phi_new_f  = jets[jetindexphoton12_new_f[0]]->Phi();

std::cout<<"begin process old jets _f new array4!!!"<<std::endl;

        jet1e_new_f    = jets[jetindexphoton12_new_f[0]]->E();

std::cout<<"begin process old jets _f new array5!!!"<<std::endl;
std::cout<<"jetindexphoton12_new_f[1]  "<<jetindexphoton12_new_f[1]<<std::endl;
std::cout<<"  jets[jetindexphoton12_new_f[1]] pt  "<<jets[jetindexphoton12_new_f[1]]->Pt()<<std::endl;
        jet2pt_new_f   = jets[jetindexphoton12_new_f[1]]->Pt();

std::cout<<"begin process old jets _f new array6!!!"<<std::endl;

        jet2eta_new_f  = jets[jetindexphoton12_new_f[1]]->Eta();

std::cout<<"begin process old jets _f new array7!!!"<<std::endl;

        jet2phi_new_f  = jets[jetindexphoton12_new_f[1]]->Phi();

std::cout<<"begin process old jets _f new array8!!!"<<std::endl;

        jet2e_new_f    = jets[jetindexphoton12_new_f[1]]->E();

std::cout<<"begin process old jets _f new btag!!!"<<std::endl;

        jet1csv_new_f  = (*ak4jets)[jetindexphoton12_new_f[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet2csv_new_f  = (*ak4jets)[jetindexphoton12_new_f[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet1icsv_new_f = (*ak4jets)[jetindexphoton12_new_f[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet2icsv_new_f = (*ak4jets)[jetindexphoton12_new_f[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

std::cout<<"begin process old jets _f new dr"<<std::endl;

        drj1a_new_f    = deltaR(jet1eta_new_f, jet1phi_new_f, photoneta_f, photonphi_f);
        drj2a_new_f    = deltaR(jet2eta_new_f, jet2phi_new_f, photoneta_f, photonphi_f);
        drj1l_new_f    = deltaR(jet1eta_new_f, jet1phi_new_f, etalep1, philep1);
        drj2l_new_f    = deltaR(jet2eta_new_f, jet2phi_new_f, etalep1, philep1);
        //drj1l2_new_f   = deltaR(jet1eta_new_f, jet1phi_new_f, etalep2, philep2);
        //drj2l2_new_f   = deltaR(jet2eta_new_f, jet2phi_new_f, etalep2, philep2);

        TLorentzVector j1p4_f;
        j1p4_f.SetPtEtaPhiE(jet1pt_new_f, jet1eta_new_f, jet1phi_new_f, jet1e_new_f);
        TLorentzVector j2p4_f;
        j2p4_f.SetPtEtaPhiE(jet2pt_new_f, jet2eta_new_f, jet2phi_new_f, jet2e_new_f);
        TLorentzVector photonp42_f;
        photonp42_f.SetPtEtaPhiE(photonet_f, photoneta_f, photonphi_f, photone_f);
        TLorentzVector vp4_f;
        vp4_f.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());

std::cout<<"begin process old jets _f new phijmet!!!"<<std::endl;

        j1metPhi_new_f = fabs(jet1phi_new_f - MET_phi_new);
        if (j1metPhi_new_f > Pi) {
            j1metPhi_new_f = 2.0 * Pi - j1metPhi_new_f;
        }
        j2metPhi_new_f = fabs(jet2phi_new_f - MET_phi_new);
        if (j2metPhi_new_f > Pi) {
            j2metPhi_new_f = 2.0 * Pi - j2metPhi_new_f;
        }
        Mjj_new_f      = (j1p4_f + j2p4_f).M();

std::cout<<"begin process old jets _f new deltaeta!!!"<<std::endl;

        deltaeta_new_f = fabs(jet1eta_new_f - jet2eta_new_f);
        zepp_new_f     = fabs((vp4_f + photonp42_f).Rapidity() - (j1p4_f.Rapidity() + j2p4_f.Rapidity()) / 2.0);
    }

std::cout<<"begin process old jets _f jec up!!!"<<std::endl;

    if (jetindexphoton12_JEC_up_f[0] > -1 && jetindexphoton12_JEC_up_f[1] > -1) {
        jet1pt_JEC_up_f   = jets[jetindexphoton12_JEC_up_f[0]]->Pt();
        jet1eta_JEC_up_f  = jets[jetindexphoton12_JEC_up_f[0]]->Eta();
        jet1phi_JEC_up_f  = jets[jetindexphoton12_JEC_up_f[0]]->Phi();
        jet1e_JEC_up_f    = jets[jetindexphoton12_JEC_up_f[0]]->E();
        jet2pt_JEC_up_f   = jets[jetindexphoton12_JEC_up_f[1]]->Pt();
        jet2eta_JEC_up_f  = jets[jetindexphoton12_JEC_up_f[1]]->Eta();
        jet2phi_JEC_up_f  = jets[jetindexphoton12_JEC_up_f[1]]->Phi();
        jet2e_JEC_up_f    = jets[jetindexphoton12_JEC_up_f[1]]->E();
        jet1csv_JEC_up_f  = (*ak4jets)[jetindexphoton12_JEC_up_f[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet2csv_JEC_up_f  = (*ak4jets)[jetindexphoton12_JEC_up_f[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet1icsv_JEC_up_f = (*ak4jets)[jetindexphoton12_JEC_up_f[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet2icsv_JEC_up_f = (*ak4jets)[jetindexphoton12_JEC_up_f[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        drj1a_JEC_up_f    = deltaR(jet1eta_JEC_up_f, jet1phi_JEC_up_f, photoneta_f, photonphi_f);
        drj2a_JEC_up_f    = deltaR(jet2eta_JEC_up_f, jet2phi_JEC_up_f, photoneta_f, photonphi_f);
        drj1l_JEC_up_f    = deltaR(jet1eta_JEC_up_f, jet1phi_JEC_up_f, etalep1, philep1);
        drj2l_JEC_up_f    = deltaR(jet2eta_JEC_up_f, jet2phi_JEC_up_f, etalep1, philep1);
        //drj1l2_JEC_up_f   = deltaR(jet1eta_JEC_up_f, jet1phi_JEC_up_f, etalep2, philep2);
        //drj2l2_JEC_up_f   = deltaR(jet2eta_JEC_up_f, jet2phi_JEC_up_f, etalep2, philep2);
        TLorentzVector j1p4_f;
        j1p4_f.SetPtEtaPhiE(jet1pt_JEC_up_f, jet1eta_JEC_up_f, jet1phi_JEC_up_f, jet1e_JEC_up_f);
        TLorentzVector j2p4_f;
        j2p4_f.SetPtEtaPhiE(jet2pt_JEC_up_f, jet2eta_JEC_up_f, jet2phi_JEC_up_f, jet2e_JEC_up_f);
        TLorentzVector photonp42_f;
        photonp42_f.SetPtEtaPhiE(photonet_f, photoneta_f, photonphi_f, photone_f);
        TLorentzVector vp4_f;
        vp4_f.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
        j1metPhi_JEC_up_f = fabs(jet1phi_JEC_up_f - MET_phi_JEC_up);
        if (j1metPhi_JEC_up_f > Pi) {
            j1metPhi_JEC_up_f = 2.0 * Pi - j1metPhi_JEC_up_f;
        }
        j2metPhi_JEC_up_f = fabs(jet2phi_JEC_up_f - MET_phi_JEC_up);
        if (j2metPhi_JEC_up_f > Pi) {
            j2metPhi_JEC_up_f = 2.0 * Pi - j2metPhi_JEC_up_f;
        }
        Mjj_JEC_up_f      = (j1p4_f + j2p4_f).M();
        deltaeta_JEC_up_f = fabs(jet1eta_JEC_up_f - jet2eta_JEC_up_f);
        zepp_JEC_up_f     = fabs((vp4_f + photonp42_f).Rapidity() - (j1p4_f.Rapidity() + j2p4_f.Rapidity()) / 2.0);
    }

std::cout<<"begin process old jets _f jec down!!!"<<std::endl;

    if (jetindexphoton12_JEC_down_f[0] > -1 && jetindexphoton12_JEC_down_f[1] > -1) {
        jet1pt_JEC_down_f   = jets[jetindexphoton12_JEC_down_f[0]]->Pt();
        jet1eta_JEC_down_f  = jets[jetindexphoton12_JEC_down_f[0]]->Eta();
        jet1phi_JEC_down_f  = jets[jetindexphoton12_JEC_down_f[0]]->Phi();
        jet1e_JEC_down_f    = jets[jetindexphoton12_JEC_down_f[0]]->E();
        jet2pt_JEC_down_f   = jets[jetindexphoton12_JEC_down_f[1]]->Pt();
        jet2eta_JEC_down_f  = jets[jetindexphoton12_JEC_down_f[1]]->Eta();
        jet2phi_JEC_down_f  = jets[jetindexphoton12_JEC_down_f[1]]->Phi();
        jet2e_JEC_down_f    = jets[jetindexphoton12_JEC_down_f[1]]->E();
        jet1csv_JEC_down_f  = (*ak4jets)[jetindexphoton12_JEC_down_f[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet2csv_JEC_down_f  = (*ak4jets)[jetindexphoton12_JEC_down_f[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet1icsv_JEC_down_f = (*ak4jets)[jetindexphoton12_JEC_down_f[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet2icsv_JEC_down_f = (*ak4jets)[jetindexphoton12_JEC_down_f[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        drj1a_JEC_down_f    = deltaR(jet1eta_JEC_down_f, jet1phi_JEC_down_f, photoneta_f, photonphi_f);
        drj2a_JEC_down_f    = deltaR(jet2eta_JEC_down_f, jet2phi_JEC_down_f, photoneta_f, photonphi_f);
        drj1l_JEC_down_f    = deltaR(jet1eta_JEC_down_f, jet1phi_JEC_down_f, etalep1, philep1);
        drj2l_JEC_down_f    = deltaR(jet2eta_JEC_down_f, jet2phi_JEC_down_f, etalep1, philep1);
        //drj1l2_JEC_down_f   = deltaR(jet1eta_JEC_down_f, jet1phi_JEC_down_f, etalep2, philep2);
        //drj2l2_JEC_down_f   = deltaR(jet2eta_JEC_down_f, jet2phi_JEC_down_f, etalep2, philep2);
        TLorentzVector j1p4_f;
        j1p4_f.SetPtEtaPhiE(jet1pt_JEC_down_f, jet1eta_JEC_down_f, jet1phi_JEC_down_f, jet1e_JEC_down_f);
        TLorentzVector j2p4_f;
        j2p4_f.SetPtEtaPhiE(jet2pt_JEC_down_f, jet2eta_JEC_down_f, jet2phi_JEC_down_f, jet2e_JEC_down_f);
        TLorentzVector photonp42_f;
        photonp42_f.SetPtEtaPhiE(photonet_f, photoneta_f, photonphi_f, photone_f);
        TLorentzVector vp4_f;
        vp4_f.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
        j1metPhi_JEC_down_f = fabs(jet1phi_JEC_down_f - MET_phi_JEC_down);
        if (j1metPhi_JEC_down_f > Pi) {
            j1metPhi_JEC_down_f = 2.0 * Pi - j1metPhi_JEC_down_f;
        }
        j2metPhi_JEC_down_f = fabs(jet2phi_JEC_down_f - MET_phi_JEC_down);
        if (j2metPhi_JEC_down_f > Pi) {
            j2metPhi_JEC_down_f = 2.0 * Pi - j2metPhi_JEC_down_f;
        }
        Mjj_JEC_down_f      = (j1p4_f + j2p4_f).M();
        deltaeta_JEC_down_f = fabs(jet1eta_JEC_down_f - jet2eta_JEC_down_f);
        zepp_JEC_down_f     = fabs((vp4_f + photonp42_f).Rapidity() - (j1p4_f.Rapidity() + j2p4_f.Rapidity()) / 2.0);
    }

std::cout<<"begin process old jets _f jer up!!!"<<std::endl;

    if (jetindexphoton12_JER_up_f[0] > -1 && jetindexphoton12_JER_up_f[1] > -1) {
        jet1pt_JER_up_f   = jets[jetindexphoton12_JER_up_f[0]]->Pt();
        jet1eta_JER_up_f  = jets[jetindexphoton12_JER_up_f[0]]->Eta();
        jet1phi_JER_up_f  = jets[jetindexphoton12_JER_up_f[0]]->Phi();
        jet1e_JER_up_f    = jets[jetindexphoton12_JER_up_f[0]]->E();
        jet2pt_JER_up_f   = jets[jetindexphoton12_JER_up_f[1]]->Pt();
        jet2eta_JER_up_f  = jets[jetindexphoton12_JER_up_f[1]]->Eta();
        jet2phi_JER_up_f  = jets[jetindexphoton12_JER_up_f[1]]->Phi();
        jet2e_JER_up_f    = jets[jetindexphoton12_JER_up_f[1]]->E();
        jet1csv_JER_up_f  = (*ak4jets)[jetindexphoton12_JER_up_f[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet2csv_JER_up_f  = (*ak4jets)[jetindexphoton12_JER_up_f[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet1icsv_JER_up_f = (*ak4jets)[jetindexphoton12_JER_up_f[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet2icsv_JER_up_f = (*ak4jets)[jetindexphoton12_JER_up_f[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        drj1a_JER_up_f    = deltaR(jet1eta_JER_up_f, jet1phi_JER_up_f, photoneta_f, photonphi_f);
        drj2a_JER_up_f    = deltaR(jet2eta_JER_up_f, jet2phi_JER_up_f, photoneta_f, photonphi_f);
        drj1l_JER_up_f    = deltaR(jet1eta_JER_up_f, jet1phi_JER_up_f, etalep1, philep1);
        drj2l_JER_up_f    = deltaR(jet2eta_JER_up_f, jet2phi_JER_up_f, etalep1, philep1);
        //drj1l2_JER_up_f   = deltaR(jet1eta_JER_up_f, jet1phi_JER_up_f, etalep2, philep2);
        //drj2l2_JER_up_f   = deltaR(jet2eta_JER_up_f, jet2phi_JER_up_f, etalep2, philep2);
        TLorentzVector j1p4_f;
        j1p4_f.SetPtEtaPhiE(jet1pt_JER_up_f, jet1eta_JER_up_f, jet1phi_JER_up_f, jet1e_JER_up_f);
        TLorentzVector j2p4_f;
        j2p4_f.SetPtEtaPhiE(jet2pt_JER_up_f, jet2eta_JER_up_f, jet2phi_JER_up_f, jet2e_JER_up_f);
        TLorentzVector photonp42_f;
        photonp42_f.SetPtEtaPhiE(photonet_f, photoneta_f, photonphi_f, photone_f);
        TLorentzVector vp4_f;
        vp4_f.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
        j1metPhi_JER_up_f = fabs(jet1phi_JER_up_f - MET_phi_JER_up);
        if (j1metPhi_JER_up_f > Pi) {
            j1metPhi_JER_up_f = 2.0 * Pi - j1metPhi_JER_up_f;
        }
        j2metPhi_JER_up_f = fabs(jet2phi_JER_up_f - MET_phi_JER_up);
        if (j2metPhi_JER_up_f > Pi) {
            j2metPhi_JER_up_f = 2.0 * Pi - j2metPhi_JER_up_f;
        }
        Mjj_JER_up_f      = (j1p4_f + j2p4_f).M();
        deltaeta_JER_up_f = fabs(jet1eta_JER_up_f - jet2eta_JER_up_f);
        zepp_JER_up_f     = fabs((vp4_f + photonp42_f).Rapidity() - (j1p4_f.Rapidity() + j2p4_f.Rapidity()) / 2.0);
    }

std::cout<<"begin process old jets _f jer down!!!"<<std::endl;

    if (jetindexphoton12_JER_down_f[0] > -1 && jetindexphoton12_JER_down_f[1] > -1) {
        jet1pt_JER_down_f   = jets[jetindexphoton12_JER_down_f[0]]->Pt();
        jet1eta_JER_down_f  = jets[jetindexphoton12_JER_down_f[0]]->Eta();
        jet1phi_JER_down_f  = jets[jetindexphoton12_JER_down_f[0]]->Phi();
        jet1e_JER_down_f    = jets[jetindexphoton12_JER_down_f[0]]->E();
        jet2pt_JER_down_f   = jets[jetindexphoton12_JER_down_f[1]]->Pt();
        jet2eta_JER_down_f  = jets[jetindexphoton12_JER_down_f[1]]->Eta();
        jet2phi_JER_down_f  = jets[jetindexphoton12_JER_down_f[1]]->Phi();
        jet2e_JER_down_f    = jets[jetindexphoton12_JER_down_f[1]]->E();
        jet1csv_JER_down_f  = (*ak4jets)[jetindexphoton12_JER_down_f[0]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet2csv_JER_down_f  = (*ak4jets)[jetindexphoton12_JER_down_f[1]].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet1icsv_JER_down_f = (*ak4jets)[jetindexphoton12_JER_down_f[0]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet2icsv_JER_down_f = (*ak4jets)[jetindexphoton12_JER_down_f[1]].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        drj1a_JER_down_f    = deltaR(jet1eta_JER_down_f, jet1phi_JER_down_f, photoneta_f, photonphi_f);
        drj2a_JER_down_f    = deltaR(jet2eta_JER_down_f, jet2phi_JER_down_f, photoneta_f, photonphi_f);
        drj1l_JER_down_f    = deltaR(jet1eta_JER_down_f, jet1phi_JER_down_f, etalep1, philep1);
        drj2l_JER_down_f    = deltaR(jet2eta_JER_down_f, jet2phi_JER_down_f, etalep1, philep1);
        //drj1l2_JER_down_f   = deltaR(jet1eta_JER_down_f, jet1phi_JER_down_f, etalep2, philep2);
        //drj2l2_JER_down_f   = deltaR(jet2eta_JER_down_f, jet2phi_JER_down_f, etalep2, philep2);
        TLorentzVector j1p4_f;
        j1p4_f.SetPtEtaPhiE(jet1pt_JER_down_f, jet1eta_JER_down_f, jet1phi_JER_down_f, jet1e_JER_down_f);
        TLorentzVector j2p4_f;
        j2p4_f.SetPtEtaPhiE(jet2pt_JER_down_f, jet2eta_JER_down_f, jet2phi_JER_down_f, jet2e_JER_down_f);
        TLorentzVector photonp42_f;
        photonp42_f.SetPtEtaPhiE(photonet_f, photoneta_f, photonphi_f, photone_f);
        TLorentzVector vp4_f;
        vp4_f.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());
        j1metPhi_JER_down_f = fabs(jet1phi_JER_down_f - MET_phi_JER_down);
        if (j1metPhi_JER_down_f > Pi) {
            j1metPhi_JER_down_f = 2.0 * Pi - j1metPhi_JER_down_f;
        }
        j2metPhi_JER_down_f = fabs(jet2phi_JER_down_f - MET_phi_JER_down);
        if (j2metPhi_JER_down_f > Pi) {
            j2metPhi_JER_down_f = 2.0 * Pi - j2metPhi_JER_down_f;
        }
        Mjj_JER_down_f      = (j1p4_f + j2p4_f).M();
        deltaeta_JER_down_f = fabs(jet1eta_JER_down_f - jet2eta_JER_down_f);
        zepp_JER_down_f     = fabs((vp4_f + photonp42_f).Rapidity() - (j1p4_f.Rapidity() + j2p4_f.Rapidity()) / 2.0);
    }
*/

    outTree_->Fill();
    delete jecAK4_;
    jecAK4_ = 0;

}

//-------------------------------------------------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------------------------------------------------//
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


//define this as a plug-in
DEFINE_FWK_MODULE(PKUTreeMaker);
