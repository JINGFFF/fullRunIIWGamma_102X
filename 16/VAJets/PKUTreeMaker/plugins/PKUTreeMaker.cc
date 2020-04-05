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
#include "include/hlt_info.h"
#include "include/weight_info.h"
#include "include/gen_info.h"
#include "include/filter_info.h"
#include "include/met_info.h"
#include "include/leptonicV_info.h"
#include "include/photon_info.h"



using namespace std;

       
//------------------------------------
//------------------------------------
void PKUTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;

	//Initalize variables with dummy values
    setDummyValues(); 

    nevent = iEvent.eventAuxiliary().event();
    run    = iEvent.eventAuxiliary().run();
    ls     = iEvent.eventAuxiliary().luminosityBlock();

    //events weight
    if (RunOnMC_) {

		weight_info(iEvent);

    }

	// store HLT Info
	hlt_info(edm::Event const & iEvent);
	
	//Gen info for photon, lepton, and jets
    if (RunOnMC_) {
		gen_photon_lepton_info(iEvent);
		gen_jet_info(iEvent);
	}

    //filter
	filter_info(iEvent);	

	//Vertice
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(VertexToken_, vertices);
    if (vertices->empty()) {
        outTree_->Fill();
        return;
    }  
	
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
    } 

    // MET info
	met_info(iEvent);
    
	// lepton, V info
	leptonicV_info(iEvent);

    if (is_leptonicVs_Empty) {
        outTree_->Fill();
        return;
    }  

    if (RunOnMC_ && ptlep1 > 10) {
        lepton_istrue = matchToTrueLep(etalep1, philep1, genParticles, dR1_, ispromptLep_);
    }

	// photon info
	photon_info(iEvent);


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

//define this as a plug-in
DEFINE_FWK_MODULE(PKUTreeMaker);
