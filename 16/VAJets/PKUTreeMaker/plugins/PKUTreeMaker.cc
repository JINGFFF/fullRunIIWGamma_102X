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
#include "include/jet_info.h"


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
	hlt_info(iEvent);
	
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

	// photon info
	photon_info(iEvent);

	// jets info
	jet_info(iEvent);

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
