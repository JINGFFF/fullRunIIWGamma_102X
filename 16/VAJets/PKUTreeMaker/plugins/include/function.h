#ifndef _function_
#define _function_
using namespace std;
       
//------------------------------------
double PKUTreeMaker::getJEC(reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_) {
    double jetCorrFactor = 1.;
    if (fabs(rawJetP4.eta()) < jetCorrEtaMax) {
        jecAK4_->setJetEta(rawJetP4.eta());
        jecAK4_->setJetPt(rawJetP4.pt());
        jecAK4_->setJetE(rawJetP4.energy());
        jecAK4_->setJetPhi(rawJetP4.phi());
        jecAK4_->setJetA(jet.jetArea());
        jecAK4_->setRho(*(rho_.product()));
        jecAK4_->setNPV(nVtx);
        jetCorrFactor = jecAK4_->getCorrection();
    }
    reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
    corrJetP4 *= jetCorrFactor;
    return jetCorrFactor;
}
//------------------------------------
double PKUTreeMaker::getJECOffset(reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_) {
    double jetCorrFactor = 1.;
    if (fabs(rawJetP4.eta()) < jetCorrEtaMax) {
        jecOffset_->setJetEta(rawJetP4.eta());
        jecOffset_->setJetPt(rawJetP4.pt());
        jecOffset_->setJetE(rawJetP4.energy());
        jecOffset_->setJetPhi(rawJetP4.phi());
        jecOffset_->setJetA(jet.jetArea());
        jecOffset_->setRho(*(rho_.product()));
        jecOffset_->setNPV(nVtx);
        jetCorrFactor = jecOffset_->getCorrection();
    }
    reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
    corrJetP4 *= jetCorrFactor;
    return jetCorrFactor;
}
//------------------------------------
void PKUTreeMaker::addTypeICorr(edm::Event const& event) {
    TypeICorrMap_.clear();
    edm::Handle<pat::JetCollection> jets_;
    event.getByToken(t1jetSrc_, jets_);
    event.getByToken(rhoToken_, rho_);
    edm::Handle<reco::VertexCollection> vertices_;
    event.getByToken(VertexToken_, vertices_);
    edm::Handle<edm::View<pat::Muon>> muons_;
    event.getByToken(t1muSrc_, muons_);
    bool                                      skipEM_                  = true;
    double                                    skipEMfractionThreshold_ = 0.9;
    bool                                      skipMuons_               = true;
    std::string                               skipMuonSelection_string = "isGlobalMuon | isStandAloneMuon";
    StringCutObjectSelector<reco::Candidate>* skipMuonSelection_       = new StringCutObjectSelector<reco::Candidate>(skipMuonSelection_string, true);
    double                                    jetCorrEtaMax_           = 9.9;
    double                                    type1JetPtThreshold_     = 15.0;
    double                                    corrEx                   = 0;
    double                                    corrEy                   = 0;
    double                                    corrSumEt                = 0;

    std::vector<JetCorrectorParameters> vPar;
    for (std::vector<std::string>::const_iterator payloadBegin = jecAK4chsLabels_.begin(), payloadEnd = jecAK4chsLabels_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    jecAK4_ = new FactorizedJetCorrector(vPar);
    vPar.clear();
    for (std::vector<std::string>::const_iterator payloadBegin = offsetCorrLabel_.begin(), payloadEnd = offsetCorrLabel_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    jecOffset_ = new FactorizedJetCorrector(vPar);
    vPar.clear();
    for (const pat::Jet& jet : *jets_) {
        double emEnergyFraction = jet.chargedEmEnergyFraction() + jet.neutralEmEnergyFraction();
        if (skipEM_ && emEnergyFraction > skipEMfractionThreshold_)
            continue;
        reco::Candidate::LorentzVector rawJetP4 = jet.correctedP4(0);
        double                         corr     = getJEC(rawJetP4, jet, jetCorrEtaMax_, jetCorrLabel_);


        if (skipMuons_) {
            const std::vector<reco::CandidatePtr>& cands = jet.daughterPtrVector();
            for (std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();
                 cand != cands.end(); ++cand) {
                const reco::PFCandidate* pfcand = dynamic_cast<const reco::PFCandidate*>(cand->get());
                const reco::Candidate*   mu     = (pfcand != 0 ? (pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
                if (mu != 0 && (*skipMuonSelection_)(*mu)) {
                    reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
                    rawJetP4 -= muonP4;
                }
            }
        }
        reco::Candidate::LorentzVector corrJetP4 = corr * rawJetP4;
        if (corrJetP4.pt() > type1JetPtThreshold_) {
            reco::Candidate::LorentzVector tmpP4              = jet.correctedP4(0);
            corr                                              = getJECOffset(tmpP4, jet, jetCorrEtaMax_, offsetCorrLabel_);
            reco::Candidate::LorentzVector rawJetP4offsetCorr = corr * rawJetP4;
            corrEx -= (corrJetP4.px() - rawJetP4offsetCorr.px());
            corrEy -= (corrJetP4.py() - rawJetP4offsetCorr.py());
            corrSumEt += (corrJetP4.Et() - rawJetP4offsetCorr.Et());
        }
    }
    TypeICorrMap_["corrEx"]    = corrEx;
    TypeICorrMap_["corrEy"]    = corrEy;
    TypeICorrMap_["corrSumEt"] = corrSumEt;
    delete jecAK4_;
    jecAK4_ = 0;
    delete jecOffset_;
    jecOffset_ = 0;
    delete skipMuonSelection_;
    skipMuonSelection_ = 0;
}
void PKUTreeMaker::addTypeICorr_user(edm::Event const& event) {
    TypeICorrMap_user_.clear();
    edm::Handle<pat::JetCollection> jets_;
    event.getByToken(t1jetSrc_user_, jets_);
    double corrEx_JEC         = 0;
    double corrEy_JEC         = 0;
    double corrSumEt_JEC      = 0;
    double corrEx_JEC_up      = 0;
    double corrEy_JEC_up      = 0;
    double corrSumEt_JEC_up   = 0;
    double corrEx_JEC_down    = 0;
    double corrEy_JEC_down    = 0;
    double corrSumEt_JEC_down = 0;

    double corrEx_JER         = 0;
    double corrEy_JER         = 0;
    double corrSumEt_JER      = 0;
    double corrEx_JER_up      = 0;
    double corrEy_JER_up      = 0;
    double corrSumEt_JER_up   = 0;
    double corrEx_JER_down    = 0;
    double corrEy_JER_down    = 0;
    double corrSumEt_JER_down = 0;


    for (const pat::Jet& jet : *jets_) {
        corrEx_JEC += jet.userFloat("corrEx_MET_JEC");
        corrEy_JEC += jet.userFloat("corrEy_MET_JEC");
        corrSumEt_JEC += jet.userFloat("corrSumEt_MET_JEC");
        corrEx_JEC_up += jet.userFloat("corrEx_MET_JEC_up");
        corrEy_JEC_up += jet.userFloat("corrEy_MET_JEC_up");
        corrSumEt_JEC_up += jet.userFloat("corrSumEt_MET_JEC_up");
        corrEx_JEC_down += jet.userFloat("corrEx_MET_JEC_down");
        corrEy_JEC_down += jet.userFloat("corrEy_MET_JEC_down");
        corrSumEt_JEC_down += jet.userFloat("corrSumEt_MET_JEC_down");
        corrEx_JER += jet.userFloat("corrEx_MET_JER");
        corrEy_JER += jet.userFloat("corrEy_MET_JER");
        corrSumEt_JER += jet.userFloat("corrSumEt_MET_JER");
        corrEx_JER_up += jet.userFloat("corrEx_MET_JER_up");
        corrEy_JER_up += jet.userFloat("corrEy_MET_JER_up");
        corrSumEt_JER_up += jet.userFloat("corrSumEt_MET_JER_up");
        corrEx_JER_down += jet.userFloat("corrEx_MET_JER_down");
        corrEy_JER_down += jet.userFloat("corrEy_MET_JER_down");
        corrSumEt_JER_down += jet.userFloat("corrSumEt_MET_JER_down");
    }
    TypeICorrMap_user_["corrEx_JEC"]         = corrEx_JEC;
    TypeICorrMap_user_["corrEy_JEC"]         = corrEy_JEC;
    TypeICorrMap_user_["corrSumEt_JEC"]      = corrSumEt_JEC;
    TypeICorrMap_user_["corrEx_JEC_up"]      = corrEx_JEC_up;
    TypeICorrMap_user_["corrEy_JEC_up"]      = corrEy_JEC_up;
    TypeICorrMap_user_["corrSumEt_JEC_up"]   = corrSumEt_JEC_up;
    TypeICorrMap_user_["corrEx_JEC_down"]    = corrEx_JEC_down;
    TypeICorrMap_user_["corrEy_JEC_down"]    = corrEy_JEC_down;
    TypeICorrMap_user_["corrSumEt_JEC_down"] = corrSumEt_JEC_down;

    TypeICorrMap_user_["corrEx_JER"]         = corrEx_JER;
    TypeICorrMap_user_["corrEy_JER"]         = corrEy_JER;
    TypeICorrMap_user_["corrSumEt_JER"]      = corrSumEt_JER;
    TypeICorrMap_user_["corrEx_JER_up"]      = corrEx_JER_up;
    TypeICorrMap_user_["corrEy_JER_up"]      = corrEy_JER_up;
    TypeICorrMap_user_["corrSumEt_JER_up"]   = corrSumEt_JER_up;
    TypeICorrMap_user_["corrEx_JER_down"]    = corrEx_JER_down;
    TypeICorrMap_user_["corrEy_JER_down"]    = corrEy_JER_down;
    TypeICorrMap_user_["corrSumEt_JER_down"] = corrSumEt_JER_down;
}
//------------------------------------
math::XYZTLorentzVector
PKUTreeMaker::getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType) {
    double leppt     = lep.Pt();
    double lepphi    = lep.Phi();
    double lepeta    = lep.Eta();
    double lepenergy = lep.Energy();

    double metpt  = MetPt;
    double metphi = MetPhi;

    double px  = metpt * cos(metphi);
    double py  = metpt * sin(metphi);
    double pz  = 0;
    double pxl = leppt * cos(lepphi);
    double pyl = leppt * sin(lepphi);
    double pzl = leppt * sinh(lepeta);
    double El  = lepenergy;
    double a   = pow(MW_, 2) + pow(px + pxl, 2) + pow(py + pyl, 2) - px * px - py * py - El * El + pzl * pzl;
    double b   = 2. * pzl;
    double A   = b * b - 4. * El * El;
    double B   = 2. * a * b;
    double C   = a * a - 4. * (px * px + py * py) * El * El;
    ///////////////////////////pz for fnal
    double M_mu = 0;
    //if(lepType==1)M_mu=0.105658367;//mu
    //if(lepType==0)M_mu=0.00051099891;//electron
    int type       = 2;  // use the small abs real root
    a              = MW_ * MW_ - M_mu * M_mu + 2.0 * pxl * px + 2.0 * pyl * py;
    A              = 4.0 * (El * El - pzl * pzl);
    B              = -4.0 * a * pzl;
    C              = 4.0 * El * El * (px * px + py * py) - a * a;
    double tmproot = B * B - 4.0 * A * C;
    if (tmproot < 0) {
        //std::cout << "Complex root detected, taking real part..." << std::endl;
        pz = -B / (2 * A);  // take real part of complex roots
    }
    else {
        double tmpsol1 = (-B + sqrt(tmproot)) / (2.0 * A);
        double tmpsol2 = (-B - sqrt(tmproot)) / (2.0 * A);
        //std::cout << " Neutrino Solutions: " << tmpsol1 << ", " << tmpsol2 << std::endl;
        if (type == 0) {
            // two real roots, pick the one closest to pz of muon
            if (TMath::Abs(tmpsol2 - pzl) < TMath::Abs(tmpsol1 - pzl)) {
                pz = tmpsol2;
            }
            else {
                pz = tmpsol1;
            }
            // if pz is > 300 pick the most central root
            if (abs(pz) > 300.) {
                if (TMath::Abs(tmpsol1) < TMath::Abs(tmpsol2)) {
                    pz = tmpsol1;
                }
                else {
                    pz = tmpsol2;
                }
            }
        }
        if (type == 1) {
            // two real roots, pick the one closest to pz of muon
            if (TMath::Abs(tmpsol2 - pzl) < TMath::Abs(tmpsol1 - pzl)) {
                pz = tmpsol2;
            }
            else {
                pz = tmpsol1;
            }
        }
        if (type == 2) {
            // pick the most central root.
            if (TMath::Abs(tmpsol1) < TMath::Abs(tmpsol2)) {
                pz = tmpsol1;
            }
            else {
                pz = tmpsol2;
            }
        }
        /*if (type == 3 ) {
         // pick the largest value of the cosine
         TVector3 p3w, p3mu;
         p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol1);
         p3mu.SetXYZ(pxl, pyl, pzl );
         double sinthcm1 = 2.*(p3mu.Perp(p3w))/MW_;
         p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol2);
         double sinthcm2 = 2.*(p3mu.Perp(p3w))/MW_;
         double costhcm1 = sqrt(1. - sinthcm1*sinthcm1);
         double costhcm2 = sqrt(1. - sinthcm2*sinthcm2);
         if ( costhcm1 > costhcm2 ) { pz = tmpsol1; otherSol_ = tmpsol2; }
         else { pz = tmpsol2;otherSol_ = tmpsol1; }
         }*/
        //end of type3
    }  //endl of if real root
    //dont correct pt neutrino
    math::XYZTLorentzVector outP4(px, py, pz, sqrt(px * px + py * py + pz * pz));
    return outP4;
}  //end neutrinoP4

//------------------------------------
bool PKUTreeMaker::hasMatchedPromptElectron(const reco::SuperClusterRef& sc, const edm::Handle<edm::View<pat::Electron>>& eleCol, const edm::Handle<reco::ConversionCollection>& convCol, const math::XYZPoint& beamspot, float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax) {
    //check if a given SuperCluster matches to at least one GsfElectron having zero expected inner hits
    //and not matching any conversion in the collection passing the quality cuts
    if (sc.isNull())
        return false;
    for (edm::View<pat::Electron>::const_iterator it = eleCol->begin(); it != eleCol->end(); ++it) {
        //match electron to supercluster
        if (it->superCluster() != sc)
            continue;
        //check expected inner hits
        if (it->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) > 0)
            continue;
        //check if electron is matching to a conversion
        if (ConversionTools::hasMatchedConversion(*it, convCol, beamspot))
            continue;
        return true;
    }
    return false;
}

/////////////////////////////////////////////////////////////
int PKUTreeMaker::matchToTrueLep(double lept_eta, double lept_phi,
                                 const edm::Handle<edm::View<reco::GenParticle>>& genParticles, double& dR, int& ispromptLep) {
    dR                                = 999;
    const reco::Candidate* closestLep = 0;

    int im = 0;
    for (size_t i = 0; i < genParticles->size(); i++) {
        const reco::Candidate* particle = &(*genParticles)[i];
        if (particle->pt() < 10)
            continue;

        if ((abs(particle->pdgId()) != 11 && abs(particle->pdgId()) != 13) || particle->status() != 1)

            continue;
        double dRtmp = deltaR(lept_eta, lept_phi, particle->eta(), particle->phi());
        if (dRtmp < dR) {
            dR         = dRtmp;
            im         = i;
            closestLep = particle;
        }
    }

    if (!(closestLep != 0 && dR < 0.3)) {
        return UNMATCHED;
    }
    ispromptLep = ((*genParticles)[im].isPromptFinalState() || (*genParticles)[im].isDirectPromptTauDecayProductFinalState());
    if (ispromptLep && dR < 0.3)
        ispromptLep = 1;
    else
        ispromptLep = 0;
    return 1;
}
/////////////////////////////////////////////////////////////
//------------------------------------
int PKUTreeMaker::matchToTruth(const reco::Photon&                              pho,
                               const edm::Handle<edm::View<reco::GenParticle>>& genParticles, bool& ISRPho, double& dR, int& isprompt) {
    //
    // Explicit loop and geometric matching method
    //
    // Find the closest status 1 gen photon to the reco photon
    dR                                   = 999;
    const reco::Candidate* closestPhoton = 0;
    //std::cout<<"genParticles->size() = "<<genParticles->size()<<std::endl;
    int im = 0;
    for (size_t i = 0; i < genParticles->size(); i++) {
        const reco::Candidate* particle = &(*genParticles)[i];
        if (particle->pt() < 10)
            continue;                                                                                   /////////////////////////////////////////////////////////////
        if ((abs(particle->pdgId()) != 11 && abs(particle->pdgId()) != 22) || particle->status() != 1)  /////////////////////////////////////////////////////////////
                                                                                                        //        if( abs(particle->pdgId()) != 22 || particle->status() != 1 )
            continue;
        double dRtmp = deltaR(pho.eta(), pho.phi(), particle->eta(), particle->phi());
        if (dRtmp < dR) {
            dR            = dRtmp;
            im            = i;
            closestPhoton = particle;
        }
    }
    // See if the closest photon (if it exists) is close enough.
    // If not, no match found.
    if (!(closestPhoton != 0 && dR < 0.3)) {
        return UNMATCHED;
        // ISRPho = false;
    }
    //     isprompt=(*genParticles)[im].isPromptFinalState();
    isprompt                        = ((*genParticles)[im].isPromptFinalState() || (*genParticles)[im].isDirectPromptTauDecayProductFinalState());
    const reco::Candidate* particle = &(*genParticles)[im];
    if (abs(particle->pdgId()) == 11 && isprompt && dR < 0.3)
        isprompt = 3;
    else if (abs(particle->pdgId()) == 22 && isprompt && dR < 0.3)
        isprompt = 2;
    else
        isprompt = 0;

    // Find ID of the parent of the found generator level photon match
    int ancestorPID    = -999;
    int ancestorStatus = -999;
    findFirstNonPhotonMother(closestPhoton, ancestorPID, ancestorStatus);
    // Allowed parens: quarks pdgId 1-5, or a gluon 21
    std::vector<int> allowedParents{-1, 1, -2, 2, -3, 3, -4, 4, -5, 5, 21, -11, 11, -13, 13, -15, 15, 23, -24, 24};
    if (!(std::find(allowedParents.begin(),
                    allowedParents.end(), ancestorPID)
          != allowedParents.end())) {
        // So it is not from g, u, d, s, c, b. Check if it is from pi0 or not.
        if (abs(ancestorPID) == 111)
            return MATCHED_FROM_PI0;
        // ISRPho =true;
        else
            //      std::cout<<"Mother = "<<abs(ancestorPID)<<" "<<closestPhoton->mother(0)->pdgId()<<" "<<closestPhoton->mother(0)->status()<<std::endl;
            //      std::cout<<"Run="<<run<<" Event="<<nevent<<" lumi="<<ls<<std::endl;
            return MATCHED_FROM_OTHER_SOURCES;
        //  ISRPho =true;
    }
    return MATCHED_FROM_GUDSCB;  //Z, W ?
                                 //   ISRPho =true;
}
//------------------------------------

void PKUTreeMaker::findFirstNonPhotonMother(const reco::Candidate* particle,
                                            int& ancestorPID, int& ancestorStatus) {
    if (particle == 0) {
        printf("SimplePhotonNtupler: ERROR! null candidate pointer, this should never happen\n");
        return;
    }
    // Is this the first non-photon parent? If yes, return, otherwise
    // go deeper into recursion
    if (abs(particle->pdgId()) == 22) {
        findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
    }
    else {
        ancestorPID    = particle->pdgId();
        ancestorStatus = particle->status();
    }
    return;
}

#endif
