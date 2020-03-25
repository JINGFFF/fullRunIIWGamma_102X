
/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 *
 *
 * Authors:
 *
 *   Chayanit Asawatangtrakuldee chayanit@cern.ch 
 *
 * Description:
 *   - Selects "loose" and "tight" electrons needed for V-boson analysis.
 *   - Saves collection of the reference vectors of electrons passing the 
 *     required electron ID.
 * History:
 *   
 *
 *****************************************************************************/
////////////////////////////////////////////////////////////////////////////////
// Includes
////////////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <memory>
#include <vector>
#include <sstream>
#include <cmath>
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////
class ElectronIdSelector : public edm::EDProducer
{
public:
  	// construction/destruction
  	ElectronIdSelector(const edm::ParameterSet& iConfig);
  	virtual ~ElectronIdSelector();
  
  	// member functions
  	void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
  	void endJob();

private:  
  	// member data
  	//edm::InputTag  src_;
  	std::string    moduleLabel_;
  	std::string    idLabel_;  
  	bool           useDetectorIsolation_;
  	bool           applyTightID_;
  	bool           applyMediumID_;
  	bool           applyLooseID_;
  	bool           applyVetoID_;
    bool           applyFakeID_;
  	unsigned int nTot_;
  	unsigned int nPassed_;
  	edm::EDGetTokenT<pat::ElectronCollection> ElectronToken_;
  	edm::EDGetTokenT<reco::VertexCollection> VertexToken_;
  	edm::EDGetTokenT<double> RhoToken_;
  	EffectiveAreas effectiveAreas_;
};



////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
ElectronIdSelector::ElectronIdSelector(const edm::ParameterSet& iConfig)
  	: moduleLabel_(iConfig.getParameter<std::string>   ("@module_label"))
  	, idLabel_(iConfig.existsAs<std::string>("idLabel") ? iConfig.getParameter<std::string>("idLabel") : "loose")
  	, useDetectorIsolation_(iConfig.existsAs<bool>("useDetectorIsolation") ? iConfig.getParameter<bool>("useDetectorIsolation") : false)
  	, nTot_(0)
  	, nPassed_(0)
  	, ElectronToken_ (consumes<pat::ElectronCollection> (iConfig.getParameter<edm::InputTag>( "src" ) ) ) 
  	, VertexToken_ (consumes<reco::VertexCollection> (iConfig.getParameter<edm::InputTag>( "vertex" ) ) )
  	, RhoToken_ (consumes<double> (iConfig.getParameter<edm::InputTag>( "rho") ) )
  	, effectiveAreas_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath() )
{
  	produces<std::vector<pat::Electron> >();

  	/// ------- Decode the ID criteria --------
  	applyTightID_ 	= false;
  	applyMediumID_ 	= false;
  	applyLooseID_ 	= false;
  	applyVetoID_ 	= false;
    applyFakeID_    = false;

  	if( (idLabel_.compare("tight")==0) || 
      	(idLabel_.compare("Tight")==0) || 
      	(idLabel_.compare("TIGHT")==0) ||
      	(idLabel_.compare("WP70")==0) ||
      	(idLabel_.compare("wp70")==0) )  
 	applyTightID_ = true;

  	else if( (idLabel_.compare("medium")==0) ||
      	(idLabel_.compare("Medium")==0) ||
      	(idLabel_.compare("MEDIUM")==0) ||
      	(idLabel_.compare("WP80")==0) ||
      	(idLabel_.compare("wp80")==0) )  
	applyMediumID_ = true;

  	else if( (idLabel_.compare("loose")==0) || 
      	(idLabel_.compare("Loose")==0) || 
      	(idLabel_.compare("LOOSE")==0) ||
      	(idLabel_.compare("WP90")==0) ||
      	(idLabel_.compare("wp90")==0) )  
	applyLooseID_ = true;

  	else if( (idLabel_.compare("veto")==0) || 
      	(idLabel_.compare("Veto")==0) || 
      	(idLabel_.compare("VETO")==0) ||
      	(idLabel_.compare("VETOid")==0) ||
      	(idLabel_.compare("VetoId")==0) )  
	applyVetoID_ = true;

    else if( (idLabel_.compare("fake")==0) ||
        (idLabel_.compare("Fake")==0) ||
        (idLabel_.compare("FAKE")==0) )
    applyFakeID_ = true;
}

 
//______________________________________________________________________________
ElectronIdSelector::~ElectronIdSelector(){}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////
 
//______________________________________________________________________________
void ElectronIdSelector::produce(edm::Event& iEvent,const edm::EventSetup& iSetup)
{

   	edm::Handle<reco::VertexCollection> vtxs;
   	iEvent.getByToken(VertexToken_, vtxs);

  	std::unique_ptr<std::vector<pat::Electron> > passingElectrons(new std::vector<pat::Electron >);
   	edm::Handle<pat::ElectronCollection > electrons;
   	iEvent.getByToken(ElectronToken_, electrons);  

  	bool* isPassing = new bool[electrons->size()];
  	double rhoVal_;
  	rhoVal_=-99.;
  	edm::Handle<double> rho;
  	iEvent.getByToken(RhoToken_,rho);
  	rhoVal_ = *rho;

  	for(unsigned int iElec=0; iElec<electrons->size(); iElec++) { 

    	isPassing[iElec]=false;
    	const pat::Electron& ele = electrons->at(iElec);

    	// -------- Make sure that the electron is within acceptance ------
    	float eta = ele.superCluster()->eta();
    	bool isEB = ele.isEB() && fabs(eta) < 1.479;
    	bool isEE = ele.isEE() && fabs(eta) > 1.479 && fabs(eta) < 2.5;

    	//bool inAcceptance = (isEB || isEE) && (ele.ecalDrivenSeed()==1);
    	float pt  = ele.pt();
    	float energy = ele.superCluster()->energy();

    	// -------- Compute Detector isolation ------
    	float eA = effectiveAreas_.getEffectiveArea(fabs(eta));
    	float pf_isolation = ( ( ele.pfIsolationVariables().sumChargedHadronPt
                                   + std::max( 0., ele.pfIsolationVariables().sumNeutralHadronEt + ele.pfIsolationVariables().sumPhotonEt - eA*rhoVal_) )
                                 / pt );
    	float isolation = 100.;
    	isolation = pf_isolation;

    	// -------- Compute ID ------
    	double sigmaIEtaIEta   = ele.full5x5_sigmaIetaIeta();
    	double dPhiIn    = fabs(ele.deltaPhiSuperClusterTrackAtVtx());
    	double dEtaIn    = fabs(ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull() ? ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta() : std::numeric_limits<float>::max());
    	double hoe     = ele.hadronicOverEm();
    	double ooemoop = fabs((1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy()));

    	// impact parameter variables
    	float d0vtx         = 0.0;
    	float dzvtx         = 0.0;
    	if (vtxs->size() > 0) {
        	reco::VertexRef vtx(vtxs, 0);    
        	d0vtx = ele.gsfTrack()->dxy(vtx->position());
        	dzvtx = ele.gsfTrack()->dz(vtx->position());
    	} 
		else {
        	d0vtx = ele.gsfTrack()->dxy();
        	dzvtx = ele.gsfTrack()->dz();
    	}

      	bool vtxFitConversion = !(ele.passConversionVeto()==1);
    	float mHits=ele.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);  

    	bool isTight  	= false;  /////// <--- equivalent to WP70
    	bool isMedium 	= false;  /////// <--- equivalent to WP80
    	bool isLoose  	= false;  /////// <--- equivalent to WP90
    	bool isVeto    	= false;  /////// <--- the loosest cut for veto
		bool isFake		= false;

    	// ---------- cut-based ID -----------------
		//2016 data  94X MC
    	isTight = (pt>20.)  &&
      			(!vtxFitConversion) &&
      			((isEB && mHits<=1 && isolation<0.0287+0.506/pt && sigmaIEtaIEta<0.0104 && dPhiIn<0.022 && dEtaIn<0.00255 && hoe<0.026+1.15/energy+0.0324*rhoVal_/energy && ooemoop<0.159 && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10 )  ||
      			(isEE && mHits<=1 && isolation<0.0445+0.963/pt && sigmaIEtaIEta<0.0353 && dPhiIn<0.0236 && dEtaIn<0.00501 && hoe<0.0188+2.06/energy + 0.183*rhoVal_/energy && ooemoop<0.0197 && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20));

    	isMedium = (pt>20.)  &&
        		(!vtxFitConversion) &&
        		((isEB && mHits<=1 && isolation<0.0478+0.506/pt && sigmaIEtaIEta<0.0106 && dPhiIn<0.0547 && dEtaIn<0.0032 && hoe<0.046+1.16/energy+0.0324*rhoVal_/energy && ooemoop<0.184 && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10) ||
         		(isEE && mHits<=1 && isolation<0.0658+0.963/pt && sigmaIEtaIEta<0.0387 && dPhiIn<0.0394 && dEtaIn<0.00632 && hoe<0.0275+2.52/energy+0.183*rhoVal_/energy && ooemoop<0.0721 && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20));

    	isLoose = (pt>20.)  &&
         		(!vtxFitConversion) &&
        		((isEB && mHits<=1 && isolation<0.112+0.506/pt && sigmaIEtaIEta<0.0112 && dPhiIn<0.0884  && dEtaIn<0.00377  && hoe<0.05+1.16/energy+0.0324*rhoVal_/energy && ooemoop<0.193 && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10) ||
         		(isEE && mHits<=1 && isolation<0.108+0.963/pt && sigmaIEtaIEta<0.0425  && dPhiIn<0.169 && dEtaIn<0.00674 && hoe<0.0441+2.54/energy+0.183*rhoVal_/energy && ooemoop<0.111 && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20));

    	isVeto = (pt>20.) &&
         		(!vtxFitConversion) &&
        		((isEB && mHits<=2 && isolation<0.198+0.506/pt && sigmaIEtaIEta<0.0126 && dPhiIn<0.148 && dEtaIn<0.00463 && hoe<0.05+1.16/energy+0.0324*rhoVal_/energy && ooemoop<0.209 && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10 ) ||
         		(isEE && mHits<=3 && isolation<0.203+0.963/pt && sigmaIEtaIEta<0.0457 && dPhiIn<0.19 && dEtaIn<0.00814 && hoe<0.05+2.54/energy+0.183*rhoVal_/energy && ooemoop<0.132 && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20));

		isFake = (!isTight) && isVeto;

/*
	//2016 data    80X
    isTight = (pt>20.)  &&
      (!vtxFitConversion) &&
      ((isEB && mHits<=1 && isolation<0.0588 && sigmaIEtaIEta<0.00998 && dPhiIn<0.0816 && dEtaIn<0.00308 && hoe<0.0414 && ooemoop<0.0129 && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10 )  || (isEE && mHits<=1 && isolation<0.0571 && sigmaIEtaIEta<0.0292 && dPhiIn<0.0394 && dEtaIn<0.00605 && hoe<0.0641 && ooemoop<0.0129 && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20));

//fail tight ID, pass veto ID
//    isTight = (pt>20.)  &&
//      (!vtxFitConversion) &&
//        ((isEB && mHits<=2 && isolation<0.175 && sigmaIEtaIEta<0.0115 && dPhiIn<0.228 && dEtaIn<0.00749 && hoe<0.356 && ooemoop<0.299  && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10 && !(mHits<=1 && isolation<0.0588 && sigmaIEtaIEta<0.00998 && dPhiIn<0.0816 && dEtaIn<0.00308 && hoe<0.0414 && ooemoop<0.0129)) ||
//         (isEE && mHits<=3 && isolation<0.159 && sigmaIEtaIEta<0.037 && dPhiIn<0.213 && dEtaIn<0.00895 && hoe<0.211 && ooemoop<0.15  && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20 && !(mHits<=1 && isolation<0.0571 && sigmaIEtaIEta<0.0292 && dPhiIn<0.0394 && dEtaIn<0.00605 && hoe<0.0641 && ooemoop<0.0129)));

    isMedium = (pt>10.)  &&
        (!vtxFitConversion) &&
        ((isEB && mHits<=1 && isolation<0.0695 && sigmaIEtaIEta<0.00998 && dPhiIn<0.103 && dEtaIn<0.00311 && hoe<0.253 && ooemoop<0.134 && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10) ||
         (isEE && mHits<=1 && isolation<0.0821 && sigmaIEtaIEta<0.0298 && dPhiIn<0.045 && dEtaIn<0.00609 && hoe<0.0878 && ooemoop<0.13 && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20));

    isLoose = (pt>10.)  &&
         (!vtxFitConversion) &&
        ((isEB && mHits<=1 && isolation<0.0994 && sigmaIEtaIEta<0.011 && dPhiIn<0.222  && dEtaIn<0.00477  && hoe<0.298 && ooemoop<0.241 && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10) ||
         (isEE && mHits<=1 && isolation<0.107 && sigmaIEtaIEta<0.0314  && dPhiIn<0.213 && dEtaIn<0.00868 && hoe<0.101 && ooemoop<0.14 && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20));

    isVeto = (pt>10.) && 
         (!vtxFitConversion) &&
        ((isEB && mHits<=2 && isolation<0.175 && sigmaIEtaIEta<0.0115 && dPhiIn<0.228 && dEtaIn<0.00749 && hoe<0.356 && ooemoop<0.299  && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10 ) ||
         (isEE && mHits<=3 && isolation<0.159 && sigmaIEtaIEta<0.037 && dPhiIn<0.213 && dEtaIn<0.00895 && hoe<0.211 && ooemoop<0.15  && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20));
*/

    	/// ------- Finally apply selection --------
    	if(applyTightID_ && isTight)   	isPassing[iElec]= true;
    	if(applyMediumID_ && isMedium) 	isPassing[iElec]= true;
    	if(applyLooseID_ && isLoose)  	isPassing[iElec]= true;
    	if(applyVetoID_ && isVeto) 		isPassing[iElec]= true;
        if(applyFakeID_ && isFake)      isPassing[iElec]= true;
    
 	}
  
	for (unsigned int iElectron = 0; iElectron < electrons -> size(); iElectron ++){     
		if(isPassing[iElectron]) passingElectrons->push_back( electrons -> at(iElectron) );       
  	}
 
  	nTot_  +=electrons->size();
  	nPassed_+=passingElectrons->size();

  	delete [] isPassing;  
  	iEvent.put(std::move(passingElectrons));
}

//______________________________________________________________________________
void ElectronIdSelector::endJob()
{
  	std::stringstream ss;
  	ss<<"nTot="<<nTot_<<" nPassed="<<nPassed_
    	<<" effPassed="<<100.*(nPassed_/(double)nTot_)<<"%\n";
  	std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++"
	   	<<"\n"<<moduleLabel_<<"(ElectronIdSelector) SUMMARY:\n"<<ss.str()
	   	<<"++++++++++++++++++++++++++++++++++++++++++++++++++"
	   	<< std::endl;
}

////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////
typedef ElectronIdSelector   			    PATElectronIdSelector;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATElectronIdSelector);
