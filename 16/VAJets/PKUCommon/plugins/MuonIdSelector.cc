/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 *
 *
 * Authors:
 *
 *   PKU 
 *
 * Description:
 *   - Selects "loose" and "tight" muons needed for V-boson analysis.
 *   - Saves collection of the reference vectors of muons passing the 
 *     required muon ID.
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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include <memory>
#include <vector>
#include <sstream>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////
class MuonIdSelector : public edm::EDProducer
{
public:
  	// construction/destruction
  	MuonIdSelector(const edm::ParameterSet& iConfig);
  	virtual ~MuonIdSelector();
  
  	// member functions
  	void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
  	void endJob();

private:  
  	// member data
  	// edm::InputTag  src_;
  	std::string    	moduleLabel_;
  	std::string    	idLabel_;  
  	bool           	applyTightID_;
  	bool           	applyLooseID_;
	bool			applyFakeID_;

  	unsigned int nTot_;
  	unsigned int nPassed_;
  	edm::EDGetTokenT<pat::MuonCollection> MuonToken_;
  	edm::EDGetTokenT<reco::VertexCollection> VertexToken_;
};

////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
MuonIdSelector::MuonIdSelector(const edm::ParameterSet& iConfig)
  	: moduleLabel_(iConfig.getParameter<std::string>   ("@module_label"))
  	, idLabel_(iConfig.existsAs<std::string>("idLabel") ? iConfig.getParameter<std::string>("idLabel") : "loose")
  	, nTot_(0)
  	, nPassed_(0)
  	, MuonToken_ (consumes<pat::MuonCollection> (iConfig.getParameter<edm::InputTag>( "src" ) ) ) 
  	, VertexToken_ (consumes<reco::VertexCollection> (iConfig.getParameter<edm::InputTag>( "vertex" ) ) )

{
  	produces<std::vector<pat::Muon> >();

  	/// ------- Decode the ID criteria --------
  	applyTightID_ 	= false;
  	applyLooseID_ 	= false;
    applyFakeID_ 	= false;

  	if( (idLabel_.compare("tight")==0) || 
      	(idLabel_.compare("Tight")==0) || 
      	(idLabel_.compare("TIGHT")==0) ||
      	(idLabel_.compare("WP70")==0) ||
      	(idLabel_.compare("wp70")==0) )  
 	applyTightID_ = true;

  	else if( (idLabel_.compare("loose")==0) || 
      	(idLabel_.compare("Loose")==0) || 
      	(idLabel_.compare("LOOSE")==0) ||
      	(idLabel_.compare("WP90")==0) ||
      	(idLabel_.compare("wp90")==0) )  
	applyLooseID_ = true;

    else if( (idLabel_.compare("Fake")==0) ||
        (idLabel_.compare("fake")==0) ||
        (idLabel_.compare("FAKE")==0) )
    applyFakeID_ = true;

}

//______________________________________________________________________________
MuonIdSelector::~MuonIdSelector(){}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////
 
//______________________________________________________________________________
void MuonIdSelector::produce(edm::Event& iEvent,const edm::EventSetup& iSetup)
{

  	/////// Pileup density "rho" in the event from fastJet pileup calculation ///// 
  	edm::Handle<reco::VertexCollection> vtxs;
  	iEvent.getByToken(VertexToken_, vtxs);
  	std::unique_ptr<std::vector<pat::Muon> > passingMuons(new std::vector<pat::Muon >);

  	edm::Handle<pat::MuonCollection > muons;
  	iEvent.getByToken(MuonToken_, muons);  

  	bool* isPassing = new bool[muons->size()];

  	for(unsigned int iMu=0; iMu<muons->size(); iMu++) { 

    	isPassing[iMu]=false;
    	const pat::Muon& mu1 = muons->at(iMu);
 
    	bool isTight  = false;  
    	bool isLoose  = false; 
		bool isFake	  = false;

		//tight ID
		if (mu1.pt()>20 && fabs(mu1.eta())<2.4 &&mu1.passed(reco::Muon::CutBasedIdTight|reco::Muon::PFIsoTight)) isTight = true;

		//fake muon ID
		if (mu1.pt()>20 && fabs(mu1.eta())<2.4 &&mu1.passed(reco::Muon::CutBasedIdTight|reco::Muon::PFIsoVeryLoose) && !mu1.passed(reco::Muon::CutBasedIdTight|reco::Muon::PFIsoTight)) isFake = true;

		//loose ID
		if (mu1.pt()>10 && fabs(mu1.eta())<2.4 &&mu1.passed(reco::Muon::CutBasedIdLoose|reco::Muon::PFIsoLoose)) isLoose = true;

    	/// ------- Finally apply selection --------
    	if(applyTightID_ && isTight)   	isPassing[iMu]= true;
    	if(applyLooseID_ && isLoose)   	isPassing[iMu]= true;
        if(applyFakeID_ && isFake)   	isPassing[iMu]= true;
     
 	}
 
	for (unsigned int iMuon = 0; iMuon < muons -> size(); iMuon ++){     
		if(isPassing[iMuon]) passingMuons->push_back( muons -> at(iMuon) );  
	}

  	nTot_  +=muons->size();
  	nPassed_+=passingMuons->size();
  	delete [] isPassing;  
  	iEvent.put(std::move(passingMuons));

}

//______________________________________________________________________________
void MuonIdSelector::endJob()
{
  	std::stringstream ss;
  	ss<<"nTot="<<nTot_<<" nPassed="<<nPassed_
    <<" effPassed="<<100.*(nPassed_/(double)nTot_)<<"%\n";
  	std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++"
	   	<<"\n"<<moduleLabel_<<"(MuonIdSelector) SUMMARY:\n"<<ss.str()
	   	<<"++++++++++++++++++++++++++++++++++++++++++++++++++"
	   	<< std::endl;

}

////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////
typedef MuonIdSelector    PATMuonIdSelector;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATMuonIdSelector);
