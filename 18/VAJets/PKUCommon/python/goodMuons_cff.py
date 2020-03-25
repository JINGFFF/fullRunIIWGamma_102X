import FWCore.ParameterSet.Config as cms

tightMuIdLabel = "tight"
looseMuIdLabel = "loose"
fakeMuIdLabel  = "fake"


goodMuons = cms.EDProducer("PATMuonIdSelector",
    src = cms.InputTag( "slimmedMuons" ),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    idLabel = cms.string(tightMuIdLabel)
)   

looseMuons = cms.EDProducer("PATMuonIdSelector",
    src = cms.InputTag( "slimmedMuons" ),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    idLabel = cms.string(looseMuIdLabel)
)   

muSequence = cms.Sequence(goodMuons+looseMuons)
