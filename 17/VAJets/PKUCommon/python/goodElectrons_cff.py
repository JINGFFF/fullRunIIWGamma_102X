import FWCore.ParameterSet.Config as cms

tightEleIdLabel = "tight"
mediumEleIdLabel = "medium"
looseEleIdLabel = "loose"
vetoEleIdLabel = "veto"
fakeEleIdLabel = "fake"

goodElectrons = cms.EDProducer("PATElectronIdSelector",
    src = cms.InputTag( "slimmedElectrons" ),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    idLabel = cms.string(tightEleIdLabel),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt")
)

vetoElectrons = cms.EDProducer("PATElectronIdSelector",
    src = cms.InputTag( "slimmedElectrons" ),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    idLabel = cms.string(vetoEleIdLabel),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt")
)


eleSequence = cms.Sequence(goodElectrons+vetoElectrons)
