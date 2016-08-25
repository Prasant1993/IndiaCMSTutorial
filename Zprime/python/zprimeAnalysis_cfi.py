import FWCore.ParameterSet.Config as cms

zpdimuon = cms.EDAnalyzer('ZPrime',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
)
