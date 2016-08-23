import FWCore.ParameterSet.Config as cms

zmumutnp = cms.EDAnalyzer('ZmumuTnP',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
)
