import FWCore.ParameterSet.Config as cms

zpdimuon = cms.EDAnalyzer('ZPrime',
    #vertices = cms.InputTag("selectedPrimaryVertices"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    packed = cms.InputTag("packedGenParticles"),
    pruned = cms.InputTag("prunedGenParticles")
)
