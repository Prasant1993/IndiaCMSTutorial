import FWCore.ParameterSet.Config as cms

zmumutnp = cms.EDAnalyzer('ZmumuTnP',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
	bits = cms.InputTag("TriggerResults","","HLT"),
	objects = cms.InputTag("selectedPatTrigger"),
	hltPathName = cms.string("HLT_IsoMu22_v3"),
	hltFilterName= cms.string("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09"),
	DeltaR= cms.double(0.1),
)
