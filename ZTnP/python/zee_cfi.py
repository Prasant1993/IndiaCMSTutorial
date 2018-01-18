import FWCore.ParameterSet.Config as cms

zeetnp = cms.EDAnalyzer('ZeeTnP',
    electrons = cms.InputTag("slimmedElectrons"),
)
