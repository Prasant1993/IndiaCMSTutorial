import FWCore.ParameterSet.Config as cms

triggerFilter = cms.EDFilter("TriggerFilter",
  verbosity = cms.untracked.int32(0),
  l1InputTag = cms.untracked.InputTag('gtDigis'),
  hltInputTag = cms.untracked.InputTag('TriggerResults','','HLT2'),
  hltPathsOfInterest = cms.vstring ( 
    'HLT_Mu50' 
  )
)
