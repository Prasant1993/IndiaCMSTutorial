import FWCore.ParameterSet.Config as cms
process = cms.Process("matchtriggerObject")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/10000/0C3755F5-173B-E611-8F95-0CC47A1DF7FE.root'
    )
)
process.matchtriggerObject=cms.EDAnalyzer("TriggerMatcher",
	electrons = cms.InputTag("slimmedElectrons"),
	muons = cms.InputTag("slimmedMuons"),
	bits = cms.InputTag("TriggerResults","","HLT2"),
	objects = cms.InputTag("selectedPatTrigger"),
	hltPath = cms.string("HLT_IsoMu22_v3"),
	hltFilter = cms.string("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09"),
	DeltaR = cms.double(0.1),
)	
process.p = cms.Path(process.matchtriggerObject)
