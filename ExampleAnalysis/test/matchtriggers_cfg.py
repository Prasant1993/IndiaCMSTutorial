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
	bits = cms.InputTag("TriggerResults","","HLT"),
	objects = cms.InputTag("selectedPatTrigger"),
	hltPath = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4"),
	hltFilter = cms.string("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2"),
	DeltaR = cms.double(0.1),
)	
process.p = cms.Path(process.matchtriggerObject)
