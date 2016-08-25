import FWCore.ParameterSet.Config as cms

process = cms.Process("ZpAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("IndiaCMSTutorial.ZTnP.zmumu_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
	fileNames = cms.untracked.vstring(
        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/10000/0C3755F5-173B-E611-8F95-0CC47A1DF7FE.root'
    )
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('zPrimeDimuon.root')
)

process.p = cms.Path(process.zpdimuon)
