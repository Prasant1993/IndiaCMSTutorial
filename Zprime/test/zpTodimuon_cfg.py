import FWCore.ParameterSet.Config as cms

process = cms.Process("ZpAnalysis")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("IndiaCMSTutorial.Zprime.zprimeAnalysis_cfi")
process.load("IndiaCMSTutorial.Zprime.triggerFilter_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_v14'

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
	fileNames = cms.untracked.vstring(
        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/449D7C8F-C63B-E611-B26C-00304867FD7B.root',
        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/881A437D-A63B-E611-96C5-0025905745B8.root',
        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/9639BE37-C63B-E611-ADE5-02163E00F456.root',
        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/AC05067C-1C3B-E611-9DA6-001E67DFF4F6.root',
        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/B8A74EF0-383B-E611-B7E4-00259048AC10.root',
        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/E0D165F0-AB3B-E611-8BF7-2C600CAFEE50.root',
        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/F2301E48-C63B-E611-BEAF-001E67C7AC24.root',
    )
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('zPrimeDimuon.root')
)

#process.triggerFilter.hltInputTag = cms.untracked.InputTag('TriggerResults','','HLT')

# Primary Vertex Selector
process.selectedPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag('offlineSlimmedPrimaryVertices'),
  cut = cms.string("!isFake && abs(z) <= 24 && abs(position.Rho) <= 2"),
  filter = cms.bool(True)                                          
)

#process.p = cms.Path(process.triggerFilter*process.zpdimuon)
process.p = cms.Path(process.triggerFilter*process.selectedPrimaryVertices*process.zpdimuon)
