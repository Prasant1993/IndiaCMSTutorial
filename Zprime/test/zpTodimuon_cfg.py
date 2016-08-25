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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_v14'

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
	fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/work/s/sroychow/public/Indiacms/prod/CMSSW_8_0_12/src/step4.root'
        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/AC05067C-1C3B-E611-9DA6-001E67DFF4F6.root'
    )
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('zPrimeDimuon.root')
)

#process.triggerFilter.hltInputTag = cms.untracked.InputTag('TriggerResults','','HLT')
process.p = cms.Path(process.triggerFilter*process.zpdimuon)
