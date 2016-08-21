import FWCore.ParameterSet.Config as cms

process = cms.Process("DemoAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("IndiaCMSTutorial.ExampleAnalysis.ExampleAnalysis_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/20000/00071E92-6F55-E611-B68C-0025905A6066.root',
        '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/20000/00684732-3155-E611-B794-0CC47A4D7640.root',
        '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/20000/022C35C6-5C55-E611-8926-0025905A60B8.root',
        '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/20000/04407AAB-5454-E611-888E-24BE05CECDD1.root',
        #'/store/mc/RunIISpring16MiniAODv2/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/56FF84F6-D84A-E611-A75F-20CF3027A5D8.root',
        #'/store/mc/RunIISpring16MiniAODv2/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/F2343506-D94A-E611-BC5D-20CF3019DF0B.root'
    )
)

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('demo.root')
                                   )

process.p = cms.Path(process.exampleminiAOD)
