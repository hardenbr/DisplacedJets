import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

#geometry and global tag sequences
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MCRUN2_72_V1A::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'root://xrootd.unl.edu//store/mc/Fall13dr/QCD_Pt-80to120_Tune4C_13TeV_pythia8/AODSIM/castor_tsg_PU40bx25_POSTLS162_V2-v1/00000/0020BEFF-7B84-E311-8519-90E6BA442F2B.root'
        'file:/afs/cern.ch/work/h/hardenbr/QCD_Pt-50to80_Tune4C_13TeV_pythia8_AOD.root'
#        'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/alca_ecalcalib/hardenbr/TEST_FILES/QCD80_120_RAW.root'
    )
)

process.analyzer = cms.EDAnalyzer('TrackAnalyzer')

process.p = cms.Path(process.analyzer)
process.analyzer.generalTracks =  cms.untracked.InputTag('generalTracks','','RECO')
