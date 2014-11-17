import FWCore.ParameterSet.Config as cms

isMC = True


process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

#geometry and global tag sequences
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MCRUN2_72_V1A::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/h/hardenbr/QCD_Pt-50to80_Tune4C_13TeV_pythia8_AOD.root'
    )
)

#calculate the track IP from the tags
process.trackIPFromJetTracks = cms.EDProducer( "TrackIPProducer",
    maximumTransverseImpactParameter = cms.double( 0.1 ),
    minimumNumberOfHits = cms.int32( 8 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    primaryVertex = cms.InputTag( "hltFastPVPixelVertices" ),
    maximumLongitudinalImpactParameter = cms.double( 0.1 ),
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    jetTracks = cms.InputTag( "offlinePrimaryVerticies","","RECO" ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    jetDirectionUsingTracks = cms.bool( False ),
    computeProbabilities = cms.bool( False ),
    useTrackQuality = cms.bool( False ),
    maximumChiSquared = cms.double( 20.0 )
)



#configure the analyzer 
process.analyzer = cms.EDAnalyzer('TrackAnalyzer')

#output configuration
process.analyzer.outputFileName =  cms.untracked.string('trackOutput.root')
process.analyzer.isMC =  cms.untracked.bool(isMC)

#tags
process.analyzer.generalTracks =  cms.untracked.InputTag('generalTracks','','RECO')
process.analyzer.ak5CaloJets =  cms.untracked.InputTag('ak5CaloJets','','RECO')

#mctags
if isMC:
    process.analyzer.ak5GenJets =  cms.untracked.InputTag('ak5GenJets','','SIM')
    process.analyzer.genMetCalo =  cms.untracked.InputTag('genMetCalo','','SIM')
    process.analyzer.genParticles =  cms.untracked.InputTag('genParticles','','SIM')



#run the anlayzer
process.p = cms.Path(process.analyzer)
