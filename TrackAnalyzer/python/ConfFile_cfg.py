import FWCore.ParameterSet.Config as cms

isMC = True


process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

#geometry and global tag sequences
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MCRUN2_72_V1A::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/h/hardenbr/QCD_Pt-50to80_Tune4C_13TeV_pythia8_AOD.root'
    )
)


# #custom ES producers for track computers
# process.DisplacedDijethltESPPromptTrackCountingESProducer = cms.ESProducer( "PromptTrackCountingESProducer",
#   maxImpactParameterSig = cms.double( 999999.0 ),
#   deltaR = cms.double( -1.0 ),
#   minimumImpactParameter = cms.double( -1.0 ),
#   maximumDecayLength = cms.double( 999999.0 ),
#   impactParameterType = cms.int32( 1 ),
#   trackQualityClass = cms.string( "any" ),
#   deltaRmin = cms.double( 0.0 ),
#   maxImpactParameter = cms.double( 0.1 ),
#   useSignedImpactParameterSig = cms.bool( True ),
#   maximumDistanceToJetAxis = cms.double( 999999.0 ),
#   nthTrack = cms.int32( -1 )
# )

# process.DisplacedDijethltESPTrackCounting2D1st = cms.ESProducer( "TrackCountingESProducer",
#   b_pT = cms.double( 0.3684 ),
#   deltaR = cms.double( -1.0 ),
#   minimumImpactParameter = cms.double( 0.05 ),
#   a_dR = cms.double( -0.001053 ),
#   min_pT = cms.double( 120.0 ),
#   maximumDistanceToJetAxis = cms.double( 9999999.0 ),
#   max_pT = cms.double( 500.0 ),
#   impactParameterType = cms.int32( 1 ),
#   trackQualityClass = cms.string( "any" ),
#   useVariableJTA = cms.bool( False ),
#   min_pT_dRcut = cms.double( 0.5 ),
#   max_pT_trackPTcut = cms.double( 3.0 ),
#   max_pT_dRcut = cms.double( 0.1 ),
#   b_dR = cms.double( 0.6263 ),
#   a_pT = cms.double( 0.005263 ),
#   maximumDecayLength = cms.double( 999999.0 ),
#   nthTrack = cms.int32( 1 ),
#   useSignedImpactParameterSig = cms.bool( False )
# )

# # second highest IP significance track computer
# process.DisplacedDijethltESPTrackCounting2D2nd= cms.ESProducer( "TrackCountingESProducer",
#   b_pT = cms.double( 0.3684 ),
#   deltaR = cms.double( -1.0 ),
#   minimumImpactParameter = cms.double( 0.05 ),
#   a_dR = cms.double( -0.001053 ),
#   min_pT = cms.double( 120.0 ),
#   maximumDistanceToJetAxis = cms.double( 9999999.0 ),
#   max_pT = cms.double( 500.0 ),
#   impactParameterType = cms.int32( 1 ),
#   trackQualityClass = cms.string( "any" ),
#   useVariableJTA = cms.bool( False ),
#   min_pT_dRcut = cms.double( 0.5 ),
#   max_pT_trackPTcut = cms.double( 3.0 ),
#   max_pT_dRcut = cms.double( 0.1 ),
#   b_dR = cms.double( 0.6263 ),
#   a_pT = cms.double( 0.005263 ),
#   maximumDecayLength = cms.double( 999999.0 ),
#   nthTrack = cms.int32( 2 ),
#   useSignedImpactParameterSig = cms.bool( False )
# )


# #calculate the track IPs from the jet tracks associator
# process.trackIPsFromJetTracks = cms.EDProducer( "TrackIPProducer",
#     maximumTransverseImpactParameter = cms.double( 0.1 ),
#     minimumNumberOfHits = cms.int32( 8 ),
#     minimumTransverseMomentum = cms.double( 1.0 ),
#     primaryVertex =  cms.InputTag( "offlinePrimaryVerticies","","RECO" ),
#     maximumLongitudinalImpactParameter = cms.double( 0.1 ),
#     computeGhostTrack = cms.bool( False ),
#     ghostTrackPriorDeltaR = cms.double( 0.03 ),
#     jetTracks = cms.InputTag( "ak5JetTracksAssociatorAtVertex","","RECO" ),
#     jetDirectionUsingGhostTrack = cms.bool( False ),
#     minimumNumberOfPixelHits = cms.int32( 2 ),
#     jetDirectionUsingTracks = cms.bool( False ),
#     computeProbabilities = cms.bool( False ),
#     useTrackQuality = cms.bool( False ),
#     maximumChiSquared = cms.double( 20.0 )
# )

# #build the jet tags from the track IP
# process.jetTagProducer = cms.EDProducer( "JetTagProducer",
#     jetTagComputer = cms.string( "DisplacedDijethltESPPromptTrackCountingESProducer" ),
#     tagInfos = cms.VInputTag( 'trackIPFromJetTracks' )
# )


################################################################################################


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
