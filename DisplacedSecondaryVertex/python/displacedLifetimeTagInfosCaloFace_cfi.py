import FWCore.ParameterSet.Config as cms

displacedLifetimeTagInfosCaloFace = cms.EDProducer( "TrackIPProducer",
    maximumTransverseImpactParameter = cms.double( 999999.0 ), #displaced .2 -> infinity
    minimumNumberOfHits = cms.int32( 0 ), #displaced 8 -> 0
    minimumTransverseMomentum = cms.double( 1.0 ), 
    primaryVertex = cms.InputTag( 'offlinePrimaryVerticesWithBS'), #with BS or without BS?
    maximumLongitudinalImpactParameter = cms.double( 999999.0 ), #displaced 17 -> infinity
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    jetTracks = cms.InputTag( "displacedAk4JetTracksAssociatorAtCaloFace" ), #displaced  hltFastPixelBLifetimeL3Associator -> jetTracks
    jetDirectionUsingGhostTrack = cms.bool( False ),
    minimumNumberOfPixelHits = cms.int32( 0 ), #displaced 2 -> 0 AOD PROBLEMS
    jetDirectionUsingTracks = cms.bool( False ),
    computeProbabilities = cms.bool( False ),
    useTrackQuality = cms.bool( False ),
    maximumChiSquared = cms.double( 20.0 ) 
)
