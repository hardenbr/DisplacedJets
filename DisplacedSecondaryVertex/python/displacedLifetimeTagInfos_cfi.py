import FWCore.ParameterSet.Config as cms

displacedLifetimeTagInfos = cms.EDProducer( "TrackIPProducer",
    maximumTransverseImpactParameter = cms.double( 1000.0 ), #displaced .2 -> infinity
    minimumNumberOfHits = cms.int32( 0 ), #displaced 8 -> 0
    minimumTransverseMomentum = cms.double( 1.0 ), 
    primaryVertex = cms.InputTag( 'offlinePrimaryVertices'), #with BS or without BS?
    maximumLongitudinalImpactParameter = cms.double( 17.0 ), #displaced 17 -> infinity
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    jetTracks = cms.InputTag( "ak5JetTracksAssociatorAtVertex" ), #displaced  hltFastPixelBLifetimeL3Associator -> jetTracks
    jetDirectionUsingGhostTrack = cms.bool( False ),
    minimumNumberOfPixelHits = cms.int32( 0 ), #displaced 2 -> 0 
    jetDirectionUsingTracks = cms.bool( False ),
    computeProbabilities = cms.bool( False ),
    useTrackQuality = cms.bool( False ),
    maximumChiSquared = cms.double( 20.0 )
)
