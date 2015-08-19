import FWCore.ParameterSet.Config as cms

displacedImpactParameterTagInfos = cms.EDProducer("TrackIPProducer",
    jetTracks = cms.InputTag("displacedAk4JetTracksAssociatorAtVertex"), #displaced 
    primaryVertex = cms.InputTag("offlinePrimaryVerticesWithBS"),
    computeProbabilities = cms.bool(True),
    computeGhostTrack = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    minimumNumberOfPixelHits = cms.int32(0), #displaced
    minimumNumberOfHits = cms.int32(0), #displaced 8 -> 0
    maximumTransverseImpactParameter = cms.double(999999999.0), #displaced
    minimumTransverseMomentum = cms.double(1.0),
    maximumChiSquared = cms.double(999999999.0), #displaced 20 -> infity
    maximumLongitudinalImpactParameter = cms.double(999999999.0),
    jetDirectionUsingTracks = cms.bool(True),
    jetDirectionUsingGhostTrack = cms.bool(False),
    useTrackQuality = cms.bool(False),
)
