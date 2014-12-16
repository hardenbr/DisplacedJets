import FWCore.ParameterSet.Config as cms

displacedImpactParameterTagInfos = cms.EDProducer("TrackIPProducer",
    jetTracks = cms.InputTag("ak5JetTracksAssociatorAtVertex"), #displaced (using calo not PF)
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    computeProbabilities = cms.bool(True),
    computeGhostTrack = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    minimumNumberOfPixelHits = cms.int32(0), #displaced
    minimumNumberOfHits = cms.int32(0), #displaced 8 -> 0
    maximumTransverseImpactParameter = cms.double(999999999.0), #displaced
    minimumTransverseMomentum = cms.double(1.0),
    maximumChiSquared = cms.double(20.0),
    maximumLongitudinalImpactParameter = cms.double(9999999.0),
    jetDirectionUsingTracks = cms.bool(False),
    jetDirectionUsingGhostTrack = cms.bool(False),
    useTrackQuality = cms.bool(False),
)
