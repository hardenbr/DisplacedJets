import FWCore.ParameterSet.Config as cms

displacedTrackVertexArbitrator = cms.EDProducer("TrackVertexArbitrator",
       beamSpot = cms.InputTag("offlineBeamSpot"),
       primaryVertices = cms.InputTag("offlinePrimaryVertices"),
       tracks = cms.InputTag("generalTracks"),
       secondaryVertices = cms.InputTag("displacedVertexMerger"),
       dLenFraction = cms.double(0.333),
       dRCut = cms.double(0.4),
       distCut = cms.double(0.04),
       sigCut = cms.double(5),
       fitterSigmacut =  cms.double(3),
       fitterTini = cms.double(256),
       fitterRatio = cms.double(0.25),
       trackMinLayers = cms.int32(0), #dijet -> 0
       trackMinPt = cms.double(1.0), #dijet .4 -> 1
       trackMinPixels = cms.int32(0) #djet 1 -> 0

)


