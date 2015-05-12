import FWCore.ParameterSet.Config as cms

displacedTrackVertexArbitrator = cms.EDProducer("TrackVertexArbitrator",
       beamSpot = cms.InputTag("offlineBeamSpot"),
       primaryVertices = cms.InputTag("offlinePrimaryVertices"),
#       tracks = cms.InputTag("displacedAssocToTracksCaloFace","displacedAssocToTracksCaloFace","ANA"),
       tracks = cms.InputTag("generalTracks"),
       secondaryVertices = cms.InputTag("displacedInclusiveVertexFinder"),
       dLenFraction = cms.double(0.333), #djet .333 -> .2
       dRCut = cms.double(1.0), # djet .4 -> 1
       distCut = cms.double(0.1), #djet .04 -> .1
       sigCut = cms.double(10), #djet 5->10
       fitterSigmacut =  cms.double(3),
       fitterTini = cms.double(256),
       fitterRatio = cms.double(0.25),
       trackMinLayers = cms.int32(0), #dijet -> 0
       trackMinPt = cms.double(.4), 
       trackMinPixels = cms.int32(0) #djet 1 -> 0
)


