import FWCore.ParameterSet.Config as cms

displacedVertexMerger = cms.EDProducer("VertexMerger",
       secondaryVertices = cms.InputTag("displacedInclusiveVertexFinderJetMatchedTracksCaloFace"),
       maxFraction = cms.double(0.7), #djet .7 -> .5
       minSignificance = cms.double(2)) #djet  2 -> 1


