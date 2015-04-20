import FWCore.ParameterSet.Config as cms

displacedVertexMerger = cms.EDProducer("VertexMerger",
       secondaryVertices = cms.InputTag("displacedInclusiveVertexFinderJetMatchedTracksCaloFace"),
       maxFraction = cms.double(0.7),
       minSignificance = cms.double(2)) 


