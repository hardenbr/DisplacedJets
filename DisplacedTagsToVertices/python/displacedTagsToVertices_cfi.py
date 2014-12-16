import FWCore.ParameterSet.Config as cms

displacedTagsToVertices = cms.EDProducer("DisplacedTagsToVertices",
        secondaryVertexTagInfo = cms.untracked.InputTag("displacedSecondaryVertexTagInfos","","ANA"),        
        outputLabel = cms.untracked.string('displacedSecondaryVertices'),
        jetPtCut = cms.untracked.double(40.0)                                                                                  
)
