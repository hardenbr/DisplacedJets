import FWCore.ParameterSet.Config as cms

displacedTagsToVertices = cms.EDProducer("DisplacedTagsToVertices",
        secondaryVertexTagInfo = cms.untracked.InputTag("displacedSecondaryVertexTagInfos","","ANA"),        
        genParticleTag = cms.untracked.InputTag("genParticles","","RECO"),        
        outputLabel = cms.untracked.string('displacedSecondaryVertices'),
        jetPtCut = cms.untracked.double(40.0),
        isSignalMC = cms.untracked.bool(True),
        doGenMatch = cms.untracked.bool(True),
        debug = cms.untracked.int32(2)
)
