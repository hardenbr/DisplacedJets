import FWCore.ParameterSet.Config as cms

displacedTagsToVerticesNoPV = cms.EDProducer("DisplacedTagsToVertices",
        secondaryVertexTagInfo = cms.untracked.InputTag("displacedSecondaryVertexTagInfosNoPV","","ANA"),        
        genParticleTag = cms.untracked.InputTag("genParticles","","RECO"),        
        outputLabel = cms.untracked.string('displacedSecondaryVerticesNoPV'),
        jetPtCut = cms.untracked.double(40.0),
        isSignalMC = cms.untracked.bool(False),
        doGenMatch = cms.untracked.bool(False),
        debug = cms.untracked.int32(100)
)
