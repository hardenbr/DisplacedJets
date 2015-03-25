import FWCore.ParameterSet.Config as cms

displacedTagsToVerticesCaloFace = cms.EDProducer("DisplacedTagsToVertices",
        secondaryVertexTagInfo = cms.untracked.InputTag("displacedSecondaryVertexTagInfosCaloFace","","ANA"),        
        genParticleTag = cms.untracked.InputTag("genParticles","","RECO"),        
        outputLabel = cms.untracked.string('displacedSecondaryVerticesCaloFace'),
        jetPtCut = cms.untracked.double(100.0),
        isSignalMC = cms.untracked.bool(False),
        doGenMatch = cms.untracked.bool(False),
        debug = cms.untracked.int32(100)
)
