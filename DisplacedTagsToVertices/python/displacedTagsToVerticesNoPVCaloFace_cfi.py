import FWCore.ParameterSet.Config as cms

displacedTagsToVerticesNoPVCaloFace = cms.EDProducer("DisplacedTagsToVertices",
        secondaryVertexTagInfo = cms.untracked.InputTag("displacedSecondaryVertexTagInfosNoPVCaloFace","","ANA"),        
        genParticleTag = cms.untracked.InputTag("genParticles","","RECO"),        
        outputLabel = cms.untracked.string('displacedSecondaryVerticesNoPVCaloFace'),
        jetPtCut = cms.untracked.double(70.0),
        isSignalMC = cms.untracked.bool(False),
        doGenMatch = cms.untracked.bool(False),
        debug = cms.untracked.int32(100)
)
