import FWCore.ParameterSet.Config as cms

displacedAssocToTracksInnerHit = cms.EDProducer("DisplacedAssocToTracks",
        jetTracksAssociation = cms.untracked.InputTag("displacedAk4JetTracksAssociatorAtCaloFace", "", ""),   
        genParticleTag = cms.untracked.InputTag("genParticles","","RECO"),        
        outputLabel = cms.untracked.string('displacedAssocToTracksInnerHit'),
        jetPtCut = cms.untracked.double(40.0),
        debug = cms.untracked.int32(0)
)
