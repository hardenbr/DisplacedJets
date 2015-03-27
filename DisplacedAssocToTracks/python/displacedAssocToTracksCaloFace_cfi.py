import FWCore.ParameterSet.Config as cms

displacedAssocToTracksCaloFace = cms.EDProducer("DisplacedAssocToTracks",
        jetTracksAssociation = cms.untracked.InputTag("ak5JetTracksAssociatorAtCaloFace","","ANA"),        
        genParticleTag = cms.untracked.InputTag("genParticles","","RECO"),        
        outputLabel = cms.untracked.string('displacedAssocToTracksCaloFace'),
        jetPtCut = cms.untracked.double(40.0),
        debug = cms.untracked.int32(100)
)
