import FWCore.ParameterSet.Config as cms

displacedAssocToTracks = cms.EDProducer("DisplacedAssocToTracks",
        jetTracksAssociation = cms.untracked.InputTag("displacedAk4JetTracksAssociatorAtVertex","",""),        
        genParticleTag = cms.untracked.InputTag("genParticles","","RECO"),        
        outputLabel = cms.untracked.string('displacedAssocToTracks'),
        jetPtCut = cms.untracked.double(40.0),
        debug = cms.untracked.int32(0)
)
