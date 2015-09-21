import FWCore.ParameterSet.Config as cms

# from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import * ##propagator

# from RecoJets.JetAssociationProducers.j2tParametersCALO_cfi import *
# from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import *

displacedAk4JetTracksAssociatorAtInnerHit = cms.EDProducer("JetTracksAssociatorAtInnerHitProducer",
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.4),
    jets = cms.InputTag("ak4CaloJetsL2L3")
)

ak4JTA_InnerHit = cms.Sequence(displacedAk4JetTracksAssociatorAtInnerHit)
