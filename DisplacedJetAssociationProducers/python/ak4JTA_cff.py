import FWCore.ParameterSet.Config as cms

from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import * ##propagator

from RecoJets.JetAssociationProducers.j2tParametersCALO_cfi import *
from DisplacedJets.DisplacedJetAssociationProducers.j2tParametersVX_cfi import *

displacedAk4JetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
    j2tParametersVX,
    jets = cms.InputTag("ak4CaloJetsL1FastL2L3")
)


displacedAk4JetTracksAssociatorAtVertexRegionalIter0124 = cms.EDProducer("JetTracksAssociatorAtVertex",
    j2tParametersVXRegionalIter0124,
    jets = cms.InputTag("ak4CaloJetsL1FastL2L3")
)

displacedAk4JetTracksAssociatorAtVertexRegionalIter012 = cms.EDProducer("JetTracksAssociatorAtVertex",
    j2tParametersVXRegionalIter012,
    jets = cms.InputTag("ak4CaloJetsL1FastL2L3")
)

displacedAk4JetTracksAssociatorAtVertexRegionalIter4 = cms.EDProducer("JetTracksAssociatorAtVertex",
    j2tParametersVXRegionalIter4,
    jets = cms.InputTag("ak4CaloJetsL1FastL2L3")
)

# ak4JetTracksAssociatorAtVertexPF = cms.EDProducer("JetTracksAssociatorAtVertex",
#     j2tParametersVX,
#     jets = cms.InputTag("ak4PFJetsCHS")
# )


# ak4JetTracksAssociatorExplicit = cms.EDProducer("JetTracksAssociatorExplicit",
#     j2tParametersVX,
#     jets = cms.InputTag("ak4PFJetsCHS")
# )

displacedAk4JetTracksAssociatorAtCaloFace = cms.EDProducer("JetTracksAssociatorAtCaloFace",
    j2tParametersCALO,
    jets = cms.InputTag("ak4CaloJetsL1FastL2L3")
)

# ak4JetExtender = cms.EDProducer("JetExtender",
#     jets = cms.InputTag("ak4CaloJets"),
#     jet2TracksAtCALO = cms.InputTag("ak4JetTracksAssociatorAtCaloFace"),
#     jet2TracksAtVX = cms.InputTag("ak4JetTracksAssociatorAtVertex"),
#     coneSize = cms.double(0.4)
# )

# ak4JTA = cms.Sequence(ak4JetTracksAssociatorAtVertexPF*
#                       ak4JetTracksAssociatorAtVertex*
#                       ak4JetTracksAssociatorAtCaloFace*ak4JetExtender)


ak4JTA_noPF = cms.Sequence(displacedAk4JetTracksAssociatorAtVertex*
                      displacedAk4JetTracksAssociatorAtCaloFace)#*ak4JetExtender)


#ak4JTAExplicit = cms.Sequence(ak4JetTracksAssociatorExplicit)
