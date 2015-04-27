import FWCore.ParameterSet.Config as cms

displacedJetSVAssociationIVF = cms.EDProducer("DisplacedJetSVAssociator", 
                                              caloJets          = cms.untracked.InputTag("ak4CaloJets","",""), 
                                              secondaryVertices = cms.untracked.InputTag("displacedInclusiveSecondaryVertices"),
                                              primaryVertices   = cms.untracked.InputTag("offlinePrimaryVertices"),
                                              outputLabel       = cms.untracked.string('displacedIVFJetAssoc'),
                                              jetPtCut          = cms.untracked.double(70.0),
                                              algoName          = cms.untracked.string("oneOverDeltaR"),
                                              debug             = cms.untracked.int32(1)
)
