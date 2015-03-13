import FWCore.ParameterSet.Config as cms

trackCountingDJTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('displacedTrackCounting2D1st'),
    tagInfos = cms.VInputTag(cms.InputTag("displacedImpactParameterTagInfos"))
)
