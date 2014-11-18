import FWCore.ParameterSet.Config as cms

#define displaced jet sequences
from DisplacedJets.DisplacedImpactParameter.impactParameter_cff import *

djtagging = cms.Sequence(displacedImpactParameterTagInfos + trackCountingDJTags)
