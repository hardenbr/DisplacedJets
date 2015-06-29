import FWCore.ParameterSet.Config as cms

#add the pat corrected met
from PhysicsTools.PatAlgos.recoLayer0.metCorrections_cff import *
from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import *

makePatMETs = cms.Sequence(
    patMETCorrections *
    patMETs
    )


