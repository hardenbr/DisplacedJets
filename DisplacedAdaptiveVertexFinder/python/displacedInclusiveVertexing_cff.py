import FWCore.ParameterSet.Config as cms

from DisplacedJets.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexFinder_cfi import *
from DisplacedJets.DisplacedAdaptiveVertexFinder.displacedVertexMerger_cfi import *
from DisplacedJets.DisplacedAdaptiveVertexFinder.displacedTrackVertexArbitrator_cfi import *

dispalcedInclusiveSecondaryVertices = vertexMerger.clone()
displacedInclusiveSecondaryVertices.secondaryVertices = cms.InputTag("trackVertexArbitrator")
displacedInclusiveSecondaryVertices.maxFraction = 0.2
displacedInclusiveSecondaryVertices.minSignificance = 10.

displacedInclusiveVertexing = cms.Sequence(displacedInclusiveVertexFinder * displacedVertexMerger * displacedTrackVertexArbitrator * displacedInclusiveSecondaryVertices)

