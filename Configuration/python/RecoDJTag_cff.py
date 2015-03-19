import FWCore.ParameterSet.Config as cms

#impact parameter and track counting
from DisplacedJets.DisplacedImpactParameter.impactParameter_cff import *

#re-produce primary vertices, lifetime, and secondary vertex tags
from DisplacedJets.DisplacedSecondaryVertex.displacedSecondaryVertex_cff import *

#produce the vertex collections with matching tracks for showing in event displays
from DisplacedJets.DisplacedTagsToVertices.displacedTagsToVertices_cff import *

#use the adaptive vertex finder
from DisplacedJets.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexing_cff import *

djtagging = cms.Sequence( displacedImpactParameterTagInfos + 
                          trackCountingDJTags + 
                          displacedOfflinePrimaryVertices + 
                          displacedLifetimeTagInfos + 
                          displacedSecondaryVertexTagInfos+
                          displacedTagsToVertices  ) 

#                          displacedInclusiveVertexing +
