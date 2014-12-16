import FWCore.ParameterSet.Config as cms

#impact parameter and track counting
from DisplacedJets.DisplacedImpactParameter.impactParameter_cff import *

#re-produce primary vertices, lifetime, and secondary vertex tags
from DisplacedJets.DisplacedSecondaryVertex.displacedSecondaryVertex_cff import *

from DisplacedJets.DisplacedTagsToVertices.displacedTagsToVertices_cff import *

djtagging = cms.Sequence( displacedImpactParameterTagInfos + 
                          trackCountingDJTags + 
                          displacedOfflinePrimaryVertices + 
                          displacedLifetimeTagInfos + 
                          displacedSecondaryVertexTagInfos +
                          displacedTagsToVertices ) 
