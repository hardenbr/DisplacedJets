import FWCore.ParameterSet.Config as cms

#jet track associator at calo face
from DisplacedJets.DisplacedJetAssociationProducers.ak5JTA_cff import *

#impact parameter and track counting
from DisplacedJets.DisplacedImpactParameter.impactParameter_cff import *

#re-produce primary vertices, lifetime, and secondary vertex tags
#from DisplacedJets.DisplacedSecondaryVertex.displacedSecondaryVertex_cff import *
from DisplacedJets.DisplacedSecondaryVertexNoPV.displacedSecondaryVertex_cff import *

#produce the vertex collections with matching tracks for showing in event displays
from DisplacedJets.DisplacedTagsToVertices.displacedTagsToVertices_cff import *
#use the adaptive vertex finder
from DisplacedJets.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexing_cff import *

djtagging = cms.Sequence( #track matching for ak5 jets 
                          ak5JTA_noPF +                          
                          displacedImpactParameterTagInfos + 
                          trackCountingDJTags + 
                          displacedOfflinePrimaryVertices + 
                          #vertex matched tracks
                          displacedLifetimeTagInfos + 
                          displacedSecondaryVertexTagInfos +
                          displacedTagsToVerticesNoPV +
                          displacedTagsToVertices +
                          #calo face matched trackssequence
                          displacedLifetimeTagInfosCaloFace + 
                          displacedSecondaryVertexTagInfosNoPVCaloFace +
                          displacedSecondaryVertexTagInfosCaloFace +
                          displacedTagsToVerticesNoPVCaloFace +
                          #inclusive vertexing
                          displacedInclusiveVertexing 
) 
