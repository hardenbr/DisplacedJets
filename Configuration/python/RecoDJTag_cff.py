import FWCore.ParameterSet.Config as cms

#jet track associator at calo face
from DisplacedJets.DisplacedJetAssociationProducers.ak5JTA_cff import *

#impact parameter and track counting
from DisplacedJets.DisplacedImpactParameter.impactParameter_cff import *

#re-produce primary vertices, lifetime, and secondary vertex tags
from DisplacedJets.DisplacedSecondaryVertex.displacedSecondaryVertex_cff import *
from DisplacedJets.DisplacedSecondaryVertexNoPV.displacedSecondaryVertex_cff import *

#produce the vertex collections with matching tracks for showing in event displays
from DisplacedJets.DisplacedTagsToVertices.displacedTagsToVertices_cff import *

#use the adaptive vertex finder
from DisplacedJets.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexing_cff import *

#saving only the tracks matched to certain high pt jets
from DisplacedJets.DisplacedAssocToTracks.displacedAssocToTracks_cff import *

djtagging = cms.Sequence( #track matching for ak5 jets 
                          ak5JTA_noPF +                 
                          #make the track collections for tracks matched to the jets
                          displacedAssocToTracks + 
                          displacedAssocToTracksCaloFace + 
                          #impact parameter info
                          displacedImpactParameterTagInfos + 
                          trackCountingDJTags + 
                          #unconstrained PV collection
                          displacedOfflinePrimaryVertices + 
                          #vertex matched tracks
                          displacedLifetimeTagInfos + 
                          displacedSecondaryVertexTagInfos +
                          displacedSecondaryVertexTagInfosNoPV + #noPV
                          displacedTagsToVerticesNoPV + #noPV
                          displacedTagsToVertices +
                          #calo face matched trackssequence
                          displacedLifetimeTagInfosCaloFace + 
                          displacedSecondaryVertexTagInfosNoPVCaloFace + #noPV 
                          displacedTagsToVerticesNoPVCaloFace + #noPV 
                          displacedSecondaryVertexTagInfosCaloFace +
                          displacedTagsToVerticesCaloFace +
                          #inclusive vertexing
                          displacedInclusiveVertexing 
) 
