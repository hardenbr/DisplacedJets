import FWCore.ParameterSet.Config as cms

#jet track associator at calo face
from DisplacedJets.DisplacedJetAssociationProducers.ak4JTA_cff import *

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

#matching at the inner track hit
from DisplacedJets.DisplacedJetInnerHitTrackAssociator.ak4JTA_InnerHit_cff import *

#alpha calculation based on primary vertices
#from DisplacedJets.DisplacedJetVertexAssociation.displacedJetVertexAssociation_cff import *

djtagging = cms.Sequence( #track matching for ak5 jets 
     ak4JTA_noPF + # jet track matching at vertex and calo face                 
     #ak4JTA_InnerHit + #jet track matching at the inner hit                  
     #make the track collections for tracks matched to the jets
     displacedAssocToTracks + 
     displacedAssocToTracksCaloFace + 
#     displacedJetVertexAssociation + 
#     displacedAssocToTracksInnerHit + 
     #impact parameter info
     #displacedImpactParameterTagInfos + 
     #trackCountingDJTags + 
     #unconstrained PV collection
     #displacedOfflinePrimaryVertices + 
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
     displacedTagsToVerticesCaloFace+
     #inclusive vertexing
     displacedInclusiveVertexing 
)
 

regionalTrackAssocations = cms.Sequence( displacedAk4JetTracksAssociatorAtVertexRegionalIter4 + displacedAk4JetTracksAssociatorAtVertexRegionalIter012 + displacedAk4JetTracksAssociatorAtVertexRegionalIter0124)
