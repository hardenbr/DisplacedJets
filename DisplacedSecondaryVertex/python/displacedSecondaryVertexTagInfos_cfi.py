import FWCore.ParameterSet.Config as cms

from DisplacedJets.DisplacedSecondaryVertex.displacedVertexTrackSelection_cfi import *
from DisplacedJets.DisplacedSecondaryVertex.displacedVertexReco_cfi import *
from DisplacedJets.DisplacedSecondaryVertex.displacedVertexCuts_cfi import *
from DisplacedJets.DisplacedSecondaryVertex.displacedVertexSelection_cfi import *

displacedSecondaryVertexTagInfos = cms.EDProducer("DisplacedSecondaryVertexProducer",
	vertexTrackSelectionBlock,
        vertexSelectionBlock,
	vertexCutsBlock,
	vertexRecoBlock,
	constraint = cms.string("None"), #dijet  -> None
	trackIPTagInfos = cms.InputTag("displacedLifetimeTagInfos"), #dijet 
                                                  minimumTrackWeight = cms.double(.5), #dijet .5 -> 0
                                                  usePVError = cms.bool(True), #displacd True -> True
                                                  trackSort = cms.string('sip2dSig'), #dijet sip3dSig -> sip2dSig
                                                  beamSpotTag = cms.InputTag('offlineBeamSpot'),                                        
                                                  useExternalSV       = cms.bool(False), #djet True
                                                  extSVCollection     = cms.InputTag('displacedOfflinePrimaryVertices'), #djet secondaryVertices -> offlineDisplacedVertices 
                                                  extSVDeltaRToJet    = cms.double(.7), 
)
