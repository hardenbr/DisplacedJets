import FWCore.ParameterSet.Config as cms

from DisplacedJets.DisplacedSecondaryVertex.displacedVertexTrackSelection_cfi import *
from DisplacedJets.DisplacedSecondaryVertex.displacedVertexReco_cfi import *
from DisplacedJets.DisplacedSecondaryVertex.displacedVertexCuts_cfi import *
from DisplacedJets.DisplacedSecondaryVertex.displacedVertexSelection_cfi import *

displacedSecondaryVertexTagInfos = cms.EDProducer("SecondaryVertexProducer",
	vertexTrackSelectionBlock,
	vertexSelectionBlock,
	vertexCutsBlock,
	vertexRecoBlock,
	constraint = cms.string("None"), #dijet  -> None
	trackIPTagInfos = cms.InputTag("displacedLifetimeTagInfos"), #dijet 
	minimumTrackWeight = cms.double(0.0), #dijet .5 -> 0
	usePVError = cms.bool(True), #displacd True -> True
	trackSort = cms.string('sip2dSig'), #dijet sip3dSig -> sip2dSig
        beamSpotTag = cms.InputTag('offlineBeamSpot'),                                        
        useExternalSV       = cms.bool(True), #djet True
        extSVCollection     = cms.InputTag('displacedOfflinePrimaryVertices'), #djet secondaryVertices -> offlineDisplacedVertices 
        extSVDeltaRToJet    = cms.double(.7), 
)
