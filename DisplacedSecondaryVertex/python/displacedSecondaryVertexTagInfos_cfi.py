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
	constraint = cms.string("None"), #displaced  -> None
	trackIPTagInfos = cms.InputTag("displacedLifetimeTagInfos"), #displaced 
	minimumTrackWeight = cms.double(0.5),
	usePVError = cms.bool(True), #displacd True -> True
	trackSort = cms.string('sip2dSig'), #displaced sip3dSig -> sip2dSig
        beamSpotTag = cms.InputTag('offlineBeamSpot'),                                        
        useExternalSV       = cms.bool(False),
        extSVCollection     = cms.InputTag('secondaryVertices'),
        extSVDeltaRToJet    = cms.double(0.3)

)
