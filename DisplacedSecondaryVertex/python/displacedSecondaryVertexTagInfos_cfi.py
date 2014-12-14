import FWCore.ParameterSet.Config as cms

from DisplacedJets.DisplacedSecondaryVertex.displacedVertexTrackSelection_cfi import *
from DisplacedJets.DisplacedSecondaryVertex.displacedVertexReco_cfi import *
from DisplacedJets.DisplacedSecondaryVertex.displacedVertexCuts_cfi import *
from DisplacedJets.DisplacedSecondaryVertex.displacedVertexSelection_cfi import *

secondaryVertexTagInfos = cms.EDProducer("SecondaryVertexProducer",
	vertexTrackSelectionBlock,
	vertexSelectionBlock,
	vertexCutsBlock,
	vertexRecoBlock,
	constraint = cms.string("BeamSpot"),
	trackIPTagInfos = cms.InputTag(""), #displaced NEEDS TO BE LIFETIME TAG INFOS
	minimumTrackWeight = cms.double(0.5),
	usePVError = cms.bool(True),
	trackSort = cms.string('sip2dSig'), #displaced 3d->2d
        beamSpotTag = cms.InputTag('offlineBeamSpot'),                                        
        useExternalSV       = cms.bool(False),
        extSVCollection     = cms.InputTag('secondaryVertices'),
        extSVDeltaRToJet    = cms.double(0.3)

)
