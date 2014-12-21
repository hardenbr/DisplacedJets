import FWCore.ParameterSet.Config as cms

from RecoBTag.ImpactParameter.variableJTA_cfi import *

vertexTrackSelectionBlock = cms.PSet(
	trackSelection = cms.PSet(
                variableJTAPars,
		totalHitsMin = cms.uint32(0), #displaced 8->0
		jetDeltaRMax = cms.double(0.4), #displaced .3 -> .4
		qualityClass = cms.string('any'),
		pixelHitsMin = cms.uint32(0), #dispalced 2 -> 0 no pixel hits
		maxDistToAxis = cms.double(.2), #displaced .2 -> 10
		maxDecayLen = cms.double(99999.9),
		sip3dSigMin = cms.double(-99999.9),
		sip3dSigMax = cms.double(99999.9),
		sip2dValMax = cms.double(99999.9),
		ptMin = cms.double(1.0),
		sip2dSigMax = cms.double(99999.9),
		sip2dSigMin = cms.double(-99999.9),
		sip3dValMax = cms.double(99999.9),
		sip3dValMin = cms.double(-99999.9),
		sip2dValMin = cms.double(-99999.9),
		normChi2Max = cms.double(99999.9),
                useVariableJTA = cms.bool(False) 
	)
)
