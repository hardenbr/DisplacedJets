import FWCore.ParameterSet.Config as cms

from RecoBTag.ImpactParameter.variableJTA_cfi import *

vertexTrackSelectionBlock = cms.PSet(
	trackSelection = cms.PSet(
                variableJTAPars,
		totalHitsMin = cms.uint32(0), #displaced 8->0
		jetDeltaRMax = cms.double(0.5), #displaced .3 -> .5
		qualityClass = cms.string('any'),
		pixelHitsMin = cms.uint32(0), #dispalced 2 -> 0 no pixel hits
		maxDistToAxis = cms.double(99.0), #displaced .2 -> 99
		maxDecayLen = cms.double(99999.9),
		sip3dSigMin = cms.double(-99999.0), 
		sip3dSigMax = cms.double(99999.9),
		sip2dValMax = cms.double(99999.9),
		ptMin = cms.double(1.0),
		sip2dSigMax = cms.double(99999.9),
		sip2dSigMin = cms.double(5.0), #displaced -9999 -> 5.0
		sip3dValMax = cms.double(99999.9),
		sip3dValMin = cms.double(-99999.0), 
		sip2dValMin = cms.double(.1), #displaced -99999->.1 
		normChi2Max = cms.double(99999.9),
                useVariableJTA = cms.bool(False) 
	)
)
