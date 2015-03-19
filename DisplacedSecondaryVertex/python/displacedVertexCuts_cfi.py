
import FWCore.ParameterSet.Config as cms

vertexCutsBlock = cms.PSet(
	vertexCuts = cms.PSet(
		fracPV = cms.double(100), #displaced ? -> 100
		distSig3dMax = cms.double(99999.9),
		distVal2dMax = cms.double(99999.9), #displaced 2.5 -> infinity
		useTrackWeights = cms.bool(False), #displaced True -> False only matters for mass cut!
		maxDeltaRToJetAxis = cms.double(2.0), #displaced .5 -> 2
		v0Filter = cms.PSet(k0sMassWindow = cms.double(-9999)), #k short filter out? displaced .05 --> -infty  format: abs(vtx_mass - k0s mass) < masswindow)
		distSig2dMin = cms.double(0.0), #displaced 3.0 -> 0.0 
		multiplicityMin = cms.uint32(2), #displaced 2 -> 2
		massMax = cms.double(99999.9), # displaced 6.5 -> infity
		distSig2dMax = cms.double(99999.9),
		distVal3dMax = cms.double(99999.9),
		minimumTrackWeight = cms.double(0.0), #displaced .5 -> 0
		distVal3dMin = cms.double(-99999.9),
		distVal2dMin = cms.double(-99999.0), #displaced .01 -> -999999
		distSig3dMin = cms.double(-99999.9)
	)
)
