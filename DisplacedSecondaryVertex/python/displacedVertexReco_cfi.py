import FWCore.ParameterSet.Config as cms

vertexRecoBlock = cms.PSet(
	vertexReco = cms.PSet(
		#seccut = cms.double(6.0), # djet 6-> 3
		#primcut = cms.double(1.0), # djet 1.8  -> 1.0
		#smoothing = cms.bool(True), #djet False -> True
		finder = cms.string('avr'),
		#minweight = cms.double(0.0), #djet .5 -> 0 
		#weightthreshold = cms.double(0.001)
	)
)
