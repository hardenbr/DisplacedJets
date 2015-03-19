import FWCore.ParameterSet.Config as cms

vertexSelectionBlock = cms.PSet(
	vertexSelection = cms.PSet(
		sortCriterium = cms.string('dist2dSignificance') #displaced dist3dError -> dist2dSignificance
	)
)
