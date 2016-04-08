import FWCore.ParameterSet.Config as cms

j2tParametersVX = cms.PSet(
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.4),
    useAssigned = cms.bool(False),
    pvSrc = cms.InputTag("offlinePrimaryVerticesWithBS") #only used if useAssigned
)

j2tParametersVXRegionalIter0124 = cms.PSet(
    tracks = cms.InputTag("hltIter4MergedWithIter012DisplacedJets"),
    coneSize = cms.double(0.4),
    useAssigned = cms.bool(False),
    pvSrc = cms.InputTag("offlinePrimaryVerticesWithBS") #only used if useAssigned
)

j2tParametersVXRegionalIter012 = cms.PSet(
    tracks = cms.InputTag("hltIter2MergedForBTag"),
    coneSize = cms.double(0.4),
    useAssigned = cms.bool(False),
    pvSrc = cms.InputTag("offlinePrimaryVerticesWithBS") #only used if useAssigned
)

j2tParametersVXRegionalIter4 = cms.PSet(
    tracks = cms.InputTag("hltDisplacedhltIter4PFlowTrackSelectionHighPurity"),
    coneSize = cms.double(0.4),
    useAssigned = cms.bool(False),
    pvSrc = cms.InputTag("offlinePrimaryVerticesWithBS") #only used if useAssigned
)


