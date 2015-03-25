import FWCore.ParameterSet.Config as cms

j2tParametersCALO = cms.PSet(
    tracks = cms.InputTag("generalTracks"),
    trackQuality = cms.string("any"), #djet goodIterative -> any
    extrapolations = cms.InputTag("trackExtrapolator"),
    coneSize = cms.double(0.5) #djet .4 -> .5
)

