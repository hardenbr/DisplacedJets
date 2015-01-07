import FWCore.ParameterSet.Config as cms

isSignalMC = False
isMC = True
doedm = False
nevents = -1

input_file = None
if isSignalMC:
    input_file = 'file:/afs/cern.ch/work/h/hardenbr/2015/DIJET/GEN_SIGNAL_TEST/dijet_700_300_ctau30.root'
else:
    input_file =  'file:/afs/cern.ch/work/h/hardenbr/QCD_470_600_AOD_40bx25.root'

process = cms.Process("ANA")
proc_label = "RECO"

# from RecoBTag.SoftLepton.softLepton_cff import *
# from RecoBTag.ImpactParameter.impactParameter_cff import *
# from RecoBTag.SecondaryVertex.secondaryVertex_cff import *
# from RecoBTau.JetTagComputer.combinedMVA_cff import *
# from RecoBTag.Configuration.RecoBTag_EventContent_cff import *
# process.load('RecoBTag/Configuration/RecoBTag_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")

#geometry and global tag sequences
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

#old global tag with wrong connect address

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'MCRUN2_72_V1A::All'

process.GlobalTag = cms.ESSource( "PoolDBESSource",
    globaltag = cms.string( "GR_H_V39" ),
    RefreshEachRun = cms.untracked.bool( True ),
    RefreshOpenIOVs = cms.untracked.bool( False ),
    toGet = cms.VPSet( 
      cms.PSet(  record = cms.string( "JetCorrectionsRecord" ),
        tag = cms.string( "JetCorrectorParametersCollection_HLT_V1_AK4Calo" ),
        connect = cms.untracked.string( "frontier://FrontierPrep/CMS_COND_PHYSICSTOOLS" ),
        label = cms.untracked.string( "AK4CaloHLT" )
      ),
      cms.PSet(  record = cms.string( "JetCorrectionsRecord" ),
        tag = cms.string( "JetCorrectorParametersCollection_HLT_trk1B_V1_AK4PF" ),
        connect = cms.untracked.string( "frontier://FrontierPrep/CMS_COND_PHYSICSTOOLS" ),
        label = cms.untracked.string( "AK4PFHLT" )
      )
    ),
    DBParameters = cms.PSet( 
      authenticationPath = cms.untracked.string( "." ),
      connectionRetrialTimeOut = cms.untracked.int32( 60 ),
      idleConnectionCleanupPeriod = cms.untracked.int32( 10 ),
      messageLevel = cms.untracked.int32( 0 ),
      enablePoolAutomaticCleanUp = cms.untracked.bool( False ),
      enableConnectionSharing = cms.untracked.bool( True ),
      enableReadOnlySessionOnUpdateConnection = cms.untracked.bool( False ),
      connectionTimeOut = cms.untracked.int32( 0 ),
      connectionRetrialPeriod = cms.untracked.int32( 10 )
    ),
    RefreshAlways = cms.untracked.bool( False ),
    connect = cms.string( "frontier://(proxyurl=http://localhost:3128)(serverurl=http://localhost:8000/FrontierOnProd)(serverurl=http://localhost:8000/FrontierOnProd)(retrieve-ziplevel=0)/CMS_CONDITIONS" ),
    ReconnectEachRun = cms.untracked.bool( True ),
    BlobStreamerName = cms.untracked.string( "TBufferBlobStreamingService" )
)
# override the GlobalTag, connection string and pfnPrefix
if 'GlobalTag' in process.__dict__:
    from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'MCRUN2_72_V1A')
    process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_CONDITIONS'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    for pset in process.GlobalTag.toGet.value():
        pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
    # fix for multi-run processing
    process.GlobalTag.RefreshEachRun = cms.untracked.bool( False )
    process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nevents))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        input_file
#        'file:/afs/cern.ch/work/h/hardenbr/QCD_Pt-50to80_Tune4C_13TeV_pythia8_AOD.root'
#        'file:/afs/cern.ch/work/h/hardenbr/HTo2LongLivedTo4L_MH_700_MFF_300_CTau30_TSG_PU40BX25_AODSIM.root'
#       'file:/afs/cern.ch/work/h/hardenbr/HTo2LongLivedTo4L_MH_700_MFF_300_CTau30_TSG_PU40BX25_AODSIM_v6.root'
#        'file:/afs/cern.ch/work/h/hardenbr/QCD_470_600_AOD_40bx25.root'
#        'file:/afs/cern.ch/work/h/hardenbr/TEST_FILES/HTo2LongLivedTo4L_MH_700_MFF_300_CTau30_TSG_PU40BX25_AODSIM_10ev.root'
    )
)



################################################################################################

#configure the analyzer 
process.analyzer = cms.EDAnalyzer('TrackAnalyzer')

#output configuration
if isSignalMC:
    process.analyzer.outputFileName = cms.untracked.string('signal.root')
else:
    process.analyzer.outputFileName = cms.untracked.string('qcd.root')
process.analyzer.isMC =  cms.untracked.bool(isMC)
process.analyzer.isSignalMC =  cms.untracked.bool(isSignalMC)

#tags
process.analyzer.generalTracks =  cms.untracked.InputTag('generalTracks', '', proc_label)
process.analyzer.ak5CaloJets =  cms.untracked.InputTag('ak5CaloJets', '', proc_label)
process.analyzer.genParticles =  cms.untracked.InputTag('genParticles', '', proc_label)

process.analyzer.trackIPTagInfoCollection = cms.untracked.InputTag('displacedImpactParameterTagInfos', '', 'ANA')
process.analyzer.secondaryVertexTagInfo = cms.untracked.InputTag('displacedSecondaryVertexTagInfos', '', 'ANA')
process.analyzer.lifetimeIPTagInfo = cms.untracked.InputTag('displacedLifetimeTagInfos', '', 'ANA')

#cuts
process.analyzer.jetPt = cms.untracked.double(40.0)
process.analyzer.jetEta = cms.untracked.double(2.0)

# Tags related to the monte carlo
if isMC:
    # Gen Information
    process.analyzer.ak5GenJets =  cms.untracked.InputTag('ak5GenJets', '', proc_label)
    process.analyzer.genMetCalo =  cms.untracked.InputTag('genMetCalo', '', proc_label)
    process.analyzer.genParticles =  cms.untracked.InputTag('genParticles', '', proc_label)

# process.edmoutput = cms.OutputModule( "PoolOutputModule",
#     fileName = cms.untracked.string( "edmoutput.root" ),
#     fastCloning = cms.untracked.bool( False ),
#     dataset = cms.untracked.PSet(
#         filterName = cms.untracked.string( "" ),
#         dataTier = cms.untracked.string( "RAW" )
#     ),
#                                       outputCommands = cms.untracked.vstring(
#         'keep *_impactParameterTagInfos_*_*',
#         'keep *_trackCountingHighEffBJetTags_*_*',
#         'keep *_trackCountingHighPurBJetTags_*_*',
#         'keep *_jetProbabilityBJetTags_*_*',
#         'keep *_jetBProbabilityBJetTags_*_*',
#         'keep *_secondaryVertexTagInfos_*_*',
#         'keep *_inclusiveSecondaryVertexFinderTagInfos_*_*',
#         'keep *_ghostTrackVertexTagInfos_*_*',
#         'keep *_simpleSecondaryVertexBJetTags_*_*',
#         'keep *_simpleSecondaryVertexHighEffBJetTags_*_*',
#         'keep *_simpleSecondaryVertexHighPurBJetTags_*_*',
#         'keep *_combinedSecondaryVertexBJetTags_*_*',
#         'keep *_combinedInclusiveSecondaryVertexV2BJetTags_*_*',
#         'keep *_ghostTrackBJetTags_*_*',
#         'keep *_softPFMuonsTagInfos_*_*',
#         'keep *_softPFElectronsTagInfos_*_*',
#         'keep *_softPFElectronBJetTags_*_*',
#         'keep *_softPFMuonBJetTags_*_*',
#         'keep *_softMuonTagInfos_*_*',
#         'keep *_softMuonBJetTags_*_*',
#         'keep *_softMuonByIP3dBJetTags_*_*',
#         'keep *_softMuonByPtBJetTags_*_*',
#         'keep *_combinedMVABJetTags_*_*',
#         'keep *_pfImpactParameterTagInfos_*_*',
#         'keep *_pfSecondaryVertexTagInfos_*_*',
#         'keep *_pfCombinedSecondaryVertexBJetTags_*_*'
#         )
#                                           )

#output anything produced in the ANA process
process.test_output = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "/afs/cern.ch/work/h/hardenbr/edmoutput.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW" )
    ),
                                      outputCommands = cms.untracked.vstring(
        'keep *_*_*_*'))#        'keep *_*_*_ANA'))

#run the displaced jet tags
process.load('DisplacedJets.Configuration.RecoDJTag_cff')

#config gen matching for the output displaced vertices 
if isSignalMC:
    process.displacedTagsToVertices.isSignalMC = cms.untracked.bool(True)
    process.displacedTagsToVertices.doGenMatch = cms.untracked.bool(True)

if not doedm:
    process.p = cms.Path(process.djtagging + process.analyzer)
else:
    process.p = cms.Path(process.djtagging)
    
    process.btag_output = cms.EndPath( process.test_output)
