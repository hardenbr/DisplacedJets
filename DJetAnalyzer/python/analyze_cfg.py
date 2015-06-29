import FWCore.ParameterSet.Config as cms


############ FLAGS #############
#globaltags
gtag = "MCRUN2_74_V9"
# run related
nevents             = 100
debugLevel          = 0
doedm               = False
# sample related
isSignalMC          = True
isMC                = True
# analysis related
doEventPreSelection = False
doJetPreSelection   = False
doApplyTrigger      = False
# matching
doGenMatch          = True
doSimVtxMatch       = False
# analysis cuts
cut_jetPt           = 40
cut_jetEta          = 2.0
# tag cateogires
shortTagThreshold    = 0.0
mediumTagThreshold  = 10
longTagThreshold    = 30
dHTWorkingPoint     = 2

#fillter for the file list
if isSignalMC and input_file_list == None :
   myfilelist = cms.untracked.vstring()
   myfilelist.extend(['file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/MC_PRODUCTION/XXTo4J_M-300_CTau_30mm/AODSIM/XXTo4J_M-300_CTau-30mm_reco_102_1_ne5.root',
                      'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/MC_PRODUCTION/XXTo4J_M-300_CTau_30mm/AODSIM/XXTo4J_M-300_CTau-30mm_reco_105_1_1MO.root' ])
if not isSignalMC:
   myfilelist.extend(['root://xrootd-cms.infn.it//store/mc/Spring14dr/QCD_Pt-470to600_Tune4C_13TeV_pythia8/AODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/000F701A-95C8-E311-89D0-002618943983.root'])   

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")

#standard sequences (from 740x driver command)
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag = cms.ESSource( "PoolDBESSource",
    globaltag = cms.string( gtag ), 
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

# # Deliver the missing payloads fro the globaltag (TO BE REMOVED)
# if 'GlobalTag' in process.__dict__:
#     from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
#     process.GlobalTag           = customiseGlobalTag(process.GlobalTag, globaltag = gtag)
#     process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_CONDITIONS'
#     process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
#     for pset in process.GlobalTag.toGet.value():
#         pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
#     # fix for multi-run processing
#     process.GlobalTag.RefreshEachRun = cms.untracked.bool( False )
#     process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nevents))

process.source = cms.Source("PoolSource", fileNames = myfilelist )

################################################################################################

#configure the analyzers
process.analyzerCALO = cms.EDAnalyzer('DJetAnalyzer')
process.analyzerCALO.debugLevel = cms.untracked.int32(debugLevel)

#output configuration
if isSignalMC:
    process.analyzerCALO.outputFileName = cms.untracked.string('signalCALO%s.root' % appendLifetime)
else:
#    process.analyzerCALO.outputFileName = cms.untracked.string('qcdCALO%s.root' % appendBkg)
   process.analyzerCALO.outputFileName = cms.untracked.string(outputfile_string)

#tree names
process.analyzerCALO.jetTreeName    = cms.untracked.string('jets')
process.analyzerCALO.trackTreeName  = cms.untracked.string('tracks')
process.analyzerCALO.vertexTreeName = cms.untracked.string('vtx')
process.analyzerCALO.genTreeName    = cms.untracked.string('genp')

# analysis dependent flags
process.analyzerCALO.applyEventPreSelection = cms.untracked.bool(doEventPreSelection)
process.analyzerCALO.applyJetPreSelection   = cms.untracked.bool(doJetPreSelection)

# MC dependent flags
process.analyzerCALO.isMC       = cms.untracked.bool(isMC)
process.analyzerCALO.isSignalMC = cms.untracked.bool(isSignalMC)
process.analyzerCALO.doGenMatch = cms.untracked.bool(doGenMatch)
process.analyzerCALO.doSimMatch = cms.untracked.bool(doSimVtxMatch)

#tags
process.analyzerCALO.generalTracks = cms.untracked.InputTag('generalTracks', '', '')
process.analyzerCALO.ak4CaloJets   = cms.untracked.InputTag('ak4CaloJets', '', '')
process.analyzerCALO.genParticles  = cms.untracked.InputTag('genParticles', '', '')

# jet tagging categories
process.analyzerCALO.shortTagThreshold  = cms.untracked.double(shortTagThreshold)
process.analyzerCALO.mediumTagThreshold = cms.untracked.double(mediumTagThreshold)
process.analyzerCALO.longTagThreshold   = cms.untracked.double(longTagThreshold)
process.analyzerCALO.dHTWorkingPoint   = cms.untracked.int32(dHTWorkingPoint)

# calo matched ip and sv info
process.analyzerCALO.secondaryVertexTagInfo   = cms.untracked.InputTag('displacedSecondaryVertexTagInfosNoPVCaloFace', '', 'ANA')
process.analyzerCALO.lifetimeIPTagInfo        = cms.untracked.InputTag('displacedLifetimeTagInfosCaloFace', '', 'ANA')

##reconstructed vertex information
process.analyzerCALO.secondaryVertex          = cms.untracked.InputTag('displacedTagsToVerticesCaloFace', '', 'ANA')
process.analyzerCALO.inclusiveVertexCand      = cms.untracked.InputTag('displacedInclusiveVertexFinder', '', 'ANA')
process.analyzerCALO.inclusiveVertexSecondary = cms.untracked.InputTag('displacedInclusiveSecondaryVertices', '', 'ANA')

# primary vertex
process.analyzerCALO.offlinePrimaryVertices = cms.untracked.InputTag('offlinePrimaryVerticesWithBS')

#cuts
process.analyzerCALO.jetPt  = cms.untracked.double(cut_jetPt)
process.analyzerCALO.jetEta = cms.untracked.double(cut_jetEta)

# Tags related to the monte carlo
if isMC:
    # Gen Information    
    process.analyzerCALO.ak5GenJets   = cms.untracked.InputTag('ak5GenJets', '', '')
    process.analyzerCALO.genMetCalo   = cms.untracked.InputTag('genMetCalo', '', '')
    process.analyzerCALO.genParticles = cms.untracked.InputTag('genParticles', '', '')
    process.analyzerCALO.simVertices = cms.untracked.InputTag('g4SimHits', '', '')

#output anything produced in the ANA process
process.test_output = cms.OutputModule( "PoolOutputModule",
                                        fileName = cms.untracked.string( "/afs/cern.ch/work/h/hardenbr/edmoutput.root" ),
                                        fastCloning = cms.untracked.bool( False ),
                                        dataset = cms.untracked.PSet(filterName = cms.untracked.string( "" ),
                                                                      dataTier = cms.untracked.string( "RAW" )
                                                                      ),
                                        outputCommands = cms.untracked.vstring('keep *_*_*_*'))

# run the displaced jet tags
process.load('DisplacedJets.Configuration.RecoDJTag_cff')
#process.load('DisplacedJets/DisplacedTriggerFilters/displacedTriggers_cff')

#create the main path to run
process.p = cms.Path()

if doApplyTrigger: #apply the triggers and run dj tagging
   process.p *= process.InclusiveTrigger * process.DisplacedTracktrigger * process.djtagging
else: #just run the tagging sequence
   process.p *= process.djtagging

if doedm: #just dump the edm output of the djtagging sequence
    process.btag_output = cms.EndPath(process.test_output)
else: # run the analysis and tree dumper 
    process.p *= process.analyzerCALO


# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs                                                                  
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs                                                           
process = customisePostLS1(process)
