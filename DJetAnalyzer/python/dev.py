import FWCore.ParameterSet.Config as cms

base = '/afs/cern.ch/user/h/hardenbr/2014/LL_DIJET/TRACKING_STUDIES/CMSSW_7_4_6_patch2/src/DisplacedJets/'

# output options (to be appended to the file name outputted)
#appendSignal       = "XX4J300mm"
appendSignal        = "gun30mm"
#appendSignal       = "600_30mm"
#appendSignal       = "dsusy500_10mm"
#appendSignal       = "emerge"
appendData          = ""
appendBkg           = "qcd470_600"

############ FLAGS #############
debugLevel          = 0
isSignalMC          = False
isMC                = False
isData              = True
#-------------- globaltags
#gtag               = "FALL1374_25V4"
#gtag               = "PHYS14_25_V1"
#gtag               = "MCRUN2_74_V9" #guns
gtag                = "74X_dataRun2_Prompt_v0" 
# -------------json
#JSON                = 'json_DCSONLY_Run2015B.txt'
JSON                = 'Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'
#--------------trigger
trigger_process     = "HLT"
# run related
nevents             = -1
doedm               = False
#--------------analysis todos
doEventPreSelection = False
doJetPreSelection   = False
doApplyTrigger      = True
dumpGeneralTracks   = False
# trees to write
writeTrackTree      = False
writeEventTree      = True
writeJetTree        = False
writeVertexTree     = False
writeGenTree        = False
#----------- matching
doGenMatch          = False
doSimVtxMatch       = False
#----------- analysis cuts
cut_jetPt           = 40
cut_jetEta          = 2.0
#----------- tag categories
shortTagThreshold   = 0.0
mediumTagThreshold  = 10
longTagThreshold    = 30
dHTWorkingPoint     = 2

######### input lists #########
input_file_list     = None
#input_file_list    = '/SignalMCLists/740patch1/mx300_100mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx300_300mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx300_30mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx600_30mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx100_30mm_aod.list'

# other samples
#input_file_list = '/SignalMCLists/DISPLACED_SUSY/displaced_susy.list'
#input_file_list = '/SignalMCLists/DISPLACED_SUSY/displaced_susy500_10.list'
#input_file_list = '/SignalMCLists/EMERGING_JETS/emerging_jets.txt'

# gun samples
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau30mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau0mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau300mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau10mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_flat_1mm_1000mm.list'

#data samples
#input_file_list = 'DataSampleLists/PD_DisplacedJet_Jul17AOD.txt'
input_file_list   = 'DataSampleLists/PD_JetHT_Jul17AOD.txt'
#input_file_list   = 'DataSampleLists/PD_HTMHT_Jul17AOD.txt'

# parse the input files to the file list
myfilelist = cms.untracked.vstring()
if input_file_list != None:
   list_from_input_list = open(base + input_file_list, "r") 
   lines = list_from_input_list.readlines()
   stripped_lines = map(lambda x: x.rstrip("\n"), lines)
   for line in stripped_lines:
      myfilelist.extend([line])

#fillter for the file list
if isSignalMC and input_file_list == None :
   myfilelist = cms.untracked.vstring()
   myfilelist.extend(['file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/MC_PRODUCTION/XXTo4J_M-300_CTau_30mm/AODSIM/XXTo4J_M-300_CTau-30mm_reco_102_1_ne5.root',
                      'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/MC_PRODUCTION/XXTo4J_M-300_CTau_30mm/AODSIM/XXTo4J_M-300_CTau-30mm_reco_105_1_1MO.root' ])
if not isSignalMC and input_file_list == None:
   myfilelist = cms.untracked.vstring()
   myfilelist.extend(['file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/data/Run2015B/DisplacedJet/AOD/PromptReco-v1/000/251/562/00000/F6834634-9A2A-E511-9F6F-02163E012402.root'])   
#   myfilelist.extend(['root://xrootd-cms.infn.it//store/mc/Spring14dr/QCD_Pt-470to600_Tune4C_13TeV_pythia8/AODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/000F701A-95C8-E311-89D0-002618943983.root'])   
#   myfilelist.extend(['root://xrootd-cms.infn.it//store/mc/Fall13dr/QCD_Pt-470to600_Tune4C_13TeV_pythia8/AODSIM/castor_tsg_PU20bx25_POSTLS162_V2-v1/00000/005646AA-EB79-E311-A788-00304867924E.root'])
#   myfilelist.extend(['root://xrootd-cms.infn.it//store/mc/Phys14DR/QCD_Pt-470to600_Tune4C_13TeV_pythia8/AODSIM/AVE20BX25_tsg_castor_PHYS14_25_V3-v1/00000/1E538A06-988E-E411-8F36-0025905B8562.root'])

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")

#standard sequences (from 740x driver command)
process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff') #new for navigation
#process.load('Configuration.StandardSequences.GeometryExtended_cff') #new for navigation
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff') #new for navigation
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = gtag

# process.GlobalTag = cms.ESSource( "PoolDBESSource",
#     globaltag = cms.string( gtag ),
#     RefreshEachRun = cms.untracked.bool( True ),
#     RefreshOpenIOVs = cms.untracked.bool( False ),
#     toGet = cms.VPSet(
#     ),
#     DBParameters = cms.PSet(
#       authenticationPath = cms.untracked.string( "." ),
#       connectionRetrialTimeOut = cms.untracked.int32( 60 ),
#       idleConnectionCleanupPeriod = cms.untracked.int32( 10 ),
#       messageLevel = cms.untracked.int32( 0 ),
#       enablePoolAutomaticCleanUp = cms.untracked.bool( False ),
#       enableConnectionSharing = cms.untracked.bool( True ),
#       enableReadOnlySessionOnUpdateConnection = cms.untracked.bool( False ),
#       connectionTimeOut = cms.untracked.int32( 0 ),
#       connectionRetrialPeriod = cms.untracked.int32( 10 )
#     ),
#     RefreshAlways = cms.untracked.bool( False ),
#     connect = cms.string( "frontier://(proxyurl=http://localhost:3128)(serverurl=http://localhost:8000/FrontierOnProd)(serverurl=http://localhost:8000/FrontierOnProd)(retrieve-ziplevel=0)/CMS_CONDITIONS" ),
#     ReconnectEachRun = cms.untracked.bool( True ),
#     BlobStreamerName = cms.untracked.string( "TBufferBlobStreamingService" ),
#     DumpStat = cms.untracked.bool( False )
# )

# process.GlobalTag = cms.ESSource( "PoolDBESSource",
#     globaltag = cms.string( gtag ), 
#     RefreshEachRun = cms.untracked.bool( True ),
#     RefreshOpenIOVs = cms.untracked.bool( False ),
#     toGet = cms.VPSet( 
#       cms.PSet(  record = cms.string( "JetCorrectionsRecord" ),
#         tag = cms.string( "JetCorrectorParametersCollection_HLT_V1_AK4Calo" ),
#         connect = cms.untracked.string( "frontier://FrontierPrep/CMS_COND_PHYSICSTOOLS" ),
#         label = cms.untracked.string( "AK4CaloHLT" )
#       ),
#       cms.PSet(  record = cms.string( "JetCorrectionsRecord" ),
#         tag = cms.string( "JetCorrectorParametersCollection_HLT_trk1B_V1_AK4PF" ),
#         connect = cms.untracked.string( "frontier://FrontierPrep/CMS_COND_PHYSICSTOOLS" ),
#         label = cms.untracked.string( "AK4PFHLT" )
#       )
#     ),
#     DBParameters = cms.PSet( 
#       authenticationPath = cms.untracked.string( "." ),
#       connectionRetrialTimeOut = cms.untracked.int32( 60 ),
#       idleConnectionCleanupPeriod = cms.untracked.int32( 10 ),
#       messageLevel = cms.untracked.int32( 0 ),
#       enablePoolAutomaticCleanUp = cms.untracked.bool( False ),
#       enableConnectionSharing = cms.untracked.bool( True ),
#       enableReadOnlySessionOnUpdateConnection = cms.untracked.bool( False ),
#       connectionTimeOut = cms.untracked.int32( 0 ),
#       connectionRetrialPeriod = cms.untracked.int32( 10 )
#     ),
#     RefreshAlways = cms.untracked.bool( False ),
#     connect = cms.string( "frontier://(proxyurl=http://localhost:3128)(serverurl=http://localhost:8000/FrontierOnProd)(serverurl=http://localhost:8000/FrontierOnProd)(retrieve-ziplevel=0)/CMS_CONDITIONS" ),
#     ReconnectEachRun = cms.untracked.bool( True ),
#     BlobStreamerName = cms.untracked.string( "TBufferBlobStreamingService" )
# )

# Deliver the missing payloads fro the globaltag (TO BE REMOVED)
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
# add a JSON for data
import FWCore.PythonUtilities.LumiList as LumiList
if isData:
   process.source.lumisToProcess = LumiList.LumiList(filename = JSON).getVLuminosityBlockRange()

################################################################################################

# add the jet corrections
process.hltAK4CaloFastJetCorrector = cms.EDProducer( "L1FastjetCorrectorProducer",
    srcRho = cms.InputTag( "fixedGridRhoFastjetAllCalo" ),
    algorithm = cms.string( "AK4CaloHLT" ),
    level = cms.string( "L1FastJet" )
)
process.hltAK4CaloRelativeCorrector = cms.EDProducer( "LXXXCorrectorProducer",
    algorithm = cms.string( "AK4CaloHLT" ),
    level = cms.string( "L2Relative" )
)
process.hltAK4CaloAbsoluteCorrector = cms.EDProducer( "LXXXCorrectorProducer",
    algorithm = cms.string( "AK4CaloHLT" ),
    level = cms.string( "L3Absolute" )
)
process.hltAK4CaloCorrector = cms.EDProducer( "ChainedJetCorrectorProducer",
    correctors = cms.VInputTag( 'hltAK4CaloFastJetCorrector','hltAK4CaloRelativeCorrector','hltAK4CaloAbsoluteCorrector' )
)
process.ak4CaloJetsCorrected = cms.EDProducer( "CorrectedCaloJetProducer",
    src = cms.InputTag( "ak4CaloJets" ),
    correctors = cms.VInputTag( 'hltAK4CaloCorrector' )
)

process.correctJets = cms.Sequence( process.hltAK4CaloFastJetCorrector + process.hltAK4CaloRelativeCorrector + process.hltAK4CaloAbsoluteCorrector + process.hltAK4CaloCorrector + process.ak4CaloJetsCorrected)

#configure the analyzers
process.analyzerVTX = cms.EDAnalyzer('DJetAnalyzer')
process.analyzerCALO = cms.EDAnalyzer('DJetAnalyzer')

process.analyzerVTX.debugLevel  = cms.untracked.int32(debugLevel)
process.analyzerCALO.debugLevel = cms.untracked.int32(debugLevel)

#output configuration
if isSignalMC:
    process.analyzerVTX.outputFileName = cms.untracked.string('signalVTX%s.root' % appendSignal)
    process.analyzerCALO.outputFileName = cms.untracked.string('signalCALO%s.root' % appendSignal)
elif not isSignalMC and isMC:
#    process.analyzerVTX.outputFileName = cms.untracked.string('qcdVTX%s.root' % appendBkg)
#    process.analyzerCALO.outputFileName = cms.untracked.string('qcdCALO%s.root' % appendBkg)
   process.analyzerCALO.outputFileName = cms.untracked.string(outputfile_string)
else:
   process.analyzerVTX.outputFileName  = cms.untracked.string('dataVTX%s.root' % appendData)
   process.analyzerCALO.outputFileName = cms.untracked.string('dataCALO%s.root' % appendData)
   

#tree names
process.analyzerCALO.jetTreeName    = cms.untracked.string('jets')
process.analyzerCALO.trackTreeName  = cms.untracked.string('tracks')
process.analyzerCALO.vertexTreeName = cms.untracked.string('vtx')
process.analyzerCALO.genTreeName    = cms.untracked.string('genp')

process.analyzerVTX.jetTreeName     = cms.untracked.string('jets')
process.analyzerVTX.trackTreeName   = cms.untracked.string('tracks')
process.analyzerVTX.vertexTreeName  = cms.untracked.string('vtx')
process.analyzerVTX.genTreeName     = cms.untracked.string('genp')

# analysis dependent flags
#  applyEventPreSelection_ = iConfig.getUntrackedParameter<bool>("applyEventPreSelection");
#  applyJetPreSelection_   = iConfig.getUntrackedParameter<bool>("applyJetPreSelection");
process.analyzerCALO.applyEventPreSelection = cms.untracked.bool(doEventPreSelection)
process.analyzerCALO.applyJetPreSelection   = cms.untracked.bool(doJetPreSelection)
process.analyzerCALO.dumpGeneralTracks      = cms.untracked.bool(dumpGeneralTracks)

# what to write out 
process.analyzerCALO.writeTrackTree  = cms.untracked.bool(writeTrackTree)
process.analyzerCALO.writeEventTree  = cms.untracked.bool(writeEventTree)
process.analyzerCALO.writeJetTree    = cms.untracked.bool(writeJetTree)
process.analyzerCALO.writeVertexTree = cms.untracked.bool(writeVertexTree)
process.analyzerCALO.writeGenTree    = cms.untracked.bool(writeGenTree)

# MC dependent flags
process.analyzerVTX.isMC        = cms.untracked.bool(isMC)
process.analyzerCALO.isMC       = cms.untracked.bool(isMC)
process.analyzerVTX.isSignalMC  = cms.untracked.bool(isSignalMC)
process.analyzerCALO.isSignalMC = cms.untracked.bool(isSignalMC)
process.analyzerVTX.doGenMatch  = cms.untracked.bool(doGenMatch)
process.analyzerCALO.doGenMatch = cms.untracked.bool(doGenMatch)
process.analyzerVTX.doSimMatch  = cms.untracked.bool(doSimVtxMatch)
process.analyzerCALO.doSimMatch = cms.untracked.bool(doSimVtxMatch)

#trigger
process.analyzerVTX.triggerResultPath  = cms.untracked.string(trigger_process)
process.analyzerCALO.triggerResultPath = cms.untracked.string(trigger_process)
process.analyzerVTX.triggerResults     = cms.untracked.InputTag('TriggerResults', '', 'HLT')
process.analyzerCALO.triggerResults    = cms.untracked.InputTag('TriggerResults', '', 'HLT')

# collection tags
process.analyzerVTX.generalTracks      = cms.untracked.InputTag('generalTracks', '', '')
process.analyzerVTX.ak4CaloJets        = cms.untracked.InputTag('ak4CaloJets', '', '')
process.analyzerVTX.genParticles       = cms.untracked.InputTag('genParticles', '', '')
process.analyzerCALO.generalTracks     = cms.untracked.InputTag('generalTracks', '', '')
process.analyzerCALO.ak4CaloJets       = cms.untracked.InputTag('ak4CaloJetsCorrected', '', '')
process.analyzerCALO.genParticles      = cms.untracked.InputTag('genParticles', '', '')

# jet tagging categories
process.analyzerCALO.shortTagThreshold  = cms.untracked.double(shortTagThreshold)
process.analyzerCALO.mediumTagThreshold = cms.untracked.double(mediumTagThreshold)
process.analyzerCALO.longTagThreshold   = cms.untracked.double(longTagThreshold)
process.analyzerCALO.dHTWorkingPoint   = cms.untracked.int32(dHTWorkingPoint)

# vertex matched ip info
process.analyzerVTX.secondaryVertexTagInfo   = cms.untracked.InputTag('displacedSecondaryVertexTagInfosNoPV', '', 'ANA')
process.analyzerVTX.lifetimeIPTagInfo        = cms.untracked.InputTag('displacedLifetimeTagInfos', '', 'ANA')

# calo matched ip and sv info
process.analyzerCALO.secondaryVertexTagInfo   = cms.untracked.InputTag('displacedSecondaryVertexTagInfosNoPV', '', 'ANA')
process.analyzerCALO.lifetimeIPTagInfo        = cms.untracked.InputTag('displacedLifetimeTagInfos', '', 'ANA') #change to caloface later

##reconstructed vertex information
process.analyzerCALO.secondaryVertex          = cms.untracked.InputTag('displacedTagsToVerticesCaloFace', '', 'ANA')
process.analyzerCALO.inclusiveVertexCand      = cms.untracked.InputTag('displacedInclusiveVertexFinder', '', 'ANA')
process.analyzerCALO.inclusiveVertexSecondary = cms.untracked.InputTag('displacedInclusiveSecondaryVertices', '', 'ANA')

process.analyzerVTX.secondaryVertex          = cms.untracked.InputTag('displacedTagsToVertices', '', 'ANA')
process.analyzerVTX.inclusiveVertexCand      = cms.untracked.InputTag('displacedInclusiveVertexFinder', '', 'ANA')
process.analyzerVTX.inclusiveVertexSecondary = cms.untracked.InputTag('displacedInclusiveSecondaryVertices', '', 'ANA')

# primary vertex
process.analyzerCALO.offlinePrimaryVertices = cms.untracked.InputTag('offlinePrimaryVerticesWithBS')
process.analyzerVTX.offlinePrimaryVertices  = cms.untracked.InputTag('offlinePrimaryVerticesWithBS')

#cuts
process.analyzerVTX.jetPt   = cms.untracked.double(cut_jetPt)
process.analyzerVTX.jetEta  = cms.untracked.double(cut_jetEta)
process.analyzerCALO.jetPt  = cms.untracked.double(cut_jetPt)
process.analyzerCALO.jetEta = cms.untracked.double(cut_jetEta)

# Tags related to the monte carlo
if isMC:
    # Gen Information
    process.analyzerVTX.ak5GenJets    = cms.untracked.InputTag('ak5GenJets', '', '')
    process.analyzerVTX.genMetCalo    = cms.untracked.InputTag('genMetCalo', '', '')
    process.analyzerVTX.genParticles  = cms.untracked.InputTag('genParticles', '', '')
    
    process.analyzerCALO.ak5GenJets   = cms.untracked.InputTag('ak5GenJets', '', '')
    process.analyzerCALO.genMetCalo   = cms.untracked.InputTag('genMetCalo', '', '')
    process.analyzerCALO.genParticles = cms.untracked.InputTag('genParticles', '', '')

    process.analyzerCALO.simVertices = cms.untracked.InputTag('g4SimHits', '', '')
    process.analyzerVTX.simVertices  = cms.untracked.InputTag('g4SimHits', '', '')        

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
#process.load('DisplacedJets.Configuration.AdditionalPATSequences_cff')
#process.load('DisplacedJets/DisplacedTriggerFilters/displacedTriggers_cff')

#create the main path to run
process.p = cms.Path()

# trigger bits to keep
process.triggerSelection = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring(
      'HLT_HT350_DisplacedDijet40_DisplacedTrack_v*',
      'HLT_HT350_DisplacedDijet80_DisplacedTrack_v*',
      'HLT_HT500_DisplacedDijet40_Inclusive_v*',
      'HLT_HT550_DisplacedDijet40_Inclusive_v*',
      'HLT_HT650_DisplacedDijet80_Inclusive_v*',
      'HLT_HT750_DisplacedDijet80_Inclusive_v*',
      'HLT_VBF_DisplacedJet40_Hadronic_v*',
      'HLT_VBF_DisplacedJet40_DisplacedTrack_v*',
      'HLT_VBF_DisplacedJet40_DisplacedTrack_v*',                                                  
      'HLT_VBF_DisplacedJet40_Hadronic_v*',                                                        
      'HLT_VBF_DisplacedJet40_TightID_DisplacedTrack_v*',                                          
      'HLT_VBF_DisplacedJet40_TightID_Hadronic_v*',                                                
      'HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack_v*',                                         
      'HLT_VBF_DisplacedJet40_VTightID_Hadronic_v*',                                               
      'HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack_v*',                                        
      'HLT_VBF_DisplacedJet40_VVTightID_Hadronic_v*',   
      'HLT_PFHT200_v*',
      'HLT_PFHT250_v*',
      'HLT_PFHT300_v*',
      'HLT_PFHT350_v*',
      'HLT_PFHT400_v*',
      'HLT_PFHT800_v*',
      'HLT_L1_TripleJet_VBF_v*'                                                                  
      ),
    hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
    l1tResults = cms.InputTag( "" ),
    l1tIgnoreMask = cms.bool( False ),
    l1techIgnorePrescales = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True )
)

# from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection

# # add the pat jets                                                                                                                                                          
# addJetCollection(
#    process,
#    labelName = 'AK4Calo',
#    jetSource = cms.InputTag('ak4CaloJets'),
#    jetCorrections = ('AK4Calo', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-1'), 
#    btagDiscriminators = ['None'] 
#    )

# #add the pat related sequences
# process.p *= process.makePatMETs

# process.MessageLogger = cms.Service( "MessageLogger",
#     suppressInfo = cms.untracked.vstring(  ),
#     debugs = cms.untracked.PSet( 
#       threshold = cms.untracked.string( "INFO" ),
#       placeholder = cms.untracked.bool( True ),
#       suppressInfo = cms.untracked.vstring(  ),
#       suppressWarning = cms.untracked.vstring(  ),
#       suppressDebug = cms.untracked.vstring(  ),
#       suppressError = cms.untracked.vstring(  )
#     ),
#     suppressDebug = cms.untracked.vstring(  ),
#     cout = cms.untracked.PSet(  placeholder = cms.untracked.bool( True ) ),
#     cerr_stats = cms.untracked.PSet( 
#       threshold = cms.untracked.string( "WARNING" ),
#       output = cms.untracked.string( "cerr" ),
#       optionalPSet = cms.untracked.bool( True )
#     ),
#     warnings = cms.untracked.PSet( 
#       threshold = cms.untracked.string( "INFO" ),
#       placeholder = cms.untracked.bool( True ),
#       suppressInfo = cms.untracked.vstring(  ),
#       suppressWarning = cms.untracked.vstring(  ),
#       suppressDebug = cms.untracked.vstring(  ),
#       suppressError = cms.untracked.vstring(  )
#     ),
#     statistics = cms.untracked.vstring( 'cerr' ),
#     cerr = cms.untracked.PSet( 
#       INFO = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
#       noTimeStamps = cms.untracked.bool( False ),
#       FwkReport = cms.untracked.PSet( 
#         reportEvery = cms.untracked.int32( 1 ),
#         limit = cms.untracked.int32( 0 )
#       ),
#       default = cms.untracked.PSet(  limit = cms.untracked.int32( 10000000 ) ),
#       Root_NoDictionary = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
#       FwkJob = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
#       FwkSummary = cms.untracked.PSet( 
#         reportEvery = cms.untracked.int32( 1 ),
#         limit = cms.untracked.int32( 10000000 )
#       ),
#       threshold = cms.untracked.string( "INFO" ),
#       suppressInfo = cms.untracked.vstring(  ),
#       suppressWarning = cms.untracked.vstring(  ),
#       suppressDebug = cms.untracked.vstring(  ),
#       suppressError = cms.untracked.vstring(  )
#     ),
#     FrameworkJobReport = cms.untracked.PSet( 
#       default = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
#       FwkJob = cms.untracked.PSet(  limit = cms.untracked.int32( 10000000 ) )
#     ),
#     suppressWarning = cms.untracked.vstring( 'hltOnlineBeamSpot',
#       'hltCtf3HitL1SeededWithMaterialTracks',
#       'hltL3MuonsOIState',
#       'hltPixelTracksForHighMult',
#       'hltHITPixelTracksHE',
#       'hltHITPixelTracksHB',
#       'hltCtfL1SeededWithMaterialTracks',
#       'hltRegionalTracksForL3MuonIsolation',
#       'hltSiPixelClusters',
#       'hltActivityStartUpElectronPixelSeeds',
#       'hltLightPFTracks',
#       'hltPixelVertices3DbbPhi',
#       'hltL3MuonsIOHit',
#       'hltPixelTracks',
#       'hltSiPixelDigis',
#       'hltL3MuonsOIHit',
#       'hltL1SeededElectronGsfTracks',
#       'hltL1SeededStartUpElectronPixelSeeds',
#       'hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJetFastPV',
#       'hltCtfActivityWithMaterialTracks' ),
#     errors = cms.untracked.PSet( 
#       threshold = cms.untracked.string( "INFO" ),
#       placeholder = cms.untracked.bool( True ),
#       suppressInfo = cms.untracked.vstring(  ),
#       suppressWarning = cms.untracked.vstring(  ),
#       suppressDebug = cms.untracked.vstring(  ),
#       suppressError = cms.untracked.vstring(  )
#     ),
#     fwkJobReports = cms.untracked.vstring( 'FrameworkJobReport' ),
#     debugModules = cms.untracked.vstring(  ),
#     infos = cms.untracked.PSet( 
#       threshold = cms.untracked.string( "INFO" ),
#       Root_NoDictionary = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
#       placeholder = cms.untracked.bool( True ),
#       suppressInfo = cms.untracked.vstring(  ),
#       suppressWarning = cms.untracked.vstring(  ),
#       suppressDebug = cms.untracked.vstring(  ),
#       suppressError = cms.untracked.vstring(  )
#     ),
#     categories = cms.untracked.vstring( 'FwkJob',
#       'FwkReport',
#       'FwkSummary',
#       'Root_NoDictionary' ),
#     destinations = cms.untracked.vstring( 'warnings',
#       'errors',
#       'infos',
#       'debugs',
#       'cout',
#       'cerr' ),
#     threshold = cms.untracked.string( "INFO" ),
#     suppressError = cms.untracked.vstring( 'hltOnlineBeamSpot',
#       'hltL3MuonCandidates',
#       'hltL3TkTracksFromL2OIState',
#       'hltPFJetCtfWithMaterialTracks',
#       'hltL3TkTracksFromL2IOHit',
#       'hltL3TkTracksFromL2OIHit' )
# )


if doApplyTrigger: #apply the triggers and run dj tagging
   process.p *= process.triggerSelection *  process.correctJets * process.djtagging
else: #just run the tagging sequence
   process.p *= process.correctJets * process.djtagging

if doedm: #just dump the edm output of the djtagging sequence
    process.btag_output = cms.EndPath(process.test_output)
else: # run the analysis and tree dumper 
    process.p *= process.analyzerCALO

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool( True )
)
