import FWCore.ParameterSet.Config as cms

#output directories
base                = '/afs/cern.ch/user/h/hardenbr/2014/LL_DIJET/TRACKING_STUDIES/CMSSW_7_4_10_patch1/src/DisplacedJets/'
#outputDir           = "/afs/cern.ch/work/h/hardenbr/2015/DIJET/DJANALYSIS/DATA_SEPT21/"
outputDir           = ""

# output options (to be appended to the file name outputted)
appendSignal        = ""
appendData          = ""
appendBkg           = ""
############ FLAGS #############
debugLevel          = 0
reportEveryNEvents  = 1000
isSignalMC          = False
isMC                = False
isData              = not isMC
doedm               = False
nevents             = -1

#-------------- globaltags
#gtag               = "74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0" #spring 15 25ns
#gtag               = "FALL1374_25V4"
#gtag               = "PHYS14_25_V1"
#gtag                = "MCRUN2_74_V9" #guns
gtag                = "74X_dataRun2_Prompt_v2"  #data
# -------------json
#JSON               = 'json_DCSONLY_Run2015B.txt'
#JSON                = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt'
JSON                = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'

#--------------trigger
trigger_process     = "" if isMC else "HLT"

#--------------analysis todos
doEventPreSelection = False
doJetPreSelection   = False
doApplySingleMu     = False
doApplyTrigger      = True  #if isData else False
dumpGeneralTracks   = False
# trees to write
writeTrackTree      = False
writeV0Tree         = True
writeEventTree      = True
writeJetTree        = True
writeVertexTree     = True
writeGenTree        = isSignalMC 
#----------- matching
doGenMatch          = True and isSignalMC
doSimVtxMatch       = False
#----------- analysis cuts
cut_jetPt           = 40
cut_jetEta          = 2.0
#----------- tag categories
shortTagThreshold   = 0.0
mediumTagThreshold  = 0.0
longTagThreshold    = 30
dHTWorkingPoint     = 2

######### input lists #########

input_file_list     = None
#----xx4j samples
#input_file_list    = '/SignalMCLists/740patch1/mx300_1000mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx300_300mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx300_100mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx300_30mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx300_1mm_aod.list'
#input_file_list     = '/SignalMCLists/740patch1/mx300_0p1mm_aod.list'

#---- xx4j mass variants 
#input_file_list    = '/SignalMCLists/740patch1/mx600_30mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx100_30mm_aod.list'

#----other samples
#input_file_list = '/SignalMCLists/DISPLACED_SUSY/displaced_susy.list'
#input_file_list = '/SignalMCLists/DISPLACED_SUSY/displaced_susy500_10.list'
#input_file_list = '/SignalMCLists/EMERGING_JETS/emerging_jets.txt'

#-----gun samples
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau30mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau0mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau300mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau10mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_flat_1mm_1000mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau0mm_bbar.list'

#----data samples
input_file_list  = 'DataSampleLists/PD_DisplacedJet_Run2015D_Sept21.txt'
#input_file_list = 'DataSampleLists/PD_DisplacedJet_Jul17AOD.txt'
#input_file_list = 'DataSampleLists/PD_JetHT_Jul17AOD.txt'
#input_file_list = 'DataSampleLists/PD_HTMHT_Jul17AOD.txt'
#input_file_list = 'DataSampleLists/PD_SingleMuon_Jul23AOD.txt'

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
   print "NO SIGNAL INPUT.....Exiting" 
   exit(1)
#   myfilelist.extend(['file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/MC_PRODUCTION/XXTo4J_M-300_CTau_30mm/AODSIM/XXTo4J_M-300_CTau-30mm_reco_102_1_ne5.root',
 #                     'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/MC_PRODUCTION/XXTo4J_M-300_CTau_30mm/AODSIM/XXTo4J_M-300_CTau-30mm_reco_105_1_1MO.root' ])
if not isSignalMC and input_file_list == None and not isData:
   myfilelist = cms.untracked.vstring()
#   myfilelist.extend(['file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/data/Run2015B/DisplacedJet/AOD/PromptReco-v1/000/251/562/00000/F6834634-9A2A-E511-9F6F-02163E012402.root'])   
   myfilelist.extend(['/store/mc/RunIISpring15DR74/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/00D76158-CCFC-E411-89EA-AC853DA06B56.root'])
if isData and input_file_list == None:
   myfilelist.extend(['/store/data/Run2015C/JetHT/AOD/PromptReco-v1/000/254/905/00000/263140AF-B34B-E511-A678-02163E0146DB.root'])
#   myfilelist.extend(['file:pickevents.root'])

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
#standard sequences (from 740x driver command)
process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = reportEveryNEvents
process.MessageLogger.warnings.suppressInfo = cms.untracked.vstring()
supressWarnings = ["TwoTrackMinimumDistance", "displacedInclusiveVertexFinder","inclusiveVertexFinder","InclusiveVertexFinder"]
for warning in supressWarnings: process.MessageLogger.suppressWarning.extend([warning])

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff') #new for navigation
#process.load('Configuration.StandardSequences.GeometryExtended_cff') #new for navigation
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff') #new for navigation
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
#global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.load('JetMETCorrections.Configuration.JetCorrectorsForReco_cff')
#process.load('JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('JetMETCorrections.Configuration.CorrectedJetProducers_cff')

process.GlobalTag.globaltag = gtag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nevents))
process.source = cms.Source("PoolSource", fileNames = myfilelist )

# add a JSON for data
import FWCore.PythonUtilities.LumiList as LumiList
if isData:
 process.source.lumisToProcess = LumiList.LumiList(filename = JSON).getVLuminosityBlockRange()

################################################################################################

#jet corrections taken from the global tag
process.correctJets = cms.Sequence( process.ak4CaloL2L3CorrectorChain * process.ak4CaloJetsL2L3)

#configure the analyzers
process.analyzerVTX = cms.EDAnalyzer('DJetAnalyzer')
process.analyzerCALO = cms.EDAnalyzer('DJetAnalyzer')

process.analyzerVTX.debugLevel  = cms.untracked.int32(debugLevel)
process.analyzerCALO.debugLevel = cms.untracked.int32(debugLevel)

#output configuration
if isSignalMC:
    process.analyzerVTX.outputFileName = cms.untracked.string('signalVTX%s.root' % appendSignal)
    process.analyzerCALO.outputFileName = cms.untracked.string('%ssignal%s.root' % (outputDir, appendSignal))
elif not isSignalMC and isMC:
#    process.analyzerVTX.outputFileName = cms.untracked.string('qcdVTX%s.root' % appendBkg)
#    process.analyzerCALO.outputFileName = cms.untracked.string('qcdCALO%s.root' % appendBkg)
   process.analyzerCALO.outputFileName = cms.untracked.string("%sdjana.root" % outputDir)
else:
   process.analyzerVTX.outputFileName  = cms.untracked.string('dataVTX%s.root' % appendData)
   process.analyzerCALO.outputFileName = cms.untracked.string('%data%s.root' % (outputDir, appendData))   


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
process.analyzerCALO.writeV0Tree    = cms.untracked.bool(writeV0Tree)
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
process.analyzerVTX.triggerResults     = cms.untracked.InputTag('TriggerResults', '', '')
process.analyzerCALO.triggerResults    = cms.untracked.InputTag('TriggerResults', '', trigger_process)

# collection tags
process.analyzerVTX.generalTracks                  = cms.untracked.InputTag('generalTracks', '', '')
process.analyzerVTX.ak4CaloJets                    = cms.untracked.InputTag('ak4CaloJets', '', '')
process.analyzerVTX.genParticles                   = cms.untracked.InputTag('genParticles', '', '')
process.analyzerCALO.generalTracks                 = cms.untracked.InputTag('generalTracks', '', '')
process.analyzerCALO.ak4CaloJets                   = cms.untracked.InputTag('ak4CaloJetsL2L3', '', '')
process.analyzerCALO.genParticles                  = cms.untracked.InputTag('genParticles', '', '')
process.analyzerCALO.caloMatchedTrackAssociation   = cms.untracked.InputTag('displacedAk4JetTracksAssociatorAtCaloFace','','ANA')
process.analyzerCALO.vertexMatchedTrackAssociation = cms.untracked.InputTag('displacedAk4JetTracksAssociatorAtVertex','','ANA')


# jet tagging categories
process.analyzerCALO.shortTagThreshold  = cms.untracked.double(shortTagThreshold)
process.analyzerCALO.mediumTagThreshold = cms.untracked.double(mediumTagThreshold)
process.analyzerCALO.longTagThreshold   = cms.untracked.double(longTagThreshold)
process.analyzerCALO.dHTWorkingPoint    = cms.untracked.int32(dHTWorkingPoint)

# jet alpha 
#process.analyzerCALO.jetVertexAssociation = cms.untracked.InputTag("displacedJetVertexAssociation","Var","ANA")

# vertex matched ip info
process.analyzerVTX.secondaryVertexTagInfo    = cms.untracked.InputTag('displacedSecondaryVertexTagInfosNoPV', '', 'ANA')
process.analyzerVTX.lifetimeIPTagInfo         = cms.untracked.InputTag('displacedLifetimeTagInfos', '', 'ANA')
# calo matched ip and sv info
process.analyzerCALO.secondaryVertexTagInfo   = cms.untracked.InputTag('displacedSecondaryVertexTagInfosNoPV', '', 'ANA')

# use the tracks matched at the inner hit for the ip tags
process.analyzerCALO.lifetimeIPTagInfo        = cms.untracked.InputTag('displacedLifetimeTagInfos', '', 'ANA') #change to final associaton later
##reconstructed vertex information
process.analyzerCALO.secondaryVertex          = cms.untracked.InputTag('displacedTagsToVertices', '', 'ANA') 
process.analyzerCALO.inclusiveVertexCand      = cms.untracked.InputTag('displacedInclusiveVertexFinder', '', 'ANA')
process.analyzerCALO.inclusiveVertexSecondary = cms.untracked.InputTag('displacedInclusiveSecondaryVertices', '', 'ANA')
process.analyzerVTX.secondaryVertex           = cms.untracked.InputTag('displacedTagsToVertices', '', 'ANA')
process.analyzerVTX.inclusiveVertexCand       = cms.untracked.InputTag('displacedInclusiveVertexFinder', '', 'ANA')
process.analyzerVTX.inclusiveVertexSecondary  = cms.untracked.InputTag('displacedInclusiveSecondaryVertices', '', 'ANA')

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
#process.displacedLifetimeTagInfos.jetTracks   = cms.InputTag( "displacedAk4JetTracksAssociatorAtInnerHit" )
process.displacedLifetimeTagInfos.jetTracks   = cms.InputTag( "displacedAk4JetTracksAssociatorAtVertex" )
#process.load('DisplacedJets.Configuration.AdditionalPATSequences_cff')
#process.load('DisplacedJets/DisplacedTriggerFilters/displacedTriggers_cff')

#create the main path to run
process.p = cms.Path()

# trigger bits to keep
process.triggerSelection = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring(
#      'HLT_HT250_DisplacedDijet40_DisplacedTrack_v*',
      'HLT_HT350_DisplacedDijet40_DisplacedTrack_v*',
      'HLT_HT350_DisplacedDijet80_DisplacedTrack_v*',
#      'HLT_HT400_DisplacedDijet40_Inclusive_v*',
      'HLT_HT500_DisplacedDijet40_Inclusive_v*',
      'HLT_HT550_DisplacedDijet40_Inclusive_v*',
      'HLT_HT650_DisplacedDijet80_Inclusive_v*',
      'HLT_HT750_DisplacedDijet80_Inclusive_v*',
      'HLT_VBF_DisplacedJet40_Hadronic_v*',
      'HLT_VBF_DisplacedJet40_DisplacedTrack_v*',
      'HLT_VBF_DisplacedJet40_DisplacedTrack_v*',                                                  
      'HLT_VBF_DisplacedJet40_TightID_DisplacedTrack_v*',                                          
      'HLT_VBF_DisplacedJet40_TightID_Hadronic_v*',                                                
      'HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack_v*',                                         
      'HLT_VBF_DisplacedJet40_VTightID_Hadronic_v*',                                               
      'HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack_v*',                                        
      'HLT_VBF_DisplacedJet40_VVTightID_Hadronic_v*',   
      'HLT_HT200_v*',
      'HLT_HT275_v*',
      'HLT_HT325_v*',
      'HLT_HT425_v*',
      'HLT_HT575_v*',
      'HLT_PFHT800_v*',
      'HLT_L1_TripleJet_VBF_v*',                                                                  
      'HLT_Mu20_v*',                                                                  
      'HLT_PFMET170_v*',                                                                  
      'HLT_PFMET170_NoiseCleaned_v*'                                                                  
      ),
    hltResults = cms.InputTag( "TriggerResults", "", trigger_process ),
    l1tResults = cms.InputTag( "" ),
    l1tIgnoreMask = cms.bool( False ),
    l1techIgnorePrescales = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True )
)


process.mcTriggerSelection = cms.EDFilter( "TriggerResultsFilter",
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
#      'HLT_L1_TripleJet_VBF_v*',                                                                  
#      'HLT_Mu20_v*',                                                                  
#      'HLT_PFMET170_v*',                                                                  
#      'HLT_PFMET170_NoiseCleaned_v*'                                                                  
#      'HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack_v*',                                         
#      'HLT_VBF_DisplacedJet40_VTightID_Hadronic_v*',                                               
#      'HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack_v*',                                        
#      'HLT_VBF_DisplacedJet40_VVTightID_Hadronic_v*',   
#      'HLT_PFHT200_v*',
#      'HLT_PFHT250_v*',
#      'HLT_PFHT300_v*',
#      'HLT_PFHT350_v*',
#      'HLT_PFHT400_v*',
#      'HLT_PFHT800_v*',
#      'HLT_L1_TripleJet_VBF_v*'                                                                  
      ),
    hltResults = cms.InputTag( "TriggerResults", "", trigger_process ),
    l1tResults = cms.InputTag( "" ),
    l1tIgnoreMask = cms.bool( False ),
    l1techIgnorePrescales = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True )
)

# trigger bits to keep
process.singleMuTrigger = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring(
      'HLT_Mu20_v*'                                                                 
      ),
    hltResults = cms.InputTag( "TriggerResults", "", trigger_process ),
    l1tResults = cms.InputTag( "" ),
    l1tIgnoreMask = cms.bool( False ),
    l1techIgnorePrescales = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True )
)

if doApplyTrigger: #apply the triggers and run dj tagging
   # if doApplySingleMu: 
   #    process.p *= process.doApplySingleMu
   if isMC:
      process.p *= process.mcTriggerSelection *  process.correctJets * process.djtagging
   else:
      process.p *= process.triggerSelection *  process.correctJets * process.djtagging
else: #just run the tagging sequence
   process.p *=  process.correctJets * process.djtagging

if doedm: #just dump the edm output of the djtagging sequence no analyzer
    process.btag_output = cms.EndPath(process.test_output)
else: # run the analysis and tree dumper using the djtagging sequence output
    process.p *= process.analyzerCALO

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool( True )
)
