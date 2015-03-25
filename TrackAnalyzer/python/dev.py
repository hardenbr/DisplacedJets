import FWCore.ParameterSet.Config as cms

isSignalMC = True
isMC = True
doedm = True
nevents = 30

input_file = None

# test files
if isSignalMC:
    input_file = 'file:/afs/cern.ch/work/h/hardenbr/2015/DIJET/GEN_SIGNAL_TEST/dijet_700_300_ctau30.root'
#    input_file = 'file:/afs/cern.ch/work/h/hardenbr/2015/DIJET/GEN_SIGNAL_TEST/dijet_700_300_ctau300.root'
#   input_file = 'file:/afs/cern.ch/work/h/hardenbr/2015/DIJET/GEN_SIGNAL_TEST/dijet_700_300_ctau3000.root'
#   input_file = 'file:/afs/cern.ch/work/h/hardenbr/2015/DIJET/GEN_SIGNAL_TEST/dijet_700_300_ctau3.root'
#   input_file = 'file:/afs/cern.ch/user/t/tkolberg/public/hepmcreco_RAW2DIGI_RECO.root'
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
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_100_1_NtX.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_101_1_j8l.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_102_1_wex.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_103_1_YLq.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_104_1_vvt.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_105_1_eMQ.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_106_1_xxr.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_107_1_KEJ.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_108_1_A93.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_10_1_6Ew.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_110_1_AI4.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_111_1_tkS.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_112_1_PrO.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_113_1_SA1.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_114_1_IkF.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_115_1_l74.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_116_1_QCM.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_117_1_LbX.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_118_1_Exh.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_119_1_R9b.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_11_1_KP8.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_120_1_5Rl.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_12_1_DAb.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_13_1_bEK.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_14_1_bte.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_15_1_Ume.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_16_1_xao.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_17_1_PKE.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_18_1_Bat.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_19_1_hg8.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_1_1_xQK.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_20_1_GDn.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_21_1_1py.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_22_1_5mE.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_23_1_N24.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_24_1_xhT.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_25_1_IFp.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_26_1_bi3.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_27_1_f7A.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_28_1_Fli.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_29_1_kTt.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_2_1_G8a.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_30_1_MNI.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_31_1_Nwr.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_32_1_ps4.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_33_1_KyN.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_34_1_AXq.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_35_1_T0f.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_36_1_tAS.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_37_1_UVa.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_38_1_dyH.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_39_1_lEE.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_3_1_9nI.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_40_1_meC.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_41_1_O1F.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_42_1_DsK.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_43_1_8Wa.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_44_1_LD6.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_45_1_Ywr.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_46_1_wC2.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_47_1_e5y.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_48_1_09b.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_49_1_VsM.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_4_1_BcI.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_50_1_tbJ.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_51_1_W5o.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_52_1_ICZ.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_53_1_TVt.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_54_1_zSc.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_55_1_nmt.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_56_1_R0u.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_57_1_lw7.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_58_1_zc3.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_59_1_Z6f.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_5_1_aEU.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_60_1_8XJ.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_61_1_at2.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_62_1_FYf.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_63_1_X7E.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_64_1_GeN.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_65_1_dW8.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_66_1_bRT.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_67_1_d4P.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_68_1_2cl.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_69_1_T2o.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_6_1_oCB.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_70_1_lIe.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_71_1_Lhy.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_72_1_LlY.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_73_1_aJU.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_74_1_m4M.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_75_1_Ibt.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_76_1_Vca.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_77_1_vcH.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_78_1_jWq.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_79_1_7V6.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_7_1_FHP.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_80_1_Yhd.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_81_1_Jyk.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_82_1_NEg.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_83_1_Acw.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_84_1_UoT.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_85_1_Tf3.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_86_1_71a.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_87_1_x4D.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_88_1_YXY.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_89_1_KHe.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_8_1_H07.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_90_1_Vms.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_91_1_yuI.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_92_1_xLk.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_93_1_Dl9.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_94_1_JnO.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_95_1_WUS.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_96_1_yZA.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_97_1_W9n.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_99_1_j4a.root',
      # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU100_40BX25_AOD/dijet_700_300_ctau100_9_1_egE.root'
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
process.analyzer.secondaryVertexTagInfo = cms.untracked.InputTag('displacedSecondaryVertexTagInfosNoPV', '', 'ANA')
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
    process.displacedTagsToVertices.doGenMatch = cms.untracked.bool(False)

if doedm:
    process.p = cms.Path(process.djtagging)    
    process.btag_output = cms.EndPath( process.test_output)
else:
    process.p = cms.Path(process.djtagging + process.analyzer)
