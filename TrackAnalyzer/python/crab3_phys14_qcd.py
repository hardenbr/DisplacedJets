from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'QCD_PHYS14'
config.General.workArea = 'CRAB3_QCD_PHYS14_4UNITS'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'dev.py'
config.JobType.outputFiles = ['qcd.root']

#config.Data.inputDataset = '/MinimumBias/Commissioning2015-v1/RAW'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 4
config.Data.totalUnits = -1 #config.Data.unitsPerJob * NJOBS
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#config.Data.runRange = '193093-193999' # '193093-194075'
config.Data.outLFNDirBase = '/store/group/phys_susy/razor/josh/RAZOR_DIJET/DJANALYSIS/QCD/'
config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_tutorial_May2015_Data_analysis'
config.Site.storageSite = 'T2_CH_CERN'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    #config.General.workArea = 'crab_projects'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    config.General.requestName = 'qcd15_30'
    config.Data.inputDataset = '/QCD_Pt-15to30_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v2/AODSIM'
    # config.Data.unitsPerJob = 2
    # =config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd30_50'
    config.Data.inputDataset = '/QCD_Pt-30to50_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v2/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd50_80'
    config.Data.inputDataset = '/QCD_Pt-50to80_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd80_120'
    config.Data.inputDataset = '/QCD_Pt-80to120_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd120_170'
    config.Data.inputDataset = '/QCD_Pt-120to170_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd170_300'
    config.Data.inputDataset = '/QCD_Pt-170to300_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd300_470'
    config.Data.inputDataset = '/QCD_Pt-300to470_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd470_600'
    config.Data.inputDataset = '/QCD_Pt-470to600_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd600_800'
    config.Data.inputDataset = '/QCD_Pt-600to800_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd800_1000'
    config.Data.inputDataset = '/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd1000_1400'
    config.Data.inputDataset = '/QCD_Pt-1000to1400_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd1400_1800'
    config.Data.inputDataset = '/QCD_Pt-1400to1800_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd1800'
    config.Data.inputDataset = '/QCD_Pt-1800_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)
