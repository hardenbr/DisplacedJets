from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName     = 'DATA_2015D_DEC26'
config.General.workArea        = 'CRAB3_DATA_RUN2015D_DEC26'
config.General.transferOutputs = True
config.General.transferLogs    = True

config.JobType.pluginName      = 'Analysis'
config.JobType.psetName        = 'dev.py'
config.JobType.outputFiles     = ['data.root']
#silver json
config.Data.lumiMask         = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.txt'
#config.Data.lumiMask         = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
#config.Data.lumiMask         = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY.txt'
#config.Data.lumiMask         = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'
#config.Data.lumiMask         = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'
#config.Data.lumiMask         = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt'
#config.Data.lumiMask         = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt'
#config.Data.lumiMask         = 'Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'

config.Data.inputDBS         = 'global'
config.Data.splitting        = 'FileBased'
config.Data.unitsPerJob      = 3
config.Data.totalUnits       = -1 #config.Data.unitsPerJob * NJOBS
#config.Data.runRange        = '193093-193999' # '193093-194075'
config.Data.outLFNDirBase    = '/store/group/phys_susy/razor/josh/RAZOR_DIJET/DJANALYSIS/RUN2015D_DEC26_JSON_SILVER_v2_NOPRESEL'
config.Data.publication      = False
config.Data.ignoreLocality   = True #seth zenz fix 
#config.Data.publishDataName = 'CRAB3_tutorial_May2015_Data_analysis'
config.Site.storageSite      = 'T2_CH_CERN'

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

#     config.General.requestName = 'PD_2015C_DisplacedJet'
# #    config.Data.inputDataset = '/DisplacedJet/Run2015B-PromptReco-v1/AOD'
#     config.Data.inputDataset = '/DisplacedJet/Run2015C-PromptReco-v1/AOD'
#     # config.Data.unitsPerJob = 2
#     #config.Data.totalUnits = 4
#     submit(config)

    config.General.requestName = 'PD_2015D_DisplacedJet'
    config.Data.inputDataset   = '/DisplacedJet/Run2015D-PromptReco-v3/AOD'
    # config.Data.unitsPerJob  = 2
    #config.Data.totalUnits    = 4
    submit(config)

    # config.General.requestName = 'PD_2015C_JetHT'
    # config.Data.inputDataset   = '/JetHT/Run2015C-PromptReco-v1/AOD'
    # # config.Data.unitsPerJob  = 2
    # #                          = config.Data.totalUnits = 4
    # submit(config)

    config.General.requestName = 'PD_2015D_JetHT'
    config.Data.inputDataset   = '/JetHT/Run2015D-PromptReco-v3/AOD'
    # config.Data.unitsPerJob  = 2
    #                          = config.Data.totalUnits = 4
    submit(config)

    # config.General.requestName = 'PD_MET'
    # config.Data.inputDataset   = '/MET/Run2015B-PromptReco-v1/AOD'
    # # config.Data.unitsPerJob  = 2
    # #                          = config.Data.totalUnits = 4
    # submit(config)

    # config.General.requestName = 'PD_2015D_SingleMuon'
    # config.Data.inputDataset   = '/SingleMuon/Run2015D-PromptReco-v3/AOD'
    # # config.Data.unitsPerJob  = 2
    # #                          = config.Data.totalUnits = 4
    # submit(config)

    # config.General.requestName = 'PD_HTMHT'
    # config.Data.inputDataset   = '/HTMHT/Run2015B-PromptReco-v1/AOD'
    # submit(config)
   
