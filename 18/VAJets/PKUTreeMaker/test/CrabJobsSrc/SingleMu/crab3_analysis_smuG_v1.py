from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = '94x_smuG_v1'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles =['Summer16_07Aug2017GH_V11_DATA_L1FastJet_AK4PFchs.txt','Summer16_07Aug2017GH_V11_DATA_L2Relative_AK4PFchs.txt','Summer16_07Aug2017GH_V11_DATA_L3Absolute_AK4PFchs.txt','Summer16_07Aug2017GH_V11_DATA_L2L3Residual_AK4PFchs.txt','Summer16_07Aug2017GH_V11_DATA_L1FastJet_AK4PFPuppi.txt','Summer16_07Aug2017GH_V11_DATA_L2Relative_AK4PFPuppi.txt','Summer16_07Aug2017GH_V11_DATA_L3Absolute_AK4PFPuppi.txt','Summer16_07Aug2017GH_V11_DATA_L2L3Residual_AK4PFPuppi.txt']

config.JobType.psetName    = 'analysis_data_G.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
#config.Data.outputPrimaryDataset = 'VBS_WGAMMA_94X'
config.Data.inputDataset = '/SingleMuon/Run2016G-17Jul2018-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 40
config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.publication = False
config.Data.outputDatasetTag = '94x_smuG_v1'

config.section_("Site")
config.Site.storageSite = 'T2_CN_Beijing'
