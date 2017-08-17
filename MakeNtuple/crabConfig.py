import os
from CRABClient.UserUtilities import getUsernameFromSiteDB
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs = True
config.General.requestName = 'MetSig_DataMoriond17_2016H-v3_v3'
config.General.workArea = config.General.requestName

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'makentuple_cfg.py'
config.JobType.outputFiles = ['ntuple.root']
config.JobType.inputFiles = ['Summer16_23Sep2016AllV4_DATA.db','Summer16_23Sep2016V4_MC.db', 'Spring16_25nsV6_DATA.db', 'Spring16_25nsV6_MC.db']

config.section_("Data")
config.Data.lumiMask = 'json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

#config.Data.inputDataset = '/DoubleMuon/Run2016B-23Sep2016-v3/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016C-23Sep2016-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016D-23Sep2016-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016E-23Sep2016-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016F-23Sep2016-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016G-23Sep2016-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016H-PromptReco-v2/MINIAOD'
config.Data.inputDataset = '/DoubleMuon/Run2016H-PromptReco-v3/MINIAOD'


config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased' #FIXME check!
config.Data.ignoreLocality = False

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.unitsPerJob = 5 #was 50 for lumibased
#config.Data.totalUnits = 10

config.section_("Site")
#config.Site.blacklist = ['T2_US_Purdue', 'T2_US_Nebraska', 'T2_US_MIT', 'T2_US_Caltech']
config.Site.storageSite = 'T2_AT_Vienna'

config.section_("User")

config.section_("Debug")

