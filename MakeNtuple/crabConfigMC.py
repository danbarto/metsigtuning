import os
from CRABClient.UserUtilities import getUsernameFromSiteDB
from WMCore.Configuration import Configuration
config = Configuration()

# Set variables
production_label = os.environ["CRAB_PROD_LABEL"]
dataset=os.environ["CRAB_DATASET"]
unitsPerJob = int(os.environ["CRAB_UNITS_PER_JOB"])
if "CRAB_TOTAL_UNITS" in os.environ: totalUnits = os.environ["CRAB_TOTAL_UNITS"]

config.section_("General")
config.General.transferLogs = True
config.General.requestName = production_label
config.General.workArea = config.General.requestName

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'makeEleNtuple_cfg.py'
config.JobType.psetName = 'makentuple_cfg.py'
config.JobType.outputFiles = ['ntuple.root']
config.JobType.inputFiles = ['Summer16_23Sep2016AllV4_DATA.db','Summer16_23Sep2016V4_MC.db', 'Spring16_25nsV6_DATA.db', 'Spring16_25nsV6_MC.db']
#config.JobType.inputFiles = ['Spring16_25nsV6_DATA.db','Spring16_25nsV6_MC.db']

config.section_("Data")
config.Data.inputDataset = dataset

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
if "IS_DATA" in os.environ:
    config.Data.lumiMask = 'json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.ignoreLocality = False

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.unitsPerJob = unitsPerJob#10
if "CRAB_TOTAL_UNITS" in os.environ: config.Data.totalUnits = int(totalUnits)#8

config.section_("Site")
#config.Site.blacklist = ['T2_US_Purdue', 'T2_US_Nebraska', 'T2_US_MIT', 'T2_US_Caltech']
config.Site.storageSite = 'T2_AT_Vienna'

config.section_("User")

config.section_("Debug")

