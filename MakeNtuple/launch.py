import imp, os, sys
from optparse import OptionParser
import re

# datasets to run as defined from run_susyMT2.cfg
# number of jobs to run per dataset decided based on splitFactor and fineSplitFactor from cfg file
# in principle one only needs to modify the following two lines:

parser = OptionParser(usage="python launch.py [options] component1 [ component2 ...]", \
                          description="Launch heppy jobs with CRAB3. Components correspond to the variables defined in heppy_samples.py (their name attributes)")
parser.add_option("--production_label", dest="production_label", help="production label", default="heppy")
parser.add_option("--remoteDir", dest="remoteDir", help="remote subdirectory", default="")
parser.add_option("--unitsPerJob", dest="unitsPerJob", help="Nr. of units (files) / crab job", type="int", default=1)
parser.add_option("--totalUnits", dest="totalUnits", help="Total nr. of units (files)", type="int", default=None)
parser.add_option("--inputDBS", dest="inputDBS", help="dbs instance", default=None)
parser.add_option("--lumiMask", dest="lumiMask", help="lumi mask (for data)", default=None)
parser.add_option("--dataset", dest="dataset", help="DAS dataset", default="/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM")
parser.add_option("--isData", action='store_true', dest="isData", help="lumi mask (for data)", default=False)
( options, args ) = parser.parse_args()


os.system("scram runtime -sh")
os.system("source /cvmfs/cms.cern.ch/crab3/crab.sh")

#os.environ["CMG_REMOTE_DIR"]  = options.remoteDir
os.environ["CRAB_UNITS_PER_JOB"] = str(options.unitsPerJob)
#os.environ["CMG_LUMI_MASK"] = options.lumiMask if options.lumiMask else "None"
if options.totalUnits:
    os.environ["CRAB_TOTAL_UNITS"] = str(options.totalUnits)

m=re.match("\/(.*)\/(.*)\/(.*)",options.dataset)
#print m.group(1), m.group(2)
if options.isData:
    os.environ["CRAB_PROD_LABEL"]  = options.production_label+"_"+m.group(1)+"_"+m.group(2)
else:
    os.environ["CRAB_PROD_LABEL"]  = options.production_label+"_"+m.group(1)#+"_"+m.group(2)

if options.isData:
    os.environ["IS_DATA"] = "True"

#print os.environ["CRAB_PROD_LABEL"]

os.environ["CRAB_DATASET"] = options.dataset

os.system("which crab")
os.system("crab submit -c crabConfigMC.py")
