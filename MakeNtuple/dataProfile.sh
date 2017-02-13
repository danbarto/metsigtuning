#!/bin/sh

PILEUP_LATEST=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
JSON=json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
LUMI=36000

if [ ! -f "$PILEUP_LATEST" ]; then
   echo "File $PILEUP_LATEST does not exist on this site, copying from lxplus"
   scp $USER@lxplus.cern.ch:$PILEUP_LATEST pileup_latest.txt
   PILEUP_LATEST=pileup_latest.txt
fi


#pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 62280 --maxPileupBin 50 --numPileupBins 50 PU_2016_${LUMI}_XSecVDown.root
#pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 65740 --maxPileupBin 100 --numPileupBins 100 PU_2016_${LUMI}_XSecDown.root
pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 69200 --maxPileupBin 75 --numPileupBins 75 PU_2016_${LUMI}_XSecCentral.root
#pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 72660 --maxPileupBin 100 --numPileupBins 100 PU_2016_${LUMI}_XSecUp.root
#pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 76120 --maxPileupBin 50 --numPileupBins 50 PU_2016_${LUMI}_XSecVUp.root
#pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 79580 --maxPileupBin 50 --numPileupBins 50 PU_2016_${LUMI}_XSecVVUp.root

