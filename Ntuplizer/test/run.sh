#INPUT=file:/home/hep/rjb3/work/tracker_driven_seeds/EleReco-Seed2p0/CMSSW_9_4_8/src/EleReco/Ntuplizer/scripts/crab/EleReco_Seed2p0.root

SUFFIX=5b
INPUT=file:/home/hep/rjb3/work/tracker_driven_seeds/EleReco-Open/CMSSW_9_4_8/src/EleReco/Ntuplizer/scripts/crab/EleReco-Open-$SUFFIX.root
OUTPUT=output.root
NEVENTS=-1

cmsRun LowPtEleNtuplizer_cfg.py inputFiles=$INPUT outputFile=$OUTPUT maxEvents=$NEVENTS
