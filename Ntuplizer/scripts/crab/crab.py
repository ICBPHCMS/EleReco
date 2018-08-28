from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

import argparse
parser = argparse.ArgumentParser(description='Arguments to crab.py file')
parser.add_argument('-n','--name', help='Name of input cfg file, output root file, etc', default='EleReco-Gsf0p5', required=False)
parser.add_argument('-d','--dataset', help='DAS dataset name', default='/BToKee_Pythia/tstreble-BToKee_Pythia_PUMix_18_03_18-c9b9e020b5bce5ee6bee9ef5f38c415a/USER', required=False)
parser.add_argument('-i','--instance', help='DBS instance', default='phys03', required=False)
args = vars(parser.parse_args())

config.General.requestName     = args.name
config.General.workArea        = args.name
config.General.transferOutputs = True
config.General.transferLogs    = True

config.JobType.pluginName = 'ANALYSIS'
config.JobType.psetName   = args.name+'_cfg.py'
config.JobType.outputFiles = [args.name+'.root']

config.Data.inputDataset         = args.dataset
config.Data.inputDBS             = args.instance
config.Data.splitting            = 'FileBased' #Automatic,EventAwareLumiBased
config.Data.unitsPerJob          = 1
config.Data.totalUnits           = -1
config.Data.outLFNDirBase        = '/store/user/%s/'%(getUsernameFromSiteDB())
config.Data.publication          = False
config.Data.outputDatasetTag     = args.name

config.Site.whitelist   = ["T2_UK_London_IC"]
config.Site.storageSite = 'T2_UK_London_IC'
