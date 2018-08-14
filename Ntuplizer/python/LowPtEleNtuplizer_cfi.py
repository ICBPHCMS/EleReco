import FWCore.ParameterSet.Config as cms

LowPtEleNtuplizer = cms.EDAnalyzer("LowPtEleNtuplizer",
                                   hepMCProduct = cms.InputTag("generatorSmeared"),
                                   pileup = cms.InputTag("addPileupInfo"),
                                   vertices = cms.InputTag("offlinePrimaryVertices"),
                                   genParticles = cms.InputTag("genParticles"),
                                   generalTracks = cms.InputTag("generalTracks"),
                                   preIds = cms.InputTag("trackerDrivenElectronSeeds","preid"), # recoPreIds_trackerDrivenElectronSeeds_preid_RECO
                                   electronSeeds = cms.InputTag("electronMergedSeeds"), # recoElectronSeeds_electronMergedSeeds__RECO
                                   gsfTracks = cms.InputTag("electronGsfTracks"),
                                   gsfElectrons = cms.InputTag("gedGsfElectrons"),
                                   recoMuons = cms.InputTag("muons"),
                                   )
