#ifndef LOWPTELENTUPLIZER_H
#define LOWPTELENTUPLIZER_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <FWCore/Utilities/interface/InputTag.h>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <TNtuple.h>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>

//#include <algorithm>
//#include <map>
//#include <utility>
//#include <TNtuple.h>
//#include <TString.h>
//#include <bitset>

/*******************************************************************************
 * Class declaration 
 */
class LowPtEleNtuplizer : public edm::EDAnalyzer {

public:

  explicit LowPtEleNtuplizer(const edm::ParameterSet&);
  virtual ~LowPtEleNtuplizer();
  
private:

  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  void init();
  float deltaPhi( float phi1, float phi2 );
  float deltaEta( float eta1, float eta2 );
  float deltaR( float eta1, float eta2, float phi1, float phi2 );

  static bool genSortByPt(const reco::GenParticle* i,
			  const reco::GenParticle* j) {return (i->pt()>j->pt());}
  static bool trackSortByPt(const reco::Track* i,
			    const reco::Track* j) {return (i->pt()>j->pt());}
  static bool recoSortByPt(const reco::RecoCandidate* i,
			   const reco::RecoCandidate* j) {return (i->pt()>j->pt());}
  
  TTree* _tree;
  
  edm::EDGetTokenT<edm::HepMCProduct> _hepMCProductTag;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> > _pileupTag;
  edm::EDGetTokenT< std::vector<reco::Vertex> > _verticesTag;
  edm::EDGetTokenT< std::vector<reco::GenParticle> > _genParticlesTag;
  edm::EDGetTokenT< std::vector<reco::Track> > _generalTracksTag;
  edm::EDGetTokenT< std::vector<reco::GsfTrack> > _gsfTracksTag;
  edm::EDGetTokenT< std::vector<reco::GsfElectron> > _gsfElectronsTag;
  edm::EDGetTokenT< std::vector<reco::Muon> > _recoMuonsTag;

  // variables

  Int_t _run;
  Int_t _lumi;
  Long64_t _event;
  Int_t _pileup;
  Int_t _vertices;

  Double_t _genLabDispl;
  Double_t _genProperDispl;

  Int_t _genParticlesN;

  Int_t _genElesN;
  std::vector<Double_t> _genElesPt;
  std::vector<Double_t> _genElesEta;
  std::vector<Double_t> _genElesPhi;

  Int_t _genMuonsN;
  std::vector<Double_t> _genMuonsPt;
  std::vector<Double_t> _genMuonsEta;
  std::vector<Double_t> _genMuonsPhi;

  Int_t _genPionsN;
  std::vector<Double_t> _genPionsPt;
  std::vector<Double_t> _genPionsEta;
  std::vector<Double_t> _genPionsPhi;

  Int_t _genKaonsN;
  std::vector<Double_t> _genKaonsPt;
  std::vector<Double_t> _genKaonsEta;
  std::vector<Double_t> _genKaonsPhi;

  Int_t _genTrksN; // general tracks
  std::vector<Double_t> _genTrksPt; // general tracks
  std::vector<Double_t> _genTrksEta; // general tracks
  std::vector<Double_t> _genTrksPhi; // general tracks

  Int_t _gsfTrksN;
  std::vector<Double_t> _gsfTrksPt;
  std::vector<Double_t> _gsfTrksEta;
  std::vector<Double_t> _gsfTrksPhi;

  Int_t _gsfElesN;
  std::vector<Double_t> _gsfElesPt;
  std::vector<Double_t> _gsfElesEta; // extrapolated for electron, consider _gsfTrkEta for eta based on just GsfTrack
  std::vector<Double_t> _gsfElesPhi; // extrapolated for electron, consider _gsfTrkPhi for phi based on just GsfTrack
  std::vector<Double_t> _gsfEles_ctfTrkPt; // closest ctf trk 
  std::vector<Double_t> _gsfEles_ctfTrkEta; // closest ctf trk 
  std::vector<Double_t> _gsfEles_ctfTrkPhi; // closest ctf trk 
  std::vector<Double_t> _gsfEles_gsfTrkPt; // closest gsf trk 
  std::vector<Double_t> _gsfEles_gsfTrkEta; // closest gsf trk
  std::vector<Double_t> _gsfEles_gsfTrkPhi; // closest gsf trk
  std::vector<Int_t> _gsfEles_trackerDriven;
  std::vector<Int_t> _gsfEles_ecalDriven;

  Int_t _recoMuonsN;
  std::vector<Double_t> _recoMuonsPt;
  std::vector<Double_t> _recoMuonsEta;
  std::vector<Double_t> _recoMuonsPhi;
  std::vector<Int_t> _recoMuons_global;
  std::vector<Int_t> _recoMuons_tracker;
  std::vector<Int_t> _recoMuons_pf;
  
  std::vector<Int_t> _genEles_genTrksIdx;
  std::vector<Double_t> _genEles_genTrksDR;
  std::vector<Double_t> _genEles_genTrksDXY;
  std::vector<Double_t> _genEles_genTrksDZ;

  std::vector<Int_t> _genEles_gsfTrksIdx;
  std::vector<Double_t> _genEles_gsfTrksDR;
  std::vector<Double_t> _genEles_gsfTrksDXY;
  std::vector<Double_t> _genEles_gsfTrksDZ;

  std::vector<Int_t> _genEles_gsfElesIdx;
  std::vector<Double_t> _genEles_gsfElesDR;
  std::vector<Double_t> _genEles_gsfElesDXY;
  std::vector<Double_t> _genEles_gsfElesDZ;

  std::vector<Int_t> _genMuons_genTrksIdx;
  std::vector<Double_t> _genMuons_genTrksDR;
  std::vector<Double_t> _genMuons_genTrksDXY;
  std::vector<Double_t> _genMuons_genTrksDZ;

  std::vector<Int_t> _genMuons_recoMuonsIdx;
  std::vector<Double_t> _genMuons_recoMuonsDR;
  std::vector<Double_t> _genMuons_recoMuonsDXY;
  std::vector<Double_t> _genMuons_recoMuonsDZ;

  std::vector<Int_t> _genKaons_genTrksIdx;
  std::vector<Double_t> _genKaons_genTrksDR;
  std::vector<Double_t> _genKaons_genTrksDXY;
  std::vector<Double_t> _genKaons_genTrksDZ;
  //Float_t
  //Bool_t

  Int_t _genEleLead_genTrksN; // number of general tracks matched to leading gen electron
  std::vector<Int_t> _genEleLead_genTrksIdx; // index of general track matched to leading gen electron 
  std::vector<Double_t> _genEleLead_genTrksDR; // deltaR of matched general track
  std::vector<Double_t> _genEleLead_genTrksDEta; // deltaEta of matched general track
  std::vector<Double_t> _genEleLead_genTrksDPhi; // deltaPhi of matched general track
  std::vector<Double_t> _genEleLead_genTrksPtRatio; // pT ratio (RECO/GEN) of matched general track
  std::vector<Double_t> _genEleLead_genTrksDXY; // deltaXY of vertex of matched general track
  std::vector<Double_t> _genEleLead_genTrksDZ; // deltaZ of vertex of matched general track
  std::vector<Double_t> _genEleLead_genTrksPxlHits; // number of pixel hits associated to matched general track
  std::vector<Double_t> _genEleLead_genTrksStrHits; // number of strip hits associated to matched general track
  Float_t _genEleLead_genPartMinDR; // minDR between leading gen electron and other gen particles 

  Int_t _genEleSub_genTrksN; // number of general tracks matched to subleading gen electron
  std::vector<Int_t> _genEleSub_genTrksIdx; // index of general track matched to subleading gen electron 
  std::vector<Double_t> _genEleSub_genTrksDR; // deltaR of matched general track
  std::vector<Double_t> _genEleSub_genTrksDEta; // deltaEta of matched general track
  std::vector<Double_t> _genEleSub_genTrksDPhi; // deltaPhi of matched general track
  std::vector<Double_t> _genEleSub_genTrksPtRatio; // pT ratio (RECO/GEN) of matched general track
  std::vector<Double_t> _genEleSub_genTrksDXY; // deltaXY of vertex of matched general track
  std::vector<Double_t> _genEleSub_genTrksDZ; // deltaZ of vertex of matched general track
  std::vector<Double_t> _genEleSub_genTrksPxlHits; // number of pixel hits associated to matched general track
  std::vector<Double_t> _genEleSub_genTrksStrHits; // number of strip hits associated to matched general track
  Float_t _genEleSub_genPartMinDR; // minDR between subleading gen electron and other gen particles 
  
};

//******************************************************************************
//******************************************************************************
//******************************************************************************

/*******************************************************************************
 * Constructor
 */
LowPtEleNtuplizer::LowPtEleNtuplizer(const edm::ParameterSet& iConfig) :
  _hepMCProductTag( consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("hepMCProductLabel")) ),
  _pileupTag( consumes< std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupLabel")) ),
  _verticesTag( consumes< std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("verticesLabel")) ),
  _genParticlesTag( consumes< std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticlesLabel")) ),
  _generalTracksTag( consumes< std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("generalTracksLabel")) ),
  _gsfTracksTag( consumes< std::vector<reco::GsfTrack> >(iConfig.getParameter<edm::InputTag>("gsfTracksLabel")) ),
  _gsfElectronsTag( consumes< std::vector<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("gsfElectronsLabel")) ),
  _recoMuonsTag( consumes< std::vector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("recoMuonsLabel")) ),
  _run(-1),
  _lumi(-1),
  _event(-1),
  _pileup(-1),
  _vertices(-1),
  _genLabDispl(-1.),
  _genProperDispl(-1.),
  _genParticlesN(0),
  _genElesN(0),
  _genElesPt(),  
  _genElesEta(),  
  _genElesPhi(),  
  _genMuonsN(0),
  _genMuonsPt(),  
  _genMuonsEta(),  
  _genMuonsPhi(),  
  _genPionsN(0),
  _genPionsPt(),  
  _genPionsEta(),  
  _genPionsPhi(),  
  _genKaonsN(0),
  _genKaonsPt(),  
  _genKaonsEta(),  
  _genKaonsPhi(),  
  _genTrksN(0),
  _genTrksPt(),  
  _genTrksEta(),  
  _genTrksPhi(),
  _gsfTrksN(0),
  _gsfTrksPt(),  
  _gsfTrksEta(),  
  _gsfTrksPhi(),
  _gsfElesN(0),
  _gsfElesPt(),  
  _gsfElesEta(),  
  _gsfElesPhi(),
  _gsfEles_ctfTrkPt(),  
  _gsfEles_ctfTrkEta(),  
  _gsfEles_ctfTrkPhi(),
  _gsfEles_gsfTrkPt(),  
  _gsfEles_gsfTrkEta(),  
  _gsfEles_gsfTrkPhi(),
  _gsfEles_trackerDriven(),
  _gsfEles_ecalDriven(),
  _recoMuonsN(0),
  _recoMuonsPt(),  
  _recoMuonsEta(),  
  _recoMuonsPhi(),
  _recoMuons_global(),
  _recoMuons_tracker(),
  _recoMuons_pf(),
  _genEles_genTrksIdx(),
  _genEles_genTrksDR(),
  _genEles_genTrksDXY(),
  _genEles_genTrksDZ(),
  _genEles_gsfTrksIdx(),
  _genEles_gsfTrksDR(),
  _genEles_gsfTrksDXY(),
  _genEles_gsfTrksDZ(),
  _genEles_gsfElesIdx(),
  _genEles_gsfElesDR(),
  _genEles_gsfElesDXY(),
  _genEles_gsfElesDZ(),
  _genMuons_genTrksIdx(),
  _genMuons_genTrksDR(),
  _genMuons_genTrksDXY(),
  _genMuons_genTrksDZ(),
  _genMuons_recoMuonsIdx(),
  _genMuons_recoMuonsDR(),
  _genMuons_recoMuonsDXY(),
  _genMuons_recoMuonsDZ(),
  _genKaons_genTrksIdx(),
  _genKaons_genTrksDR(),
  _genKaons_genTrksDXY(),
  _genKaons_genTrksDZ(),
  _genEleLead_genTrksN(0),
  _genEleLead_genTrksIdx(),
  _genEleLead_genTrksDR(),
  _genEleLead_genTrksDEta(),
  _genEleLead_genTrksDPhi(),
  _genEleLead_genTrksPtRatio(),
  _genEleLead_genTrksDXY(),
  _genEleLead_genTrksDZ(),
  _genEleLead_genTrksPxlHits(),
  _genEleLead_genTrksStrHits(),
  _genEleLead_genPartMinDR(-1.),
  _genEleSub_genTrksN(0),
  _genEleSub_genTrksIdx(),
  _genEleSub_genTrksDR(),
  _genEleSub_genTrksDEta(),
  _genEleSub_genTrksDPhi(),
  _genEleSub_genTrksPtRatio(),
  _genEleSub_genTrksDXY(),
  _genEleSub_genTrksDZ(),
  _genEleSub_genTrksPxlHits(),
  _genEleSub_genTrksStrHits(),
  _genEleSub_genPartMinDR(-1.)
{
  init();
}

/*******************************************************************************
 * Constructor
 */
LowPtEleNtuplizer::~LowPtEleNtuplizer()
{}

/*******************************************************************************
 * 
 */
void LowPtEleNtuplizer::beginJob()
{

  // Create TTree file 
  edm::Service<TFileService> fs;
  _tree = fs->make<TTree>("tree","tree");
  
  // Branches
  _tree->Branch("RunNumber",&_run);
  _tree->Branch("LumiSection",&_lumi);
  _tree->Branch("EventNumber",&_event,"EventNumber/L");
  _tree->Branch("Pileup",&_pileup);
  _tree->Branch("PrimaryVertices",&_vertices);
  
  _tree->Branch("genLabDispl",&_genLabDispl);
  _tree->Branch("genProperDispl",&_genProperDispl);
  
  _tree->Branch("genParticles_N",&_genParticlesN,"genParticles_N/I");
  
  _tree->Branch("genEles_N",&_genElesN,"genEles_N/I");
  _tree->Branch("genEles_Pt",&_genElesPt);
  _tree->Branch("genEles_Eta",&_genElesEta);
  _tree->Branch("genEles_Phi",&_genElesPhi);

  _tree->Branch("genMuons_N",&_genMuonsN,"genMuons_N/I");
  _tree->Branch("genMuons_Pt",&_genMuonsPt);
  _tree->Branch("genMuons_Eta",&_genMuonsEta);
  _tree->Branch("genMuons_Phi",&_genMuonsPhi);

  _tree->Branch("genPions_N",&_genPionsN,"genPions_N/I");
  _tree->Branch("genPions_Pt",&_genPionsPt);
  _tree->Branch("genPions_Eta",&_genPionsEta);
  _tree->Branch("genPions_Phi",&_genPionsPhi);

  _tree->Branch("genKaons_N",&_genKaonsN,"genKaons_N/I");
  _tree->Branch("genKaons_Pt",&_genKaonsPt);
  _tree->Branch("genKaons_Eta",&_genKaonsEta);
  _tree->Branch("genKaons_Phi",&_genKaonsPhi);

  _tree->Branch("genTrks_N",&_genTrksN,"genTrks_N/I");
  _tree->Branch("genTrks_Pt",&_genTrksPt);
  _tree->Branch("genTrks_Eta",&_genTrksEta);
  _tree->Branch("genTrks_Phi",&_genTrksPhi);
  
  _tree->Branch("gsfTrks_N",&_gsfTrksN,"gsfTrks_N/I");
  _tree->Branch("gsfTrks_Pt",&_gsfTrksPt);
  _tree->Branch("gsfTrks_Eta",&_gsfTrksEta);
  _tree->Branch("gsfTrks_Phi",&_gsfTrksPhi);
  
  _tree->Branch("gsfEles_N",&_gsfElesN,"gsfEles_N/I");
  _tree->Branch("gsfEles_Pt",&_gsfElesPt);
  _tree->Branch("gsfEles_Eta",&_gsfElesEta);
  _tree->Branch("gsfEles_Phi",&_gsfElesPhi);
  _tree->Branch("gsfEles_ctfTrk_Pt", &_gsfEles_ctfTrkPt);
  _tree->Branch("gsfEles_ctfTrk_Eta",&_gsfEles_ctfTrkEta);
  _tree->Branch("gsfEles_ctfTrk_Phi",&_gsfEles_ctfTrkPhi);
  _tree->Branch("gsfEles_gsfTrk_Pt", &_gsfEles_gsfTrkPt);
  _tree->Branch("gsfEles_gsfTrk_Eta",&_gsfEles_gsfTrkEta);
  _tree->Branch("gsfEles_gsfTrk_Phi",&_gsfEles_gsfTrkPhi);
  _tree->Branch("gsfEles_trackerDriven",&_gsfEles_trackerDriven);
  _tree->Branch("gsfEles_ecalDriven",&_gsfEles_ecalDriven);
  
  _tree->Branch("recoMuons_N",&_recoMuonsN,"recoMuons_N/I");
  _tree->Branch("recoMuons_Pt",&_recoMuonsPt);
  _tree->Branch("recoMuons_Eta",&_recoMuonsEta);
  _tree->Branch("recoMuons_Phi",&_recoMuonsPhi);
  _tree->Branch("recoMuons_global",&_recoMuons_global);
  _tree->Branch("recoMuons_tracker",&_recoMuons_tracker);
  _tree->Branch("recoMuons_pf",&_recoMuons_pf);
  
  _tree->Branch("genEles_genTrks_Idx",&_genEles_genTrksIdx);
  _tree->Branch("genEles_genTrks_DR",&_genEles_genTrksDR);
  _tree->Branch("genEles_genTrks_DXY",&_genEles_genTrksDXY);
  _tree->Branch("genEles_genTrks_DZ",&_genEles_genTrksDZ);
  _tree->Branch("genEles_gsfTrks_Idx",&_genEles_gsfTrksIdx);
  _tree->Branch("genEles_gsfTrks_DR",&_genEles_gsfTrksDR);
  _tree->Branch("genEles_gsfTrks_DXY",&_genEles_gsfTrksDXY);
  _tree->Branch("genEles_gsfTrks_DZ",&_genEles_gsfTrksDZ);
  _tree->Branch("genEles_gsfEles_Idx",&_genEles_gsfElesIdx);
  _tree->Branch("genEles_gsfEles_DR",&_genEles_gsfElesDR);
  _tree->Branch("genEles_gsfEles_DXY",&_genEles_gsfElesDXY);
  _tree->Branch("genEles_gsfEles_DZ",&_genEles_gsfElesDZ);

  _tree->Branch("genMuons_genTrks_Idx",&_genMuons_genTrksIdx);
  _tree->Branch("genMuons_genTrks_DR",&_genMuons_genTrksDR);
  _tree->Branch("genMuons_genTrks_DXY",&_genMuons_genTrksDXY);
  _tree->Branch("genMuons_genTrks_DZ",&_genMuons_genTrksDZ);
  _tree->Branch("genMuons_recoMuons_Idx",&_genMuons_recoMuonsIdx);
  _tree->Branch("genMuons_recoMuons_DR",&_genMuons_recoMuonsDR);
  _tree->Branch("genMuons_recoMuons_DXY",&_genMuons_recoMuonsDXY);
  _tree->Branch("genMuons_recoMuons_DZ",&_genMuons_recoMuonsDZ);
  
  _tree->Branch("genKaons_genTrks_Idx",&_genKaons_genTrksIdx);
  _tree->Branch("genKaons_genTrks_DR",&_genKaons_genTrksDR);
  _tree->Branch("genKaons_genTrks_DXY",&_genKaons_genTrksDXY);
  _tree->Branch("genKaons_genTrks_DZ",&_genKaons_genTrksDZ);

  _tree->Branch("genEleLead_genTrks_N",&_genEleLead_genTrksN,"genEleLead_genTrks_N/I");
  _tree->Branch("genEleLead_genTrks_Idx",&_genEleLead_genTrksIdx);
  _tree->Branch("genEleLead_genTrks_DR",&_genEleLead_genTrksDR);
  _tree->Branch("genEleLead_genTrks_DEta",&_genEleLead_genTrksDEta);
  _tree->Branch("genEleLead_genTrks_DPhi",&_genEleLead_genTrksDPhi);
  _tree->Branch("genEleLead_genTrks_PtRatio",&_genEleLead_genTrksPtRatio);
  _tree->Branch("genEleLead_genTrks_DXY",&_genEleLead_genTrksDXY);
  _tree->Branch("genEleLead_genTrks_DZ",&_genEleLead_genTrksDZ);
  _tree->Branch("genEleLead_genTrks_PxlHits",&_genEleLead_genTrksPxlHits);
  _tree->Branch("genEleLead_genTrks_StrHits",&_genEleLead_genTrksStrHits);
  _tree->Branch("genEleLead_genPart_MinDR",&_genEleLead_genPartMinDR,"genEleLead_genPart_MinDR/I");

  _tree->Branch("genEleSub_genTrks_N",&_genEleSub_genTrksN,"genEleSub_genTrks_N/I");
  _tree->Branch("genEleSub_genTrks_Idx",&_genEleSub_genTrksIdx);
  _tree->Branch("genEleSub_genTrks_DR",&_genEleSub_genTrksDR);
  _tree->Branch("genEleSub_genTrks_DEta",&_genEleSub_genTrksDEta);
  _tree->Branch("genEleSub_genTrks_DPhi",&_genEleSub_genTrksDPhi);
  _tree->Branch("genEleSub_genTrks_PtRatio",&_genEleSub_genTrksPtRatio);
  _tree->Branch("genEleSub_genTrks_DXY",&_genEleSub_genTrksDXY);
  _tree->Branch("genEleSub_genTrks_DZ",&_genEleSub_genTrksDZ);
  _tree->Branch("genEleSub_genTrks_PxlHits",&_genEleSub_genTrksPxlHits);
  _tree->Branch("genEleSub_genTrks_StrHits",&_genEleSub_genTrksStrHits);
  _tree->Branch("genEleSub_genPart_MinDR",&_genEleSub_genPartMinDR,"genEleSub_genPart_MinDR/I");

  return;
}

/*******************************************************************************
 * 
 */
void LowPtEleNtuplizer::endJob()
{}

/*******************************************************************************
 * 
 */
void LowPtEleNtuplizer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{}

/*******************************************************************************
 * 
 */
void LowPtEleNtuplizer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{}

/*******************************************************************************
 * 
 */
void LowPtEleNtuplizer::init() {

  _run = -1;
  _lumi = -1;
  _event = -1;
  _pileup = -1;
  
  _genLabDispl = -1.;
  _genProperDispl = -1.;
  
  _genParticlesN = 0;

  _genElesN = 0;
  _genElesPt.clear();
  _genElesEta.clear();
  _genElesPhi.clear();
  
  _genMuonsN = 0;
  _genMuonsPt.clear();
  _genMuonsEta.clear();
  _genMuonsPhi.clear();
  
  _genPionsN = 0;
  _genPionsPt.clear();
  _genPionsEta.clear();
  _genPionsPhi.clear();
  
  _genKaonsN = 0;
  _genKaonsPt.clear();
  _genKaonsEta.clear();
  _genKaonsPhi.clear();
  
  _genTrksN = 0;
  _genTrksPt.clear();
  _genTrksEta.clear();
  _genTrksPhi.clear();
  
  _gsfTrksN = 0;
  _gsfTrksPt.clear();
  _gsfTrksEta.clear();
  _gsfTrksPhi.clear();
  
  _gsfElesN = 0;
  _gsfElesPt.clear();
  _gsfElesEta.clear();
  _gsfElesPhi.clear();
  _gsfEles_ctfTrkPt.clear();
  _gsfEles_ctfTrkEta.clear();
  _gsfEles_ctfTrkPhi.clear();
  _gsfEles_gsfTrkPt.clear();
  _gsfEles_gsfTrkEta.clear();
  _gsfEles_gsfTrkPhi.clear();
  _gsfEles_trackerDriven.clear();
  _gsfEles_ecalDriven.clear();
  
  _recoMuonsN = 0;
  _recoMuonsPt.clear();
  _recoMuonsEta.clear();
  _recoMuonsPhi.clear();
  _recoMuons_global.clear();
  _recoMuons_tracker.clear();
  _recoMuons_pf.clear();
  
  _genEles_genTrksIdx.clear();
  _genEles_genTrksDR.clear();
  _genEles_genTrksDXY.clear();
  _genEles_genTrksDZ.clear();
  _genEles_gsfTrksIdx.clear();
  _genEles_gsfTrksDR.clear();
  _genEles_gsfTrksDXY.clear();
  _genEles_gsfTrksDZ.clear();
  _genEles_gsfElesIdx.clear();
  _genEles_gsfElesDR.clear();
  _genEles_gsfElesDXY.clear();
  _genEles_gsfElesDZ.clear();
  
  _genMuons_genTrksIdx.clear();
  _genMuons_genTrksDR.clear();
  _genMuons_genTrksDXY.clear();
  _genMuons_genTrksDZ.clear();
  _genMuons_recoMuonsIdx.clear();
  _genMuons_recoMuonsDR.clear();
  _genMuons_recoMuonsDXY.clear();
  _genMuons_recoMuonsDZ.clear();
  
  _genKaons_genTrksIdx.clear();
  _genKaons_genTrksDR.clear();
  _genKaons_genTrksDXY.clear();
  _genKaons_genTrksDZ.clear();
  
  _genEleLead_genTrksN = 0;
  _genEleLead_genTrksIdx.clear();
  _genEleLead_genTrksDR.clear();
  _genEleLead_genTrksDEta.clear();
  _genEleLead_genTrksDPhi.clear();
  _genEleLead_genTrksPtRatio.clear();
  _genEleLead_genTrksDXY.clear();
  _genEleLead_genTrksDZ.clear();
  _genEleLead_genTrksPxlHits.clear();
  _genEleLead_genTrksStrHits.clear();
  _genEleLead_genPartMinDR = -1.    ;
  
  _genEleSub_genTrksN = 0;
  _genEleSub_genTrksIdx.clear();
  _genEleSub_genTrksDR.clear();
  _genEleSub_genTrksDEta.clear();
  _genEleSub_genTrksDPhi.clear();
  _genEleSub_genTrksPtRatio.clear();
  _genEleSub_genTrksDXY.clear();
  _genEleSub_genTrksDZ.clear();
  _genEleSub_genTrksPxlHits.clear();
  _genEleSub_genTrksStrHits.clear();
  _genEleSub_genPartMinDR = -1.;
  
}

/*******************************************************************************
 * 
 */
float LowPtEleNtuplizer::deltaEta( float eta1, float eta2 ) {
  return abs(eta1 - eta2);
}

/*******************************************************************************
 * 
 */
float LowPtEleNtuplizer::deltaPhi( float phi1, float phi2 ) {
  float dphi = phi1 - phi2;
  while ( dphi >  1.*M_PI ) { dphi -= 2.*M_PI; }
  while ( dphi < -1.*M_PI ) { dphi += 2.*M_PI; }
  return abs(dphi);
}

/*******************************************************************************
 * 
 */
float LowPtEleNtuplizer::deltaR( float eta1, float eta2, float phi1, float phi2 ) {
  float deta = deltaEta(eta1,eta2);
  float dphi = deltaPhi(phi1,phi2);
  return sqrt(deta*deta + dphi*dphi);
}

/*******************************************************************************
 * 
 */
void LowPtEleNtuplizer::analyze(const edm::Event& iEvent, 
				const edm::EventSetup& eSetup)
{

  init();

  _run = iEvent.id().run();
  _lumi = iEvent.luminosityBlock();
  _event = iEvent.id().event();

  std::cout << "EventNumber: " << _event << std::endl;

  //////////
  // get gen event

//  edm::Handle<edm::HepMCProduct> hepMCProductHandle;
//  try { iEvent.getByToken(_hepMCProductTag, hepMCProductHandle); }
//  catch (...) { std::cout << "Problem with hepMCProductHandle" << std::endl; }

//  const HepMC::GenEvent* evt = 0;
//  if(hepMCProductHandle.isValid()) {
//    evt = hepMCProductHandle.product()->GetEvent();
//  } else { std::cout << "Problem with hepMCProductHandle.isValid()" << std::endl; }

  //////////
  // get pileup

  edm::Handle< std::vector<PileupSummaryInfo> > pileupHandle;
  try { iEvent.getByToken(_pileupTag, pileupHandle); }
  catch (...) { std::cout << "Problem with pileupHandle" << std::endl; }
  if ( pileupHandle.isValid() ) {
    for ( const auto& pu : *pileupHandle.product() ) {
      if ( pu.getBunchCrossing() == 0 ) { 
	_pileup = pu.getPU_NumInteractions();
	break;
      }
    }
  } else { std::cout << "Problem with pileupHandle.isValid()" << std:: endl; }

  //////////
  // get vertices

  edm::Handle<reco::VertexCollection> verticesHandle;
  try { iEvent.getByToken(_verticesTag, verticesHandle); }
  catch (...) { std::cout << "Problem with verticesHandle" << std::endl; }
  if ( verticesHandle.isValid() ) {
    _vertices = verticesHandle->size();
  } else { std::cout << "Problem with verticesHandle.isValid()" << std:: endl; }

  //////////
  // get gen particles

  edm::Handle< std::vector<reco::GenParticle> > genParticlesHandle;
  try { iEvent.getByToken(_genParticlesTag, genParticlesHandle); } 
  catch (...) { std::cout << "Problem with genParticlesHandle" << std::endl; }
  //std::sort(genParticlesHandle->begin(),genParticlesHandle->end(),LowPtEleNtuplizer::genSortByPt);

  // extract only interesting gen particles (from B decays)
  std::vector<const reco::GenParticle*> genParticles;
  if(genParticlesHandle.isValid()) {
    for ( const auto& gen : *genParticlesHandle.product() ) { 
      if ( gen.pt() > 0.1 && abs(gen.eta()) < 2.4 && // acceptance (pT > 0.1 GeV, |eta| < 2.4)
	   ( abs(gen.pdgId()) == 11 || // electron
	     abs(gen.pdgId()) == 13 || // muon
	     abs(gen.pdgId()) == 211 || // charged pion
	     abs(gen.pdgId()) == 321 ) && // charged kaon
	   gen.numberOfMothers() > 0 && abs(gen.mother()->pdgId()) == 521 ) { // B+/- hadron 
	genParticles.push_back( &gen );
	if ( abs(gen.pdgId()) == 11 ) { 
	  //@@ horrible! assumes 1 entry per B->K(*)ee event
	  // http://www.phys.ufl.edu/~avery/course/4390/f2015/lectures/relativistic_kinematics_1.pdf
	  _genLabDispl = 0.01 * sqrt( std::pow(gen.vx()-gen.mother()->vx(),2.) +
				      std::pow(gen.vy()-gen.mother()->vy(),2.) +
				      std::pow(gen.vz()-gen.mother()->vz(),2.) ); // cm -> m
	  _genProperDispl = _genLabDispl * (gen.mother()->mass()/gen.mother()->p());  // Lorentz factor
	}
      }
    }
  } else { std::cout << "Problem with genParticlesHandle.isValid()" << std::endl; }
  std::sort(genParticles.begin(),genParticles.end(),LowPtEleNtuplizer::genSortByPt);
  _genParticlesN = genParticles.size();

  // sort gen particles into electron, muons, pions, kaons 
  std::vector<const reco::GenParticle*> genElectrons;
  std::vector<const reco::GenParticle*> genMuons;
  std::vector<const reco::GenParticle*> genPions;
  std::vector<const reco::GenParticle*> genKaons;
  std::vector<const reco::GenParticle*> genHadrons;
//  std::vector<int> genElectrons;
//  std::vector<int> genMuons;
//  std::vector<int> genPions;
//  std::vector<int> genKaons;
//  std::vector<int> genHadrons;
  int gen_idx = 0;
  for ( const auto& gen: genParticles ) { 
    if ( abs(gen->pdgId()) == 11 ) { 
      genElectrons.push_back(gen);//gen_idx);
      _genElesPt.push_back(gen->pt());
      _genElesEta.push_back(gen->eta());
      _genElesPhi.push_back(gen->phi());
    } else if ( abs(gen->pdgId()) == 13 ) { 
      genMuons.push_back(gen);//gen_idx);
      _genMuonsPt.push_back(gen->pt());
      _genMuonsEta.push_back(gen->eta());
      _genMuonsPhi.push_back(gen->phi());
    } else if ( abs(gen->pdgId()) == 211 ) {
      genPions.push_back(gen);//gen_idx);
      genHadrons.push_back(gen);//gen_idx);
      _genPionsPt.push_back(gen->pt());
      _genPionsEta.push_back(gen->eta());
      _genPionsPhi.push_back(gen->phi());
    } else if ( abs(gen->pdgId()) == 321 ) {
      genKaons.push_back(gen);//gen_idx);
      genHadrons.push_back(gen);//gen_idx);
      _genKaonsPt.push_back(gen->pt());
      _genKaonsEta.push_back(gen->eta());
      _genKaonsPhi.push_back(gen->phi());
    }
    ++gen_idx;
  }
  _genElesN  = _genElesPt.size();
  _genMuonsN = _genMuonsPt.size();
  _genPionsN = _genPionsPt.size();
  _genKaonsN = _genKaonsPt.size();

  //////////
  // get tracks 

  edm::Handle< std::vector<reco::Track> > generalTracksHandle;
  try { iEvent.getByToken(_generalTracksTag, generalTracksHandle); } 
  catch (...) { std::cout << "Problem with generalTracksHandle" << std::endl; }

  std::vector<const reco::Track*> allTracks;
  if(generalTracksHandle.isValid()) {
    for ( const auto& trk : *generalTracksHandle.product() ) { 
      if ( trk.quality(reco::TrackBase::qualityByName("highPurity")) ) {
	allTracks.push_back( &trk );
      }
    }
  } else { std::cout << "Problem with generalTracksHandle.isValid()" << std::endl; }
  std::sort(allTracks.begin(),allTracks.end(),LowPtEleNtuplizer::trackSortByPt);

  for ( const auto& trk: allTracks ) { 
    _genTrksPt.push_back(trk->pt());
    _genTrksEta.push_back(trk->eta());
    _genTrksPhi.push_back(trk->phi());
  }
  _genTrksN = _genTrksPt.size();
  
  //////////
  // get electron seeds

//  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksHandle;
//  try { iEvent.getByToken(_gsfTracksTag, gsfTracksHandle); } 
//  catch (...) { std::cout << "Problem with gsfTracksHandle" << std::endl; }
//
//  std::vector<const reco::GsfTrack*> gsfTracks;
//  if(gsfTracksHandle.isValid()) {
//    for ( const auto& gsf : *gsfTracksHandle.product() ) { 
//      gsfTracks.push_back( &gsf );
//    }
//  } else { std::cout << "Problem with gsfTracksHandle.isValid()" << std::endl; }
//  std::sort(gsfTracks.begin(),gsfTracks.end(),LowPtEleNtuplizer::trackSortByPt);
//
//  for ( const auto& gsf: gsfTracks ) { 
//    _gsfTrksPt.push_back(gsf->pt());
//    _gsfTrksEta.push_back(gsf->eta());
//    _gsfTrksPhi.push_back(gsf->phi());
//  }
//  _gsfTrksN = _gsfTrksPt.size();

  //////////
  // get gsf tracks 

  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksHandle;
  try { iEvent.getByToken(_gsfTracksTag, gsfTracksHandle); } 
  catch (...) { std::cout << "Problem with gsfTracksHandle" << std::endl; }

  std::vector<const reco::GsfTrack*> gsfTracks;
  if(gsfTracksHandle.isValid()) {
    for ( const auto& gsf : *gsfTracksHandle.product() ) { 
      gsfTracks.push_back( &gsf );
    }
  } else { std::cout << "Problem with gsfTracksHandle.isValid()" << std::endl; }
  std::sort(gsfTracks.begin(),gsfTracks.end(),LowPtEleNtuplizer::trackSortByPt);

  for ( const auto& gsf: gsfTracks ) { 
    _gsfTrksPt.push_back(gsf->pt());
    _gsfTrksEta.push_back(gsf->eta());
    _gsfTrksPhi.push_back(gsf->phi());
  }
  _gsfTrksN = _gsfTrksPt.size();

  //////////
  // get gsf electrons 

  edm::Handle< std::vector<reco::GsfElectron> > gsfElectronsHandle;
  try { iEvent.getByToken(_gsfElectronsTag, gsfElectronsHandle); } 
  catch (...) { std::cout << "Problem with gsfElectronsHandle" << std::endl; }

  std::vector<const reco::GsfElectron*> gsfElectrons;
  if(gsfElectronsHandle.isValid()) {
    for ( const auto& ele : *gsfElectronsHandle.product() ) { 
      gsfElectrons.push_back( &ele );
    }
  } else { std::cout << "Problem with gsfElectronsHandle.isValid()" << std::endl; }
  std::sort(gsfElectrons.begin(),gsfElectrons.end(),LowPtEleNtuplizer::recoSortByPt);

  for ( const auto& ele: gsfElectrons ) { 
    _gsfElesPt.push_back(ele->pt());
    _gsfElesEta.push_back(ele->eta());
    _gsfElesPhi.push_back(ele->phi());
    _gsfEles_ctfTrkPt.push_back(ele->closestTrack().isNonnull() ? ele->closestTrack()->pt() : -999.);
    _gsfEles_ctfTrkEta.push_back(ele->closestTrack().isNonnull() ? ele->closestTrack()->eta() : -999.);
    _gsfEles_ctfTrkPhi.push_back(ele->closestTrack().isNonnull() ? ele->closestTrack()->phi() : -999.);
    _gsfEles_gsfTrkPt.push_back(ele->gsfTrack().isNonnull() ? ele->gsfTrack()->pt() : -999.);
    _gsfEles_gsfTrkEta.push_back(ele->gsfTrack().isNonnull() ? ele->gsfTrack()->eta() : -999.);
    _gsfEles_gsfTrkPhi.push_back(ele->gsfTrack().isNonnull() ? ele->gsfTrack()->phi() : -999.);
    _gsfEles_trackerDriven.push_back(ele->trackerDrivenSeed() ? 1 : 0);
    _gsfEles_ecalDriven.push_back(ele->ecalDrivenSeed() ? 1 : 0);
  }
  _gsfElesN = _gsfElesPt.size();

  //////////
  // get global muons 

  edm::Handle< std::vector<reco::Muon> > recoMuonsHandle;
  try { iEvent.getByToken(_recoMuonsTag, recoMuonsHandle); } 
  catch (...) { std::cout << "Problem with recoMuonsHandle" << std::endl; }

  std::vector<const reco::Muon*> recoMuons;
  if(recoMuonsHandle.isValid()) {
    for ( const auto& muon : *recoMuonsHandle.product() ) { 
      recoMuons.push_back( &muon );
    }
  } else { std::cout << "Problem with recoMuonsHandle.isValid()" << std::endl; }
  std::sort(recoMuons.begin(),recoMuons.end(),LowPtEleNtuplizer::recoSortByPt);
  
  for ( const auto& muon: recoMuons ) { 
    _recoMuonsPt.push_back(muon->pt());
    _recoMuonsEta.push_back(muon->eta());
    _recoMuonsPhi.push_back(muon->phi());
    _recoMuons_global.push_back(muon->isGlobalMuon() ? 1 : 0);
    _recoMuons_tracker.push_back(muon->isTrackerMuon() ? 1 : 0);
    _recoMuons_pf.push_back(muon->isPFMuon() ? 1 : 0);
  }
  _recoMuonsN = _recoMuonsPt.size();
  
  //////////
  // matching to gen electrons 

  std::vector<int> trk_matched_to_genele;
  std::vector<int> gsf_matched_to_genele;
  std::vector<int> ele_matched_to_genele;
  std::vector<int> trk_matched_to_lead_genele;
  int iele_gen = 0;
  for ( const auto& iter: genElectrons ) { 
    const reco::GenParticle* gen = iter;//genParticles[iter];
    
    // general tracks 
    float trk_dr = 1.e6;
    int trk_idx = -1;
    int itrk = 0;
    for ( const auto& trk: allTracks ) { 
      float dr = deltaR(gen->eta(),trk->eta(),gen->phi(),trk->phi());
      if ( dr < trk_dr ) {
	trk_dr = dr;
	trk_idx = itrk;
      }
      ++itrk;
    }
    if ( trk_idx >= 0 && std::find(trk_matched_to_genele.begin(),
				   trk_matched_to_genele.end(),
				   trk_idx) == trk_matched_to_genele.end() ) {
      trk_matched_to_genele.push_back(trk_idx);
      _genEles_genTrksIdx.push_back(trk_idx);
      _genEles_genTrksDR.push_back(trk_dr);
      float dx = allTracks[trk_idx]->vx() - gen->vx();
      float dy = allTracks[trk_idx]->vy() - gen->vy();
      float dz = allTracks[trk_idx]->vz() - gen->vz();
      _genEles_genTrksDXY.push_back( sqrt(dx*dx + dy*dy) );
      _genEles_genTrksDZ.push_back( fabs(dz) );
    } else {
      _genEles_genTrksIdx.push_back(-1);
      _genEles_genTrksDR.push_back(1.e6);
      _genEles_genTrksDXY.push_back(1.e6);
      _genEles_genTrksDZ.push_back(1.e6);
    }

    // gsf tracks 
    float gsf_dr = 1.e6;
    int gsf_idx = -1;
    int igsf = 0;
    for ( const auto& gsf: gsfTracks ) { 
      float dr = deltaR(gen->eta(),gsf->eta(),gen->phi(),gsf->phi());
      if ( dr < gsf_dr ) {
	gsf_dr = dr;
	gsf_idx = igsf;
      }
      ++igsf;
    }
    if ( gsf_idx >= 0 && std::find(gsf_matched_to_genele.begin(),
				   gsf_matched_to_genele.end(),
				   gsf_idx) == gsf_matched_to_genele.end() ) {
      gsf_matched_to_genele.push_back(gsf_idx);
      _genEles_gsfTrksIdx.push_back(gsf_idx);
      _genEles_gsfTrksDR.push_back(gsf_dr);
      float dx = gsfTracks[gsf_idx]->vx() - gen->vx();
      float dy = gsfTracks[gsf_idx]->vy() - gen->vy();
      float dz = gsfTracks[gsf_idx]->vz() - gen->vz();
      _genEles_gsfTrksDXY.push_back( sqrt(dx*dx + dy*dy) );
      _genEles_gsfTrksDZ.push_back( fabs(dz) );
    } else {
      _genEles_gsfTrksIdx.push_back(-1);
      _genEles_gsfTrksDR.push_back(1.e6);
      _genEles_gsfTrksDXY.push_back(1.e6);
      _genEles_gsfTrksDZ.push_back(1.e6);
    }

    // gsf electrons 
    float ele_dr = 1.e6;
    int ele_idx = -1;
    int iele = 0;
    for ( const auto& ele: gsfElectrons ) { 
      float dr = deltaR(gen->eta(),ele->eta(),gen->phi(),ele->phi());
      if ( dr < ele_dr ) {
	ele_dr = dr;
	ele_idx = iele;
      }
      ++iele;
    }
    if ( ele_idx >= 0 && std::find(ele_matched_to_genele.begin(),
				   ele_matched_to_genele.end(),
				   ele_idx) == ele_matched_to_genele.end() ) {
      ele_matched_to_genele.push_back(ele_idx);
      _genEles_gsfElesIdx.push_back(ele_idx);
      _genEles_gsfElesDR.push_back(ele_dr);
      float dx = gsfElectrons[ele_idx]->vx() - gen->vx();
      float dy = gsfElectrons[ele_idx]->vy() - gen->vy();
      float dz = gsfElectrons[ele_idx]->vz() - gen->vz();
      _genEles_gsfElesDXY.push_back( sqrt(dx*dx + dy*dy) );
      _genEles_gsfElesDZ.push_back( fabs(dz) );
    } else {
      _genEles_gsfElesIdx.push_back(-1);
      _genEles_gsfElesDR.push_back(1.e6);
      _genEles_gsfElesDXY.push_back(1.e6);
      _genEles_gsfElesDZ.push_back(1.e6);
    }
    
    // match tracks to lead gen electron
    if ( iele_gen == 0 ) {
      float max_dr = 0.3;
      int itrk = 0;
      for ( const auto& trk: allTracks ) { 
	float deta = deltaEta(gen->eta(),trk->eta());
	float dphi = deltaPhi(gen->phi(),trk->phi());
	float dr = deltaR(gen->eta(),trk->eta(),gen->phi(),trk->phi());
	if ( dr < max_dr ) {
	  _genEleLead_genTrksIdx.push_back(itrk);
	  _genEleLead_genTrksDR.push_back(dr);
	  _genEleLead_genTrksDEta.push_back(deta);
	  _genEleLead_genTrksDPhi.push_back(dphi);
	  _genEleLead_genTrksPtRatio.push_back(gen->pt() > 0. ? trk->pt() / gen->pt() : 0.);
	  float dx = trk->vx() - gen->vx();
	  float dy = trk->vy() - gen->vy();
	  float dz = trk->vz() - gen->vz();
	  _genEleLead_genTrksDXY.push_back( sqrt(dx*dx + dy*dy) );
	  _genEleLead_genTrksDZ.push_back( fabs(dz) );
	  _genEleLead_genTrksPxlHits.push_back( trk->hitPattern().numberOfValidPixelHits() );
	  _genEleLead_genTrksPxlHits.push_back( trk->hitPattern().numberOfValidStripHits() );
	}
	++itrk;
      }
      _genEleLead_genTrksN = _genEleLead_genTrksIdx.size();
      // calc dr for nearest gen particle 
      float gen_dr = 1.e6;
      int iother_gen = 0;
      for ( const auto& iter: genElectrons ) { 
	const reco::GenParticle* other = iter;//genParticles[iter];
	if ( iother_gen != iele_gen ) {
	  float dr = deltaR(gen->eta(),other->eta(),gen->phi(),other->phi());
	  if ( dr < gen_dr ) {
	    gen_dr = dr;
	  }
	}
      }
      for ( const auto& iter: genHadrons ) { 
	const reco::GenParticle* other = iter;//genParticles[iter];
	float dr = deltaR(gen->eta(),other->eta(),gen->phi(),other->phi());
	if ( dr < gen_dr ) {
	  gen_dr = dr;
	}
      }
      _genEleLead_genPartMinDR = gen_dr;
    } // lead gen electron

    // match tracks to subleading gen electron
    if ( iele_gen == 1 ) {
      // match tracks to to subleading gen electron 
      float max_dr = 0.3;
      int itrk = 0;
      for ( const auto& trk: allTracks ) { 
	float deta = deltaEta(gen->eta(),trk->eta());
	float dphi = deltaPhi(gen->phi(),trk->phi());
	float dr = deltaR(gen->eta(),trk->eta(),gen->phi(),trk->phi());
	if ( dr < max_dr ) {
	  _genEleSub_genTrksIdx.push_back(itrk);
	  _genEleSub_genTrksDR.push_back(dr);
	  _genEleSub_genTrksDEta.push_back(deta);
	  _genEleSub_genTrksDPhi.push_back(dphi);
	  _genEleSub_genTrksPtRatio.push_back(gen->pt() > 0. ? trk->pt() / gen->pt() : 0.);
	  float dx = trk->vx() - gen->vx();
	  float dy = trk->vy() - gen->vy();
	  float dz = trk->vz() - gen->vz();
	  _genEleSub_genTrksDXY.push_back( sqrt(dx*dx + dy*dy) );
	  _genEleSub_genTrksDZ.push_back( fabs(dz) );
	  _genEleSub_genTrksPxlHits.push_back( trk->hitPattern().numberOfValidPixelHits() );
	  _genEleSub_genTrksPxlHits.push_back( trk->hitPattern().numberOfValidStripHits() );
	}
	++itrk;
      }
      _genEleSub_genTrksN = _genEleSub_genTrksIdx.size();
      // calc dr for nearest gen particle 
      float gen_dr = 1.e6;
      int iother_gen = 0;
      for ( const auto& iter: genElectrons ) { 
	const reco::GenParticle* other = iter;//genParticles[iter];
	if ( iother_gen != iele_gen ) {
	  float dr = deltaR(gen->eta(),other->eta(),gen->phi(),other->phi());
	  if ( dr < gen_dr ) {
	    gen_dr = dr;
	  }
	}
      }
      for ( const auto& iter: genHadrons ) { 
	const reco::GenParticle* other = iter;//genParticles[iter];
	float dr = deltaR(gen->eta(),other->eta(),gen->phi(),other->phi());
	if ( dr < gen_dr ) {
	  gen_dr = dr;
	}
      }
      _genEleSub_genPartMinDR = gen_dr;
    } // subleading gen electron

    ++iele_gen; // increment counter
  }

  //////////
  // matching to gen muons 

  std::vector<int> trk_matched_to_genmuon;
  std::vector<int> muon_matched_to_genmuon;
  for ( const auto& iter: genMuons ) { 
    const reco::GenParticle* gen = iter;//genParticles[iter];
    
    // general tracks 
    float trk_dr = 1.e6;
    int trk_idx = -1;
    int itrk = 0;
    for ( const auto& trk: allTracks ) { 
      float dr = deltaR(gen->eta(),trk->eta(),gen->phi(),trk->phi());
      if ( dr < trk_dr ) {
	trk_dr = dr;
	trk_idx = itrk;
      }
      ++itrk;
    }
    if ( trk_idx >= 0 && std::find(trk_matched_to_genmuon.begin(),
				   trk_matched_to_genmuon.end(),
				   trk_idx) == trk_matched_to_genmuon.end() ) {
      trk_matched_to_genmuon.push_back(trk_idx);
      _genMuons_genTrksIdx.push_back(trk_idx);
      _genMuons_genTrksDR.push_back(trk_dr);
      float dx = allTracks[trk_idx]->vx() - gen->vx();
      float dy = allTracks[trk_idx]->vy() - gen->vy();
      float dz = allTracks[trk_idx]->vz() - gen->vz();
      _genMuons_genTrksDXY.push_back( sqrt(dx*dx + dy*dy) );
      _genMuons_genTrksDZ.push_back( fabs(dz) );
    } else {
      _genMuons_genTrksIdx.push_back(-1);
      _genMuons_genTrksDR.push_back(1.e6);
      _genMuons_genTrksDXY.push_back(1.e6);
      _genMuons_genTrksDZ.push_back(1.e6);
    }

    // reco muons 
    float muon_dr = 1.e6;
    int muon_idx = -1;
    int imuon = 0;
    for ( const auto& muon: recoMuons ) { 
      float dr = deltaR(gen->eta(),muon->eta(),gen->phi(),muon->phi());
      if ( dr < muon_dr ) {
	muon_dr = dr;
	muon_idx = imuon;
      }
      ++imuon;
    }
    if ( muon_idx >= 0 && std::find(muon_matched_to_genmuon.begin(),
				    muon_matched_to_genmuon.end(),
				    muon_idx) == muon_matched_to_genmuon.end() ) {
      muon_matched_to_genmuon.push_back(muon_idx);
      _genMuons_recoMuonsIdx.push_back(muon_idx);
      _genMuons_recoMuonsDR.push_back(muon_dr);
      float dx = recoMuons[muon_idx]->vx() - gen->vx();
      float dy = recoMuons[muon_idx]->vy() - gen->vy();
      float dz = recoMuons[muon_idx]->vz() - gen->vz();
      _genMuons_recoMuonsDXY.push_back( sqrt(dx*dx + dy*dy) );
      _genMuons_recoMuonsDZ.push_back( fabs(dz) );
    } else {
      _genMuons_recoMuonsIdx.push_back(-1);
      _genMuons_recoMuonsDR.push_back(1.e6);
      _genMuons_recoMuonsDXY.push_back(1.e6);
      _genMuons_recoMuonsDZ.push_back(1.e6);
    }
    
  }


  //////////
  // matching to gen kaons

  std::vector<int> trk_matched_to_genkaon;
  for ( const auto& iter: genKaons ) {
    const reco::GenParticle* gen = iter;//genParticles[iter];

    // general tracks
    float trk_dr = 1.e6;
    int trk_idx = -1;
    int itrk = 0;
    for ( const auto& trk: allTracks ) {
      float dr = deltaR(gen->eta(),trk->eta(),gen->phi(),trk->phi());
      if ( dr < trk_dr ) {
	trk_dr = dr;
	trk_idx = itrk;
      }
      ++itrk;
    }
    if ( trk_idx >= 0 && std::find(trk_matched_to_genkaon.begin(),
				   trk_matched_to_genkaon.end(),
				   trk_idx) == trk_matched_to_genkaon.end() ) {
      trk_matched_to_genkaon.push_back(trk_idx);
      _genKaons_genTrksIdx.push_back(trk_idx);
      _genKaons_genTrksDR.push_back(trk_dr);
      float dx = allTracks[trk_idx]->vx() - gen->vx();
      float dy = allTracks[trk_idx]->vy() - gen->vy();
      float dz = allTracks[trk_idx]->vz() - gen->vz();
      _genKaons_genTrksDXY.push_back( sqrt(dx*dx + dy*dy) );
      _genKaons_genTrksDZ.push_back( fabs(dz) );
    } else {
      _genKaons_genTrksIdx.push_back(-1);
      _genKaons_genTrksDR.push_back(1.e6);
      _genKaons_genTrksDXY.push_back(1.e6);
      _genKaons_genTrksDZ.push_back(1.e6);
    }

  }

  // populate branches
  _tree->Fill();
  
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(LowPtEleNtuplizer);

#endif //LOWPTELENTUPLIZER_H
