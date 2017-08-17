// -*- C++ -*-
//
// Package:    METSigTuning/MakeNtuple
// Class:      MakeNtuple
// 
/**\class MakeNtuple MakeNtuple.cc METSigTuning/MakeNtuple/src/MakeNtuple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nathan Mirman
//         Created:  Wed, 26 Nov 2014 14:22:08 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
//#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"


#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"

#include "Math/LorentzVector.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>

#define NUMPUBINS 75

//
// class declaration
//

class MakeEleNtuple : public edm::EDAnalyzer {
   public:
      explicit MakeEleNtuple(const edm::ParameterSet&);
      ~MakeEleNtuple();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      std::vector<pat::Jet> cleanJets(double, double,
            std::vector<pat::Jet>&, std::vector<reco::Candidate::LorentzVector>&);
      edm::EDGetTokenT<edm::View<reco::Candidate> > inputToken_;
      edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
      std::vector< edm::EDGetTokenT<edm::View<reco::Candidate> > > lepTokens_;
      edm::EDGetTokenT<edm::View<pat::MET> > metToken_;
      edm::EDGetTokenT<GenEventInfoProduct> genToken_;
      edm::EDGetTokenT<LHEEventProduct> lheToken_;
      edm::EDGetTokenT< std::vector<pat::Muon> > muonToken_;
      edm::EDGetTokenT< std::vector<pat::Electron> > electronToken_;
      edm::EDGetTokenT<edm::View<reco::Vertex> > verticesToken_;

      // electron id
      edm::EDGetTokenT<edm::ValueMap<bool>> electronVetoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool>> electronLooseIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool>> electronMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool>> electronTightIdMapToken_;

      Bool_t runOnMC_;
      //edm::InputTag addPileupInfo_;
      edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
      edm::EDGetTokenT<double> rhoToken_;
      std::string jetSFType_;
      std::string jetResPtType_;
      std::string jetResPhiType_;
      

      edm::LumiReWeighting LumiWeights_;

      //std::string pfjetCorrectorL1_;
      //std::string pfjetCorrectorL123_;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

      static const double PU2015_MCf[NUMPUBINS];
      static const double PU2015_Dataf[NUMPUBINS];

      TTree *results_tree;
      TFile *OutFile__file;
      std::string OutputFileName_;

      Long64_t run, event, lumi;

      std::vector<double> lep_pt, lep_energy, lep_phi, lep_eta;
      std::vector<double> muon_pt, muon_energy, muon_phi, muon_eta;
      std::vector<int> muon_charge;
      std::vector<double> ele_pt, ele_energy, ele_phi, ele_eta, ele_loose, ele_medium, ele_tight;
      std::vector<int> ele_charge;
      std::vector<double> jet_pt, jet_energy, jet_phi, jet_eta;
      std::vector<double> cand_pt, cand_energy, cand_phi, cand_eta, cand_x, cand_y;
      std::vector<double> jet_sigmapt, jet_sigmaphi;
      std::vector<double> jet_corrL1, jet_corrL123;
      std::vector<bool> jet_passid;
      std::vector<double> jet_sf;
      double c_xx, c_xy, c_yy,x_tot, x_bar, y_tot, y_bar;
      double dimuon_mass;
      double met_pt, met_energy, met_phi, met_eta, met_sumpt, met_sig, alt_sumpt;

      int nvertices, met_sumpt_inputs, nmuons, nMediumElectrons, nLooseElectrons, nTightElectrons;
      double weight_pu;
      double mcweight, mcweightSum;

      double jetThreshold;

      int events_total, events_pass;

      double T_nvertices;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MakeEleNtuple::MakeEleNtuple(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   inputToken_ = consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("src"));
   jetToken_ = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"));
   std::vector<edm::InputTag> srcLeptonsTags = iConfig.getParameter< std::vector<edm::InputTag> >("leptons");
   for(std::vector<edm::InputTag>::const_iterator it=srcLeptonsTags.begin();it!=srcLeptonsTags.end();it++) {
      lepTokens_.push_back( consumes<edm::View<reco::Candidate> >( *it ) );
   }
   metToken_                        = consumes<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("met"));

   muonToken_           = consumes< std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
   electronToken_       = consumes< std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
   verticesToken_       = consumes< edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"));
   genToken_            = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"));
   lheToken_            = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheprod"));
   runOnMC_             = iConfig.getUntrackedParameter<Bool_t>("runOnMC");
   //addPileupInfo_ = iConfig.getUntrackedParameter<edm::InputTag>("pileup");
   pileupToken_         = consumes<edm::View<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>("pileup"));
   rhoToken_            = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));

   electronVetoIdMapToken_    = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"));
   electronLooseIdMapToken_   = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"));
   electronMediumIdMapToken_  = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("electronMediumIdMap"));
   electronTightIdMapToken_   = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("electronTightIdMap"));

   jetSFType_           = iConfig.getParameter<std::string>("srcJetSF");
   jetResPtType_        = iConfig.getParameter<std::string>("srcJetResPt");
   jetResPhiType_       = iConfig.getParameter<std::string>("srcJetResPhi");

   //pfjetCorrectorL1_  = iConfig.getUntrackedParameter<std::string>("pfjetCorrectorL1");
   //pfjetCorrectorL123_ = iConfig.getUntrackedParameter<std::string>("pfjetCorrectorL123");

   jetThreshold = 15;


   OutputFileName_ = "ntuple.root";

}


MakeEleNtuple::~MakeEleNtuple()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

const double MakeEleNtuple::PU2015_MCf[NUMPUBINS] = {
   // Updated to PU for Summer16 from: SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py 
   // pileup distribution for Spring2016 MC
   // obtained at https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVPileUpDescription#Startup2015
   // from file SimGeneral/MixingModule/python/mix_2015_25ns_FallMC_matchData_PoissonOOTPU_cfi.py
    1.78653e-05,
    2.56602e-05,
    5.27857e-05,
    8.88954e-05,
    0.000109362,
    0.000140973,
    0.000240998,
    0.00071209,
    0.00130121,
    0.00245255,
    0.00502589,
    0.00919534,
    0.0146697,
    0.0204126,
    0.0267586,
    0.0337697,
    0.0401478,
    0.0450159,
    0.0490577,
    0.0524855,
    0.0548159,
    0.0559937,
    0.0554468,
    0.0537687,
    0.0512055,
    0.0476713,
    0.0435312,
    0.0393107,
    0.0349812,
    0.0307413,
    0.0272425,
    0.0237115,
    0.0208329,
    0.0182459,
    0.0160712,
    0.0142498,
    0.012804,
    0.011571,
    0.010547,
    0.00959489,
    0.00891718,
    0.00829292,
    0.0076195,
    0.0069806,
    0.0062025,
    0.00546581,
    0.00484127,
    0.00407168,
    0.00337681,
    0.00269893,
    0.00212473,
    0.00160208,
    0.00117884,
    0.000859662,
    0.000569085,
    0.000365431,
    0.000243565,
    0.00015688,
    9.88128e-05,
    6.53783e-05,
    3.73924e-05,
    2.61382e-05,
    2.0307e-05,
    1.73032e-05,
    1.435e-05,
    1.36486e-05,
    1.35555e-05,
    1.37491e-05,
    1.34255e-05,
    1.33987e-05,
    1.34061e-05,
    1.34211e-05,
    1.34177e-05,
    1.32959e-05,
    1.33287e-05
};

const double MakeEleNtuple::PU2015_Dataf[NUMPUBINS] = {
   // obtained with pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 69200 --maxPileupBin 75 --numPileupBins 75 PU_2016_${LUMI}_XSecCentral.root
   // 'true' distribution for 2016 RunA-H dataset
   // obtained with pileupCalc.py (04.21.2016)
    238797.04313,
    837542.889595,
    2308427.00478,
    3124754.34372,
    4476191.30543,
    5995911.21323,
    7000895.56315,
    12891652.1791,
    35261733.5341,
    78701230.1513,
    176945810.65,
    360089516.856,
    602766485.926,
    876519381.282,
    1174474378.51,
    1489058683.57,
    1759351600.68,
    1943925871.14,
    2049172352.74,
    2101581720.96,
    2132787284.84,
    2149099182.3,
    2128985927.95,
    2062648995.45,
    1962883779.25,
    1841872058.94,
    1704135500.51,
    1554522565.28,
    1399489422.58,
    1243532663.99,
    1088821328.46,
    937304754.427,
    792044050.775,
    656717691.231,
    534466802.768,
    427126771.855,
    335105585.322,
    257724603.25,
    193751376.913,
    141830894.568,
    100671434.732,
    69013863.5534,
    45540080.2836,
    28847477.6909,
    17506316.1645,
    10162639.6709,
    5637781.04925,
    2987281.99949,
    1512002.26145,
    731845.412615,
    339821.984793,
    152545.357339,
    67404.8209077,
    30489.6912883,
    15152.1058363,
    8975.91109751,
    6496.15482279,
    5434.80517454,
    4889.95760724,
    4521.71627396,
    4208.46442542,
    3909.7628178,
    3614.27410537,
    3320.72226327,
    3031.09566114,
    2748.23668947,
    2474.97653154,
    2213.81721565,
    1966.81521602,
    1735.5463732,
    1521.10908823,
    1324.14905961,
    1144.89778819,
    983.220158682,
    838.667561122
};

// ------------ method called for each event  ------------
void
MakeEleNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // clear all vectors
   lep_pt.clear();
   lep_energy.clear();
   lep_phi.clear();
   lep_eta.clear();
   muon_pt.clear();
   muon_charge.clear();
   muon_energy.clear();
   muon_phi.clear();
   muon_eta.clear();

   ele_pt.clear();
   ele_charge.clear();
   ele_energy.clear();
   ele_phi.clear();
   ele_eta.clear();
   ele_loose.clear();
   ele_medium.clear();
   ele_tight.clear();
   
   jet_pt.clear();
   jet_passid.clear();
   jet_energy.clear();
   jet_phi.clear();
   jet_eta.clear();

   cand_pt.clear();
   cand_energy.clear();
   cand_phi.clear();
   cand_eta.clear();
   cand_x.clear();
   cand_y.clear();

   jet_sigmapt.clear();
   jet_sigmaphi.clear();
   jet_sf.clear();

   run = iEvent.id().run();
   event = iEvent.id().event();
   lumi = iEvent.id().luminosityBlock();

   Handle<View<reco::Candidate> > input;
   iEvent.getByToken(inputToken_, input);

   // offline primary vertices
   edm::Handle<edm::View<reco::Vertex> > vertices;
   iEvent.getByToken(verticesToken_, vertices);
   nvertices = int(vertices->size());

   // leptons
   std::vector<reco::CandidatePtr> footprint;
   std::vector<reco::Candidate::LorentzVector> leptons;
   for ( std::vector<EDGetTokenT<View<reco::Candidate> > >::const_iterator srcLeptons_i = lepTokens_.begin();
         srcLeptons_i != lepTokens_.end(); ++srcLeptons_i ) {
      Handle<reco::CandidateView> leptons_i;
      iEvent.getByToken(*srcLeptons_i, leptons_i);
      for ( reco::CandidateView::const_iterator lepton = leptons_i->begin();
            lepton != leptons_i->end(); ++lepton ) {
         // cut on lepton pt
         if( lepton->pt() > 10 ){
            leptons.push_back(lepton->p4());
            for( unsigned int n=0; n < lepton->numberOfSourceCandidatePtrs(); n++){
               if( lepton->sourceCandidatePtr(n).isNonnull() and lepton->sourceCandidatePtr(n).isAvailable() ){
                  footprint.push_back(lepton->sourceCandidatePtr(n));
               }
            }
         }
      }
   }

   edm::Handle<edm::ValueMap<bool>> veto_id_decisions;
   edm::Handle<edm::ValueMap<bool>> loose_id_decisions;
   edm::Handle<edm::ValueMap<bool>> medium_id_decisions;
   edm::Handle<edm::ValueMap<bool>> tight_id_decisions;
   iEvent.getByToken(electronVetoIdMapToken_, veto_id_decisions);
   iEvent.getByToken(electronLooseIdMapToken_, loose_id_decisions);
   iEvent.getByToken(electronMediumIdMapToken_, medium_id_decisions);
   iEvent.getByToken(electronTightIdMapToken_, tight_id_decisions);


   // electrons for event selection (not for tuning)!
   Handle< std::vector<pat::Electron> > electrons;
   iEvent.getByToken(electronToken_, electrons);
   nMediumElectrons = 0;
   nLooseElectrons = 0;
   nTightElectrons = 0;
   //std::cout << "Trying to get electrons" << std::endl;
   for ( std::vector<pat::Electron>::const_iterator electron = electrons->begin();
         electron != electrons->end(); ++electron ) {

      const edm::Ptr<pat::Electron> elPtr (electrons, electron - electrons->begin() );
      bool isLoose = (*loose_id_decisions) [elPtr];
      bool isMedium = (*medium_id_decisions) [elPtr];
      bool isTight = (*tight_id_decisions) [elPtr];

      //if (isLoose && true ) std::cout << "Loose" << electron->pt() << std::endl;
      //if (isMedium && true ) std::cout << "Medium" << electron->pt() << std::endl;
      //if (isTight && true ) std::cout << "Tight" << electron->pt() << std::endl;


      //bool medEleId = electron->isMediumElectron((*(vertices->begin())));
      //bool looseEleId = electron->isLooseElectron((*(vertices->begin())));

      double e_dr04chHad      = electron->dr04TkSumPt();
      double e_dr04photons    = electron->dr04EcalRecHitSumEt();
      double e_dr04neutHad    = electron->dr04HcalTowerSumEt();

      bool eleIso = (e_dr04chHad + e_dr04neutHad + e_dr04photons)/electron->pt() < 0.12;//can be changed to 0.15 before
      if ( electron->pt() > 20 and fabs(electron->eta()) < 2.4 and eleIso and isLoose){
            //std::cout << "Found an electron" << std::endl;    
            ele_pt.push_back( electron->pt() );
            ele_eta.push_back( electron->eta() );
            ele_phi.push_back( electron->phi() );
            ele_charge.push_back( electron->charge() );
            ele_loose.push_back( isLoose );
            ele_medium.push_back( isMedium );
            ele_tight.push_back( isTight );
            nLooseElectrons++;
            if (isMedium) nMediumElectrons++;
            if (isTight) nTightElectrons++;
            //if isMedium: 
      //  nMediumElectrons++;
      }
   }

   // muons (for event selection)
   Handle< std::vector<pat::Muon> > muons;
   iEvent.getByToken(muonToken_, muons);
   nmuons = 0;
   double charge = 1.0;
   for ( std::vector<pat::Muon>::const_iterator muon = muons->begin();
         muon != muons->end(); ++muon ) {

      bool muId = muon->isTightMuon((*(vertices->begin())));
      //bool muId = muon->isLooseMuon();

      double dr04chHad = muon->pfIsolationR04().sumChargedHadronPt;
      double dr04neutHad = muon->pfIsolationR04().sumNeutralHadronEt;
      double dr04photons = muon->pfIsolationR04().sumPhotonEt;

      bool muIso = (dr04chHad + dr04neutHad + dr04photons)/muon->pt() < 0.12;//can be changed to 0.15 before

      if( muon->pt() > 20 and fabs(muon->eta()) < 2.4 and muId and muIso ){
         muon_pt.push_back( muon->pt() );
         muon_energy.push_back( muon->energy() );
         muon_phi.push_back( muon->phi() );
         muon_eta.push_back( muon->eta() );
         muon_charge.push_back( muon->charge() );
         nmuons++;
         charge *= muon->charge();
      }
   }
   dimuon_mass = 0.;
   if( nmuons == 2 ){
      TLorentzVector mu1, mu2;
      mu1.SetPtEtaPhiE( muon_pt[0], muon_eta[0], muon_phi[0], muon_energy[0] );
      mu2.SetPtEtaPhiE( muon_pt[1], muon_eta[1], muon_phi[1], muon_energy[1] );
      dimuon_mass = (mu1+mu2).M();
   }

   // jet energy corrections
   //const JetCorrector* corrL1  = JetCorrector::getJetCorrector (pfjetCorrectorL1_, iSetup);
   //const JetCorrector* corrL123 = JetCorrector::getJetCorrector (pfjetCorrectorL123_, iSetup);
       
   // jets
   Handle<View<pat::Jet>> inputJets;
   iEvent.getByToken( jetToken_, inputJets );
   std::vector<pat::Jet> jets;
   for(View<pat::Jet>::const_iterator jet = inputJets->begin(); jet != inputJets->end(); ++jet) {
      jets.push_back( *jet );
   }

   // disambiguate jets and leptons
   std::vector<pat::Jet> cleanjets = cleanJets(jetThreshold, 0.4, jets, leptons);

   // loop over jets to disambiguate candidates
   for(std::vector<pat::Jet>::const_iterator jet = cleanjets.begin(); jet != cleanjets.end(); ++jet) {
      for( unsigned int n=0; n < jet->numberOfSourceCandidatePtrs(); n++){
         if( jet->sourceCandidatePtr(n).isNonnull() and jet->sourceCandidatePtr(n).isAvailable() ){
            footprint.push_back(jet->sourceCandidatePtr(n));
         }
      }
   }

   // met
   edm::Handle<edm::View<pat::MET> > metHandle;

   iEvent.getByToken(metToken_, metHandle);

   const pat::MET& met                      = (*metHandle)[0];


   // put met into a 4-vector, implement type-1 corrections
   /*
   reco::Candidate::LorentzVector pmetvec = met.shiftedP4_74x(pat::MET::METUncertainty(12),pat::MET::Raw);
   TLorentzVector metvec;
   metvec.SetPxPyPzE(pmetvec.px(), pmetvec.py(), pmetvec.pz(), pmetvec.energy());
   for(std::vector<pat::Jet>::const_iterator jet = cleanjets.begin(); jet != cleanjets.end(); ++jet){
      TLorentzVector jettemp;
      jettemp.SetPxPyPzE(jet->px(), jet->py(), jet->pz(), jet->energy());
      double cL1 = corrL1->correction( *jet, iEvent, iSetup );
      double cL123 = corrL123->correction( *jet, iEvent, iSetup );
      metvec -= (cL123-cL1)*jettemp;
   }
   */


   met_pt = met.pt();
   met_sig = met.metSignificance();
   met_energy = met.energy();
   met_phi = met.phi();
   met_eta = met.eta();


   /*
   std::cout << "MET (new) = " << met_pt << std::endl;

   edm::Handle<edm::View<reco::MET> > metHandle2;
   iEvent.getByLabel("slimmedMETs", metHandle2);
   const reco::MET& met2 = (*metHandle2)[0];

   double met2_pt = met2.pt();
   std::cout << "MET (old) = " << met2_pt << std::endl;

   edm::Handle<edm::View<reco::MET> > metHandle3;
   iEvent.getByLabel("pfMetRERUN", metHandle3);
   const reco::MET& met3 = (*metHandle3)[0];

   double met3_pt = met3.pt();
   std::cout << "MET (uncorr) = " << met3_pt << std::endl;

   std::cout << " ************************* " << std::endl;

   fflush(stdout);
   */

   // candidates
   std::vector<reco::Candidate::LorentzVector> candidates;
   for(View<reco::Candidate>::const_iterator cand = input->begin();
         cand != input->end(); ++cand) {
      unsigned int iter = cand - input->begin();
      if (std::find(footprint.begin(), footprint.end(),
               reco::CandidatePtr(input,iter)) != footprint.end()) {
         continue;
      }
      candidates.push_back( cand->p4() );
   }

   // resolutions
   /*
   std::string path = "METSigTuning/MakeNtuple/data/";
   std::string resEra_ = "Summer15_25nsV6";
   std::string resAlg_ = "AK4PFchs";
   std::string ptResFileName  = path + resEra_ + "_MC_PtResolution_" +resAlg_+".txt";
   std::string phiResFileName = path + resEra_ + "_MC_PhiResolution_" +resAlg_+".txt";
   std::string sfResFileName  = path + resEra_ + "_DATAMCSF_" +resAlg_+".txt";

   std::string phiResFileName_old = path+"Spring10_PhiResolution_AK5PF.txt";

   // old framework
   FileInPath fpt(ptResFileName);
   FileInPath fphi(phiResFileName);
   FileInPath fsf(sfResFileName);

   FileInPath fphi_old(phiResFileName_old);

   // accessJERs using text files
   //JME::JetResolution resolution_pt = JME::JetResolution(fpt.fullPath().c_str());
   //JME::JetResolutionScaleFactor resolution_sf = JME::JetResolutionScaleFactor(fsf.fullPath().c_str());
   //JME::JetResolution resolution_phi = JME::JetResolution(fphi.fullPath().c_str());
   */

   // access JERs using DB
   JME::JetResolution resolution_pt = JME::JetResolution::get(iSetup, jetResPtType_);
   JME::JetResolution resolution_phi = JME::JetResolution::get(iSetup, jetResPhiType_);
   JME::JetResolutionScaleFactor resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, jetSFType_);

   //JetResolution *phiRes_ = new JetResolution(fphi_old.fullPath().c_str(),false);

   //
   // begin ttree variables
   //

   if( runOnMC_ ){
      Handle<GenEventInfoProduct> gi;
      iEvent.getByToken(genToken_, gi);

      //Handle<LHEEventProduct> li;
      //iEvent.getByToken(lheToken_, li) ;

      mcweight = gi->weight();
      //mcweightSum = li->originalXWGTUP();
      //std::cout << mcweight << " " << mcweightSum << std::endl;
      //fflush(stdout);
   }

   weight_pu = 1.0;
   if( runOnMC_ ){
      View<PileupSummaryInfo>::const_iterator PVI;
      Handle<View<PileupSummaryInfo> > PupInfo;
      iEvent.getByToken(pileupToken_, PupInfo);
   //Handle<View<pat::Jet>> inputJets;
   //Event.getByToken( jetToken_, inputJets );
   //std::vector<pat::Jet> jets;
   //for(View<pat::Jet>::const_iterator jet = inputJets->begin(); jet != inputJets->end(); ++jet) {

      float Tnvtx = -1.0;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

         int BX = PVI->getBunchCrossing();

         if(BX == 0) { 
            Tnvtx = PVI->getTrueNumInteractions(); 
            continue;
         }

         weight_pu = LumiWeights_.weight( Tnvtx );
         T_nvertices = Tnvtx;
      }
   }

   // calculate met_sumpt
   // candidates already have clean jets removed
   met_sumpt_inputs = 0;
   met_sumpt = 0;
   for( std::vector<reco::Candidate::LorentzVector>::const_iterator cand = candidates.begin();
         cand != candidates.end(); ++cand){
      met_sumpt += cand->Pt();
      cand_pt.push_back( cand->Pt() );
      //cand_energy.push_back( cand->energy() );
      cand_phi.push_back( cand->Phi() );
      cand_eta.push_back( cand->Eta() );
      met_sumpt_inputs++;
   }

   /*
   // loop over leptons
   for ( std::vector<reco::Candidate::LorentzVector>::const_iterator lepton = leptons.begin();
         lepton != leptons.end(); ++lepton ) {
      lep_pt.push_back( lepton->Pt() );
      lep_energy.push_back( lepton->E() );
      lep_phi.push_back( lepton->Phi() );
      lep_eta.push_back( lepton->Eta() );
   }
   */
   edm::Handle<double> rho;
   iEvent.getByToken(rhoToken_, rho);

   // loop over jets
   for(std::vector<pat::Jet>::const_iterator jet = cleanjets.begin(); jet != cleanjets.end(); ++jet) {
      double jpt  = jet->pt();
      double jeta = jet->eta();

      // jet energy resolutions
      //double jeta_res = (fabs(jeta) < 9.9) ? jeta : 9.89; // JetResolutions defined for |eta|<9.9
      //TF1* fPtEta    = ptRes_ -> parameterEta("sigma",jeta_res);
      //TF1* fPhiEta   = phiRes_-> parameterEta("sigma",jeta_res);
      //double sigmapt = fPtEta->Eval(jpt);
      //double sigmaphi = fPhiEta->Eval(jpt);
      //delete fPtEta;
      //delete fPhiEta;

      // new framework
      JME::JetParameters parameters;
      parameters.setJetPt(jpt).setJetEta(jeta).setRho(*rho);
      double sigmapt = resolution_pt.getResolution(parameters);
      double sigmaphi = resolution_phi.getResolution(parameters);
      double sf = resolution_sf.getScaleFactor(parameters);

      // split into high-pt and low-pt sector
      if( jpt > jetThreshold ){
         // high-pt jets enter into the covariance matrix via JER

         jet_pt.push_back( jet->pt() );
         jet_energy.push_back( jet->energy() );
         jet_phi.push_back( jet->phi() );
         jet_eta.push_back( jet->eta() );
         jet_sigmapt.push_back( sigmapt );
         jet_sigmaphi.push_back( sigmaphi );
         jet_sf.push_back( sf );

         // jet id
         bool id_24 = true;
         bool id_30 = true;
         bool id_forw = true;
         if( fabs(jet->eta()) <= 3.0 ){
            id_30 = jet->neutralHadronEnergyFraction() < 0.99 and jet->neutralEmEnergyFraction() < 0.99
               and (jet->chargedMultiplicity() + jet->neutralMultiplicity()) > 1;
            if( fabs(jet->eta()) <= 2.4 ){
               id_24 = jet->chargedHadronEnergyFraction() > 0 and jet->chargedMultiplicity() > 0
                  and jet->chargedEmEnergyFraction() < 0.99;
            }
         } else {
            id_forw = jet->neutralEmEnergyFraction() < 0.9 and jet->neutralMultiplicity() > 10;
         }
         jet_passid.push_back( id_24 and id_30 and id_forw );

         //jet_corrL1.push_back( corrL1->correction(*jet, iEvent, iSetup) );
         //jet_corrL123.push_back( corrL123->correction(*jet, iEvent, iSetup) );

      }else{

         // add the (corrected) jet to the met_sumpt
         // (was subtracted previously)
         cand_pt.push_back( jet->pt() );
         cand_energy.push_back( jet->energy() );
         cand_phi.push_back( jet->phi() );
         cand_eta.push_back( jet->eta() );
         met_sumpt += jpt;
         met_sumpt_inputs++;

      }
   }
   int i = 0;
   x_tot = 0;
   y_tot = 0;
   for ( std::vector<double>::const_iterator candpt = cand_pt.begin(); candpt != cand_pt.end(); ++candpt ) {
             //std::cout << (*candpt) << std::endl;
             cand_x.push_back(cand_pt.at(i) * cos(cand_phi.at(i)));
             cand_y.push_back(cand_pt.at(i) * sin(cand_phi.at(i)));
             x_tot += cand_pt.at(i) * cos(cand_phi.at(i));
             y_tot += cand_pt.at(i) * sin(cand_phi.at(i));

             i++;
   }
   
   y_bar = y_tot/cand_y.size();
   x_bar = x_tot/cand_x.size();
   //std::cout << x_bar << " " << x_tot << std::endl;
   
   i = 0;
   c_xx = c_xy = c_yy = 0;
   for ( std::vector<double>::const_iterator candpt = cand_pt.begin(); candpt != cand_pt.end(); ++candpt ) {
             c_xx += (cand_x.at(i) - x_bar)*(cand_x.at(i) - x_bar);
             c_xy += (cand_x.at(i) - x_bar)*(cand_y.at(i) - y_bar);
             c_yy += (cand_y.at(i) - y_bar)*(cand_y.at(i) - y_bar);
             i++;
   }
   
   //std::cout << c_xx << " " << c_xy << " " << c_yy << std::endl;
   c_xx = (cand_pt.size()-1.)/cand_pt.size() * c_xx;
   c_xy = (cand_pt.size()-1.)/cand_pt.size() * c_xy;
   c_yy = (cand_pt.size()-1.)/cand_pt.size() * c_yy;
   
   alt_sumpt = 0;
   for ( std::vector<double>::const_iterator candpt = cand_pt.begin();
         candpt != cand_pt.end(); ++candpt ) {
            alt_sumpt += *candpt;
        }
   //iterate over candidates to get sum_pt
   //std::cout << "Inputs to met_sumpt " << met_sumpt_inputs << std::endl;
   //std::cout << "met_sumpt " << met_sumpt << std::endl;
   //std::cout << "alt_sumpt " << alt_sumpt << std::endl;
   events_total++;
   //bool pass_selection = (nmuons == 2) and (dimuon_mass > 60) and (dimuon_mass < 120);
   bool pass_selection = nMediumElectrons>0; //= (nmuons == 2) and (dimuon_mass > 80) and (dimuon_mass < 100);
   if( pass_selection ){
      results_tree -> Fill();
      events_pass++;
      //std::cout << cleanjets.size() << " / " << jets.size() << std::endl;
      //std::cout << "Muon charge prod = " << charge; 
      //if( charge > 0 ) std::cout << "****************************************************************";
      //std::cout << std::endl;
      //fflush(stdout);

   }

   //delete ptRes_;
   //delete phiRes_;

}

   std::vector<pat::Jet>
MakeEleNtuple::cleanJets(double ptThreshold, double dRmatch,
      std::vector<pat::Jet>& jets, std::vector<reco::Candidate::LorentzVector>& leptons)
{
   double dR2match = dRmatch*dRmatch;
   std::vector<pat::Jet> retVal;
   for ( std::vector<pat::Jet>::const_iterator jet = jets.begin();
         jet != jets.end(); ++jet ) {
      bool isOverlap = false;
      for ( std::vector<reco::Candidate::LorentzVector>::const_iterator lepton = leptons.begin();
            lepton != leptons.end(); ++lepton ) {
         TLorentzVector ljet, llep;
         ljet.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
         llep.SetPtEtaPhiE( lepton->pt(), lepton->eta(), lepton->phi(), lepton->energy() );
         if ( pow(ljet.DeltaR( llep ),2) < dR2match ) isOverlap = true;
      }
      //if ( jet->pt() > ptThreshold && !isOverlap ){
      if ( !isOverlap ){
         retVal.push_back(*jet);
      }
   }

   return retVal;
}


// ------------ method called once each job just before starting event loop  ------------
void 
MakeEleNtuple::beginJob()
{

   events_total=0, events_pass=0;

   OutFile__file  = new TFile( OutputFileName_.c_str(), "RECREATE" );

   // pileup reweighting
   std::vector< float > PU2015_MC;
   std::vector< float > PU2015_Data;

   for( int i=0; i<NUMPUBINS; i++) {
      PU2015_MC.push_back( PU2015_MCf[i] );
      PU2015_Data.push_back( PU2015_Dataf[i] );
   }
   LumiWeights_ = edm::LumiReWeighting( PU2015_MC, PU2015_Data);

   results_tree = new TTree("events", "events");
   results_tree -> Branch("run", &run, "run/I");
   results_tree -> Branch("lumi", &lumi, "lumi/I");
   results_tree -> Branch("event", &event, "event/I");
   results_tree -> Branch("mcweight", &mcweight, "mcweight/D");
   //results_tree -> Branch("mcweightSum", &mcweightSum, "mcweightSum/D");

   results_tree -> Branch("muon_pt", &muon_pt);
   results_tree -> Branch("muon_energy", &muon_energy);
   results_tree -> Branch("muon_phi", &muon_phi);
   results_tree -> Branch("muon_eta", &muon_eta);
   results_tree -> Branch("muon_charge", &muon_charge);
   //results_tree -> Branch("muon_OS", &charge);

   results_tree -> Branch("ele_pt", &ele_pt);
   results_tree -> Branch("ele_energy", &ele_energy);
   results_tree -> Branch("ele_phi", &ele_phi);
   results_tree -> Branch("ele_eta", &ele_eta);
   results_tree -> Branch("ele_charge", &ele_charge);
   results_tree -> Branch("ele_loose", &ele_loose);
   results_tree -> Branch("ele_medium", &ele_medium);
   results_tree -> Branch("ele_tight", &ele_tight);



   /*
   results_tree -> Branch("lep_pt", &lep_pt);
   results_tree -> Branch("lep_energy", &lep_energy);
   results_tree -> Branch("lep_phi", &lep_phi);
   results_tree -> Branch("lep_eta", &lep_eta);
   */

   results_tree -> Branch("jet_pt", &jet_pt);
   results_tree -> Branch("jet_energy", &jet_energy);
   results_tree -> Branch("jet_phi", &jet_phi);
   results_tree -> Branch("jet_eta", &jet_eta);
   results_tree -> Branch("jet_sigmapt", &jet_sigmapt);
   results_tree -> Branch("jet_sf", &jet_sf);
   results_tree -> Branch("jet_sigmaphi", &jet_sigmaphi);
   results_tree -> Branch("jet_passid", &jet_passid);
   //results_tree -> Branch("jet_corrL1", &jet_corrL1);
   //results_tree -> Branch("jet_corrL123", &jet_corrL123);

   //results_tree -> Branch("cand_pt", &cand_pt);
   //results_tree -> Branch("cand_phi", &cand_phi);

   results_tree -> Branch("cov_xx", &c_xx);
   results_tree -> Branch("cov_xy", &c_xy);
   results_tree -> Branch("cov_yy", &c_yy);
   
   results_tree -> Branch("met_pt", &met_pt);
   results_tree -> Branch("met_sig", &met_sig);
   results_tree -> Branch("met_energy", &met_energy);
   results_tree -> Branch("met_phi", &met_phi);
   results_tree -> Branch("met_eta", &met_eta);
   results_tree -> Branch("met_sumpt", &met_sumpt);

   results_tree -> Branch("nmuons", &nmuons);
   results_tree -> Branch("nLooseElectrons", &nLooseElectrons);
   results_tree -> Branch("nMediumElectrons", &nMediumElectrons);
   results_tree -> Branch("nTightElectrons", &nTightElectrons);

   results_tree -> Branch("dimuon_mass", &dimuon_mass);
   results_tree -> Branch("nvertices", &nvertices);
   results_tree -> Branch("weight_pu", &weight_pu);
   results_tree -> Branch("nvert_true", &T_nvertices);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MakeEleNtuple::endJob() 
{
   OutFile__file -> Write();
   OutFile__file -> Close();
   std::cout << " *** NUMBER OF EVENTS PASSING SELECTION *** " << std::endl;
   std::cout << " ------> " << events_pass << " / " << events_total << " = " << double(events_pass)/events_total << std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MakeNtuple::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MakeNtuple::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MakeNtuple::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MakeNtuple::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MakeEleNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeEleNtuple);
