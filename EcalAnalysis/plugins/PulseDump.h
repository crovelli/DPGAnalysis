#ifndef __Analyzer_PulseDump_H__
#define __Analyzer_PulseDump_H__

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

// Geometry                                                                                                
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"

#include "TTree.h"
#include <memory>

class PulseDump : public edm::one::EDAnalyzer<> {

 public:
  PulseDump(const edm::ParameterSet&);
  ~PulseDump() { }

  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

 private:
  void FillRecHitData(const EcalDataFrame& dataFrame,const EBRecHitCollection* hitsEB,const EERecHitCollection* hitsEE);
  void FillRecHitMC(const EcalDataFrame& dataFrame,float mc1px,float mc1py,float mc1pz, float mc1ene,float mc2px,float mc2py,float mc2pz, float mc2ene,const CaloSubdetectorGeometry*geometry);
  uint32_t setFlagBits(const std::vector<std::vector<uint32_t> >& map,const uint32_t& status );

  edm::ESHandle<EcalChannelStatus> chStatus;
  std::vector<int> v_chstatus_;
  uint32_t flagmask_; 
  std::vector<std::vector<uint32_t> > v_DB_reco_flags_; 

  // uint32_t flagKsatur_;  
  // uint32_t flagLErecov_; 

  std::string mcProducer_;
  
  edm::EDGetTokenT<EBDigiCollection> ebDigiCollectionToken_;
  edm::EDGetTokenT<EEDigiCollection> eeDigiCollectionToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalHitEBToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalHitEEToken_;
  const CaloSubdetectorTopology* theSubdetTopologyEB_;

  TTree *_tree; 

  int entry;

  int run_;
  int event_;
  bool barrel_;
  int ietaix_;
  int iphiiy_;
  int iz_;
  double pulse_[10];
  unsigned int gain_[10];
  unsigned int firstGZ_;
  unsigned int firstG1_;
  unsigned int firstG6_;
  int lastUnsatSample_;
  float mcEta_, mcEt_;
  bool kSaturated_; 
  bool kLErecovered_; 
  bool kOutOfTime_;
  float hitTime_;
  float swissCross_;

};

DEFINE_FWK_MODULE(PulseDump);
#endif
