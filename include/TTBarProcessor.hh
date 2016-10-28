#ifndef TTBarProcessor_h
#define TTBarProcessor_h 1
#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
//#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Vertex.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/PIDHandler.h>
#include <IMPL/ReconstructedParticleImpl.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <TFile.h>
#include <TTree.h>

#include "MathOperator.hh"
#include "VertexChargeOperator.hh"
#include "QQBarMCOperator.hh"
#include "TopQuark.hh"
#include "RecoJet.hh"
#include "TreeWriter.hh"
#include "TreeStructures.hh"
using namespace lcio ;
using namespace marlin ;


namespace TTBarProcessor 
{
	enum ANALYSIS_TYPE
	{
		TTBarSemileptonic = 0,
		TTBarHadronic = 1,
		BBbar = 2
	};
	class TTBarProcessor : public Processor 
	{
	  
	 public:
	  
	  virtual Processor*  newProcessor() { return new TTBarProcessor ; }
	  
	  
	  TTBarProcessor() ;
	  
	  /** Called at the begin of the job before anything is read.
	   * Use to initialize the processor, e.g. book histograms.
	   */
	  virtual void init() ;
	  
	  /** Called for every run.
	   */
	  virtual void processRunHeader( LCRunHeader* run ) ;
	  
	  /** Called for every event - the working horse.
	   */
	  virtual void processEvent( LCEvent * evt ) ; 
	  
	  
	  virtual void check( LCEvent * evt ) ; 
	  
	  
	  /** Called after data processing for clean up.
	   */
	  virtual void end() ;
	  
	  std::vector< Vertex * > * convert(const std::vector< LCObject * > & objs);
	  std::vector< RecoJet * > * getBTagJets(std::vector< RecoJet * > * alljets, std::vector< RecoJet * > * wjets = NULL);
	  std::vector< RecoJet * > * getJets(LCCollection * jetcol, LCCollection *jetrelcol);
	  std::vector< EVENT::MCParticle * > AnalyseGenerator(QQBarMCOperator & opera);
	  std::vector< EVENT::MCParticle * > AnalyseGeneratorBBBar(QQBarMCOperator & opera);

	  std::vector< RecoJet * > * formW(std::vector< RecoJet * > * wjets);
	  std::vector< TopQuark * > * composeTops(std::vector< RecoJet * > * bjets, std::vector< RecoJet * > * wjets);

	  void AnalyseTTBarSemiLeptonic(LCEvent * evt);
	  void AnalyseTTBarHadronic(LCEvent * evt);
	  void AnalyseBBBar(LCEvent * evt);

	  void Match(std::vector< EVENT::MCParticle * > & mctops, TopQuark * topHadronic,  TopQuark * top2 =NULL );
	  void MatchB(std::vector< RecoJet * > *jets, std::vector< EVENT::MCParticle * > & mcbs, LCCollection * mcvtxcol = NULL);
	  void MatchB(std::vector< EVENT::MCParticle * > & mcbs, TopQuark * topHadronic, TopQuark * top2 =NULL, LCCollection * mcvtxcol = NULL);
	  void PrintJet(RecoJet * jet);
	  void PrintJets(std::vector< RecoJet * > *jets);
	  void ComputeCharge(std::vector< RecoJet * > *jets, VertexChargeOperator & vtxOperator);
	  void ComputeCharge(TopQuark * top, TopQuark * top2);
	  void ComputeChargeLepton(TopQuark * top, TopQuark * top2);
	  void __ComputeChargeCheat(TopQuark * top, TopQuark * top2);
	  void test(TopQuark * top, TopQuark * top2);
	  float getChi2(TopQuark * c);
	  void DecideOnAsymmetry(TopQuark * top, TopQuark * top2);
	  void ClearVariables();
	  void PrintParticle(EVENT::ReconstructedParticle * particle);
	  void PrintParticle(EVENT::MCParticle * particle);
	  void ComputeChargeTVCM(TopQuark * top, TopQuark * top2, VertexChargeOperator & vtxOperator);
	  EVENT::ReconstructedParticle * findPhoton(LCCollection * pfocol);
//	  bool sortByBtag(const RecoJet &lhs, const RecoJet &rhs) {return lhs.GetBTag() < rhs.GetBTag(); }
	 protected:
	  std::string intToStr(int * number);
	  std::vector<int> getOpposite(int i, int j);
	  /** Input collection name.
	   */
	  int _type;
	  std::string _hfilename;
	  ANALYSIS_TYPE _analysisType ;
	  std::string _colName ;
	  std::string _MCColName ;
	  std::string _JetsColName ;
	  std::string _JetsRelColName ;
	  std::string _JetsVtxColName ;
	  std::string _IsoLeptonColName;
	  std::string _MCVtxColName ;
	  std::string _colRelName;
	  int _ePolarization;
	  float _massCutparameter;
	  float _Ebeamparameter;
	  float _lowBTagCutparameter;
	  float _highBTagCutparameter;
	  float _WMassparameter;
	  float _WMassSigmaparameter;
	  float _TopMassparameter;
	  float _TopMassSigmaparameter;
	  float _EBeamparameter;
	  float _EBeamSigmaparameter;
	  float _PStarparameter;
	  float _PStarSigmaparameter;
	  float _CosbWparameter;
	  float _CosbWSigmaparameter;
	  float _GammaTparameter;
	  float _GammaTSigmaparameter;
	  
	  TFile * _hfile;
	  TTree * _hTree;
	  TTree * _hGenTree;
	  TTree * _hSumTree;

	  float _totalEnergy;
	  float _missedEnergy;

	  int _nRun ;
	  SummaryData _summary;
	  StatsData _stats;
	} ;
	bool sortByBtag(RecoJet *lhs, RecoJet *rhs) {return lhs->GetBTag() > rhs->GetBTag(); }
	bool sortByCostheta(RecoJet *lhs, RecoJet *rhs) 
	{
	  	return std::abs(lhs->GetCostheta()) >  std::abs(rhs->GetCostheta()); 
	}
		
} /* TTbarAnalisys */

#endif



