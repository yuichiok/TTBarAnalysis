#ifndef TTBarProcessor_h
#define TTBarProcessor_h 1
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
#include "MCOperator.hh"
#include "TopQuark.hh"
#include "RecoJet.hh"
#include "TreeWriter.hh"
#include "TreeStructures.hh"
using namespace lcio ;
using namespace marlin ;


namespace TTbarAnalysis 
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
	  std::vector< EVENT::MCParticle * > AnalyseGenerator(MCOperator & opera);

	  void AnalyseTTBarSemiLeptonic(LCEvent * evt);
	  void AnalyseBBBar(LCEvent * evt);

	  void Match(std::vector< EVENT::MCParticle * > & mctops, TopQuark * topHadronic,  TopQuark * top2 =NULL );
	  void MatchB(std::vector< EVENT::MCParticle * > & mcbs, TopQuark * topHadronic, TopQuark * top2 =NULL, LCCollection * mcvtxcol = NULL);
	  void PrintJet(RecoJet * jet);
	  void PrintJets(std::vector< RecoJet * > *jets);
	  void ComputeCharge(TopQuark * top, TopQuark * top2);
	  void ComputeChargeLepton(TopQuark * top, TopQuark * top2);
	  void __ComputeChargeCheat(TopQuark * top, TopQuark * top2);
	  void test(TopQuark * top, TopQuark * top2);
	  float getChi2(TopQuark * c);
	  void DecideOnAsymmetry(TopQuark * top, TopQuark * top2);
	  void ClearVariables();
	  void PrintParticle(EVENT::ReconstructedParticle * particle);
	  void ComputeChargeTVCM(TopQuark * top, TopQuark * top2, VertexChargeOperator & vtxOperator);
	 protected:
	  std::string intToStr(int * number);
	  /** Input collection name.
	   */
	  std::string _hfilename;
	  std::string _colName ;
	  std::string _MCColName ;
	  std::string _JetsColName ;
	  std::string _JetsRelColName ;
	  std::string _JetsVtxColName ;
	  std::string _IsoLeptonColName;
	  std::string _MCVtxColName ;
	  std::string _colRelName;
	  int _ePolarization;
	  float _Ebeamparameter;
	  float _lowBTagCutparameter;
	  float _highBTagCutparameter;
	  float _WMassparameter;
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
		
} /* TTbarAnalisys */

#endif



