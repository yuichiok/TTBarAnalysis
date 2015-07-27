#ifndef TTBarProcessor_h
#define TTBarProcessor_h 1
#include <iostream>
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
#include "MCOperator.hh"
#include "TopQuark.hh"
#include "RecoJet.hh"
using namespace lcio ;
using namespace marlin ;


namespace TTbarAnalysis 
{
	class TTBarProcessor : public Processor {
	  
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
	  void AnalyseGenerator(MCOperator & opera);
	  void PrintJet(RecoJet * jet);
	  void PrintJets(std::vector< RecoJet * > *jets);
	  void ComputeCharge(TopQuark * top, TopQuark * top2);
	  float getChi2(TopQuark * c);
	  void ClearVariables();
	 protected:

	  /** Input collection name.
	   */
	  std::string _hfilename;
	  std::string _colName ;
	  std::string _MCColName ;
	  std::string _JetsColName ;
	  std::string _JetsRelColName ;
	  std::string _JetsVtxColName ;
	  std::string _IsoLeptonColName;
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
	  
	  TFile * _hfile;
	  TTree * _hTree;
	  int _mctag;
	  int _recotag;
	  float _W1mass;
	  float _W1momentum;
	  float _W1costheta;
	  float _W2mass;
	  float _W2momentum;
	  float _W2costheta;
	  float _Top1charge;
	  float _Top1btag;
	  float _Top1bmomentum;
	  int _Top1bntracks;
	  float _Top1mass;
	  float _Top1momentum;
	  float _Top1costheta;
	  int _Top2charge;
	  int _Top2bntracks;
	  float _Top2bmomentum;
	  float _Top2btag;
	  float _Top2mass;
	  float _Top2momentum;
	  float _Top2costheta;
	  float _qCostheta[2];
	  float _chiHad;

	  float _MCTopmomentum;
	  float _MCTopcostheta;
	  float _MCTopBarmomentum;
	  float _MCTopBarcostheta;

	  int _nRun ;
	  int _nEvt ;
	} ;
		
} /* TTbarAnalisys */

#endif



