#include "TTBarProcessor.hh"

using namespace lcio ;
using namespace marlin ;
using std::vector;
using std::string;
using std::abs;
using EVENT::Track;
using IMPL::ReconstructedParticleImpl;
using EVENT::ParticleID;
using IMPL::ParticleIDImpl;
namespace TTBarProcessor
{
	TTBarProcessor aTTBarProcessor ;


	TTBarProcessor::TTBarProcessor() : Processor("TTbarAnalisys") 
	{
	  
		// modify processor description

		// register steering parameters: name, description, class-variable, default value

		registerProcessorParameter( "ROOTFileName",
				      "Output ROOT File Name",
				      _hfilename,
				      string("TTBarProcessor.root") );
		registerProcessorParameter( "ElectronPolarisation",
				      "Helicity of e- beam",
				      _ePolarization,
				      0 );
		registerProcessorParameter( "MassCut",
				      "Mass cut on qqbar",
				      _massCutparameter,
				      (float)200. );
		registerProcessorParameter( "AnalysisType",
				      "Analysis Type",
				      _type,
				      _type);
		registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE , 
				   "PFOCollection",
				   "Name of the Calorimeter hit collection"  ,
				   _colName ,
				   string("PandoraPFOs") ) ;

		registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"JetCollectionName" , 
			"Name of the Jet collection"  ,
			_JetsColName ,
			std::string("FinalJets")
		);
		registerInputCollection( LCIO::LCRELATION,
			"JetRelCollectionName" , 
			"Name of the PrimaryVertex collection"  ,
			_JetsRelColName ,
			std::string("FinalJets_rel")
		);
		registerInputCollection( LCIO::VERTEX,
			"GenVtxCollectionName" , 
			"Name of the PrimaryVertex collection"  ,
			_MCVtxColName ,
			std::string("MCVertex")
		);
		registerInputCollection( LCIO::MCPARTICLE,
			"MCCollectionName" , 
			"Name of the MC collection"  ,
			_MCColName ,
			std::string("MCParticlesSkimmed")
		);
		registerInputCollection( LCIO::MCPARTICLE,
			"IsoLeptonCollectionName" , 
			"Name of the isolepton collection"  ,
			_IsoLeptonColName ,
			std::string("SelectedLepton")
		);
		registerInputCollection( LCIO::LCRELATION,
		    "RecoMcTruthCollectionName" , 
		    "Name of the RecoMcTruth collection"  ,
		    _colRelName ,
		    std::string("RecoMCTruthLink")
		);
		_type = 0;
		
		std::cout << "Type: " <<  _type << "\n";
		_analysisType = static_cast<ANALYSIS_TYPE>(_type);
		std::cout << "Type: " <<  _analysisType << "\n";

		_lowBTagCutparameter = 0.3;
		_highBTagCutparameter = 0.8;
		_WMassparameter = 80.3;
		_TopMassparameter = 174;
		_TopMassSigmaparameter = 6.3;
		_EBeamparameter = 250.0;
		_EBeamSigmaparameter = 8;
		_PStarparameter = 68.0;
		_PStarSigmaparameter = 5;
		_CosbWparameter = 0.23;
		_CosbWSigmaparameter = 0.14;
		_GammaTparameter = 1.435;
		_GammaTSigmaparameter = 0.05;
	}
	
	void TTBarProcessor::init() 
	{ 
		// usually a good idea to
		printParameters() ;
		std::cout << "_analysisType: " << _analysisType << "\n";
		_nRun = 0 ;
		TreeWriter writer;
		_hfile = new TFile( _hfilename.c_str(), "RECREATE", _hfilename.c_str() ) ;
		_hSumTree = new TTree( "Summary", "tree" );
		writer.InitializeSummaryTree(_hSumTree, _summary);
		_hGenTree = new TTree( "GenTree", "tree" );
		_hGenTree->Branch("qMCcostheta", _stats._qMCcostheta, "qMCcostheta[2]/F");
		_hGenTree->Branch("MCMass", &_stats._MCMass, "MCMass/F");
		_hGenTree->Branch("MCPDG", &_stats._MCPDG, "MCPDG/F");
		_hGenTree->Branch("MCquarkAngle", &_stats._MCquarkAngle, "MCMass/F");
		_hGenTree->Branch("MCBcostheta", _stats._qMCBcostheta, "MCBcostheta[2]/F");
		_hGenTree->Branch("MCLeptonPDG", &_stats._MCLeptonPDG, "MCLeptonPDG/I");
		_hTree = new TTree( "Stats", "tree" );
		switch(_analysisType)
		{
			case TTBarSemileptonic:
				writer.InitializeStatsTree(_hTree, _stats);
				break;
			case BBbar:
				writer.InitializeStatsBBBarTree(_hTree, _stats);
				break;
			case TTBarHadronic:
				writer.InitializeStatsTree(_hTree, _stats);
				break;
		}
		/*_hBkgTree = new TTree( "BkgTree", "tree" );
		_hBkgTree->Branch("f1PDG", &_bkgPDGs[0], "f1PDG/I");
		_hBkgTree->Branch("f2PDG", &_bkgPDGs[1], "f2PDG/I");
		_hBkgTree->Branch("f3PDG", &_bkgPDGs[2], "f3PDG/I");
		_hBkgTree->Branch("f4PDG", &_bkgPDGs[3], "f4PDG/I");*/
		//AnalyseTTBarSemiLeptonic(evt);

		
		//_hTree->Branch("detectedTotal", &_stats._detectedTotal, "detectedTotal/I");
	}
	void TTBarProcessor::processRunHeader( LCRunHeader* run) 
	{ 
		_nRun++ ;
	} 
	vector< MCParticle * > TTBarProcessor::AnalyseGeneratorBBBar(QQBarMCOperator & opera)
	{
		vector <MCParticle *> genfinalstate;// = opera.GetPairParticles(5);
		for (unsigned int i = 5; i > 0; i--) 
		{
			genfinalstate = opera.GetPairParticles(i);
			if (genfinalstate.size() == 2) 
			{
				break;
			}
		}
		for (unsigned int i = 0; i < genfinalstate.size(); i++) 
		{
			PrintParticle(genfinalstate[i]);
		}
		if (genfinalstate.size() == 2) 
		{
			MCParticle * Z = opera.CombineParticles(genfinalstate[0],genfinalstate[1]);
			_stats._MCMass = Z->getMass();
			_stats._MCPDG = abs(genfinalstate[0]->getPDG());
			_stats._MCquarkAngle =  MathOperator::getAngleBtw(genfinalstate[0]->getMomentum(), genfinalstate[1]->getMomentum());
			_stats._MCPt = 0.0;
			for (unsigned int i = 0; i < 3; i++) 
			{
				_stats._MCPt += std::pow(genfinalstate[0]->getMomentum()[i]+genfinalstate[1]->getMomentum()[i],2);
			}
			_stats._MCPt = std::sqrt(_stats._MCPt);
			for (unsigned int i = 0; i < genfinalstate.size(); i++) 
			{
				vector<float> direction = MathOperator::getDirection(genfinalstate[i]->getMomentum());
				int charge = genfinalstate[i]->getCharge() / std::abs(genfinalstate[i]->getCharge());
				_stats._qMCcostheta[i] = - std::cos( MathOperator::getAngles(direction)[1]) * charge;
			}
			if (Z->getMass() > _massCutparameter && _stats._MCPDG == 5) 
			{
				_summary._nAfterGenMassCuts++;
			}
			std::cout << "Generated mass: " << Z->getMass() << " Sp: " << _stats._MCPt << " angle: " << _stats._MCquarkAngle << "\n";
			if (_stats._MCPDG == 5) 
			{
				_summary._nGenUsed++;
			}
			_hGenTree->Fill();
			delete Z;
		}
		/*vector <MCParticle *> bkgfinalstate = opera.GetFinalStateBkg();
		for (unsigned int i = 0; i < bkgfinalstate.size(); i++) 
		{
			_bkgPDGs[i] = bkgfinalstate[i]->getPDG();
		}
		_hBkgTree->Fill();*/
		return genfinalstate;
		
	}
	vector< MCParticle * > TTBarProcessor::AnalyseGenerator(QQBarMCOperator & opera)
	{
		std::vector< EVENT::MCParticle * > mctops = opera.GetTopPairParticles(_stats._MCTopBangle, _stats._MCTopcosWb);

		if (mctops.size() < 2) 
		{
			_stats._mctag = 0;
			return mctops;
		}
		_summary._nGenUsed++;
		MCParticle * top = mctops[0];
		MCParticle * topbar = mctops[1];
		vector<float> direction = MathOperator::getDirection(top->getMomentum());
		vector<float> directionbar = MathOperator::getDirection(topbar->getMomentum());
		_stats._MCTopmass = top->getMass();
		_stats._MCTopBarmass = topbar->getMass();
		_stats._MCTopmomentum = MathOperator::getModule(top->getMomentum());
		_stats._MCTopBarmomentum = MathOperator::getModule(topbar->getMomentum());
		_stats._MCTopcostheta = std::cos( MathOperator::getAngles(direction)[1] );
		_stats._MCTopBarcostheta = std::cos( MathOperator::getAngles(directionbar)[1] );
		_stats._MCquarkAngle =  MathOperator::getAngleBtw(top->getMomentum(), topbar->getMomentum());
		std::cout << "Top costheta: " << _stats._MCTopcostheta << "\n";
		std::cout << "TopBar costheta: " << _stats._MCTopBarcostheta << "\n";
		if (_stats._MCTopmass > 170.0 && _stats._MCTopmass < 181.0 && _stats._MCTopBarmass > 170.0 && _stats._MCTopBarmass < 181.0) 
		{
			_stats._mctag = 1;
		}
		else 
		{
			_stats._mctag = 0;
		}
		vector< EVENT::MCParticle * > mcb = opera.GetBquarkPair();
		vector<float> bdirection = MathOperator::getDirection(mcb[0]->getMomentum());
		vector<float> bdirectionbar = MathOperator::getDirection(mcb[1]->getMomentum());
		_stats._qMCBcostheta[0] =  std::cos( MathOperator::getAngles(bdirection)[1] );
		_stats._qMCBcostheta[1] =  std::cos( MathOperator::getAngles(bdirectionbar)[1] );
		_stats._qMCcostheta[0] = _stats._MCTopcostheta;//(std::abs(_stats._qMCBcostheta[0]) < 0.9 )? _stats._MCTopcostheta: -2;
		_stats._qMCcostheta[1] = -_stats._MCTopBarcostheta;//(std::abs(_stats._qMCBcostheta[1]) < 0.9 )? -_stats._MCTopBarcostheta : -2;
		_hGenTree->Fill();
		vector <MCParticle *> final = opera.GetFinalState();
		for (unsigned int i = 0; i < final.size(); i++) 
		{
			if (abs(final[i]->getPDG()) == 11 || abs(final[i]->getPDG()) == 13 || abs(final[i]->getPDG()) == 15) 
			{
				_stats._MCLeptonMomentum = MathOperator::getModule(final[i]->getMomentum());
				vector<float> ldirection = MathOperator::getDirection(final[i]->getMomentum());
				_stats._MCLeptonCostheta = std::cos(MathOperator::getAngles(ldirection)[1]);
				_stats._MCLeptonPDG = final[i]->getPDG();
			}
		}
		if (abs(_stats._MCLeptonPDG) == 15) 
		{
			MCParticle * newlepton = opera.GetTauLepton();
			if (newlepton) 
			{
				_stats._MCLeptonMomentum = MathOperator::getModule(newlepton->getMomentum());
				vector<float> ldirection = MathOperator::getDirection(newlepton->getMomentum());
				_stats._MCLeptonCostheta = std::cos(MathOperator::getAngles(ldirection)[1]);
			}
			else 
			{
				_stats._MCLeptonMomentum = -1;
				_stats._MCLeptonCostheta = -1;
			}
		}
		MCParticle * neutrino = opera.GetNeutrino();
		if (neutrino) 
		{
			_stats._MCNeutrinoEnergy = neutrino->getEnergy();
		}
		return mctops;
	}
	void TTBarProcessor::processEvent( LCEvent * evt )
	{
		switch(_analysisType)
		{
			case TTBarSemileptonic:
				AnalyseTTBarSemiLeptonic(evt);
				break;
			case BBbar:
				AnalyseBBBar(evt);
				break;
			case TTBarHadronic:
				AnalyseTTBarHadronic(evt);
				break;
		}
		//AnalyseTTBarSemiLeptonic(evt);
	}

	void TTBarProcessor::AnalyseTTBarSemiLeptonic( LCEvent * evt )
	{
		LCCollection * mcvtxcol = NULL;
		try
		{
			 mcvtxcol = evt->getCollection(_MCVtxColName);
		}
		catch(DataNotAvailableException &e)
		{
			std::cout << e.what() <<"\n";
		}
		
		try
		{
			std::cout << "***************Analysis*" <<_summary._nEvt++<<"*****************\n";

			LCCollection * isoleptoncol = evt->getCollection(_IsoLeptonColName);
			LCCollection * jetcol = evt->getCollection(_JetsColName);
			LCCollection * jetrelcol = evt->getCollection(_JetsRelColName);
			LCCollection * mccol = evt->getCollection(_MCColName);
			//LCCollection * mcvtxcol = evt->getCollection(_MCVtxColName);
			QQBarMCOperator opera(mccol);
			VertexChargeOperator vtxOperator(evt->getCollection(_colName),evt->getCollection(_colRelName));
			vector < MCParticle * > mctops = AnalyseGenerator(opera);
			vector< RecoJet * > * jets = getJets(jetcol, jetrelcol);
			vector< RecoJet * > * wjets = new vector< RecoJet * >();
			vector< RecoJet * > * bjets = getBTagJets(jets, wjets);
			LCCollection * pfocol = evt->getCollection(_colName);
			_stats._Thrust = pfocol->getParameters().getFloatVal("majorThrustValue");
			std::cout << "B jets: \n";
			PrintJets(bjets);
			std::cout << "W jets: \n";
			PrintJets(wjets);
			if (!isoleptoncol || isoleptoncol->getNumberOfElements() == 0) 
			{
				std::cout << "No IsoLepton!\n";
				if (isoleptoncol) 
				{
					int IsTau = isoleptoncol->parameters().getIntVal("isLeptonFromTau");
					if (IsTau) 
					{
						std::cout << " IsoLepton is tau! \n";
					}
				}
				else 
				{
					std::cout << " No IsoLepton collection!!! \n";
				}
				return;
			}
			ReconstructedParticle * wLeptonic = NULL;
			std::cout << "\nIsoLeptons:\n";
			for (int i = 0; i < isoleptoncol->getNumberOfElements(); i++) 
			{
				ReconstructedParticle * particle = dynamic_cast< ReconstructedParticle * >(isoleptoncol->getElementAt(i));
				if (particle->getType() != 22) 
				{
					wLeptonic = particle;
				}
				PrintParticle(particle);
			}
			if (!wLeptonic) 
			{
				std::cout << "No IsoLepton!!! \n";
				return;
			}
			_summary._nAfterLeptonCuts++;
			if (!bjets || bjets->size() != 2 ) 
			{
				std::cout << "No B jets!!! \n";
				return;
			}
			if ((bjets->at(0)->GetBTag() <  _lowBTagCutparameter && bjets->at(1)->GetBTag() <  _lowBTagCutparameter ) ||bjets->at(0)->GetBTag() < _highBTagCutparameter ) 
			{
				std::cout << "No high B jets!!! \n";
				return;
			}
			_summary._nAfterBtagCuts++;
			float chimin = 100000.0;
			RecoJet * wHadronic = new TopQuark(wjets->at(0), wjets->at(1));
			vector<float> direction = MathOperator::getDirection(wHadronic->getMomentum());
			_stats._W1costheta = std::cos (MathOperator::getAngles(direction)[1]);
			//if (isoleptoncol->getNumberOfElements() > 1) 
			std::cout << '\n';
			vector<float> directionL = MathOperator::getDirection(wLeptonic->getMomentum());
			_stats._W2costheta = std::cos (MathOperator::getAngles(directionL)[1]);
			TopQuark * topHadronic = NULL;
			TopQuark * topLeptonic = NULL;
			int chosen = -1;
			for (unsigned int i = 0; i < bjets->size(); i++) 
			{
				TopQuark * candidate = new TopQuark(bjets->at(i), wHadronic);
				float chi2 = getChi2(candidate);
				if (chi2 < chimin) 
				{
					topHadronic = candidate;
					chimin = chi2;
					chosen = i;
					_stats._Top1cosWb = std::cos( MathOperator::getAngle(candidate->GetB()->getMomentum(), candidate->GetW()->getMomentum()) );
				}
			}
			if (chosen && wLeptonic) 
			{
				topLeptonic = new TopQuark(bjets->at(0), wLeptonic);
			}
			if (!chosen && wLeptonic) 
			{
				topLeptonic = new TopQuark(bjets->at(1), wLeptonic);
			}
			/*if (chosen) 
			{
				topLeptonic = new TopQuark(bjets->at(0));
			}
			if (!chosen) 
			{
				topLeptonic = new TopQuark(bjets->at(1));
			}*/
			_stats._W1mass = wHadronic->getMass();
			_stats._W1momentum = MathOperator::getModule(wHadronic->getMomentum());
			_stats._W2momentum = MathOperator::getModule(wLeptonic->getMomentum());
			_stats._Top1mass = topHadronic->getMass();
			_stats._Top1energy = topHadronic->getEnergy();
			_stats._chiHad = getChi2(topHadronic);
			float momentum[3];
			for (unsigned int i = 0; i < 3; i++) 
			{
				momentum[i] = topHadronic->getMomentum()[i]+topLeptonic->GetB()->getMomentum()[i];
			}
			_stats._hadMass = std::sqrt(pow(topHadronic->getEnergy()+topLeptonic->GetB()->getEnergy(),2) - momentum[0]* momentum[0]-momentum[1]*momentum[1]-momentum[2]*momentum[2]);
			vector< EVENT::MCParticle * > mcbquarks = opera.GetBquarkPair();
			vector< EVENT::MCParticle * > mcws = opera.GetWPair();
			Match(mctops, mcbquarks, mcws, topHadronic);
			MatchB(mcbquarks, topHadronic, topLeptonic, mcvtxcol);
			ComputeCharge(topHadronic, topLeptonic);
			ComputeChargeTVCM(topHadronic, topLeptonic, vtxOperator);

			DecideOnAsymmetry(topHadronic, topLeptonic);
			if (_stats._methodUsed && wHadronic->getMass() < 110 && topHadronic->getMass() < 200) 
			{
				_summary._nAfterMassCuts++;
			}
			test(topHadronic, topLeptonic, bjets, wjets, wLeptonic);
			//__ComputeChargeCheat(topHadronic, topLeptonic);
			_hTree->Fill();
			ClearVariables();
		}
		catch(DataNotAvailableException &e)
		{
			std::cout << e.what() <<"\n";
		}
	}
	void TTBarProcessor::ComputeChargeTVCM(TopQuark * top, TopQuark * top2, VertexChargeOperator & vtxOperator)
	{
		vtxOperator.GetAsymmetryTVCM(top, top2);
		_stats._Top1Kaon = (top->GetComputedCharge().ByTVCM )? 1:0;
		_stats._Top2Kaon = (top2->GetComputedCharge().ByTVCM )? 1:0;
		_stats._UsedBTVCM = vtxOperator.GetResultingB();
		vector< ReconstructedParticle * > kaons1 = vtxOperator.GetKaons(top);
		vector< ReconstructedParticle * > kaons2 = vtxOperator.GetKaons(top2);
		_stats._Top1KaonNumber = kaons1.size();
		_stats._Top2KaonNumber = kaons2.size();
		if (top->GetComputedCharge().ByTVCM || top2->GetComputedCharge().ByTVCM) 
		{
			_summary._nKaons++;
		}
		//_summary._nKaons += vtxOperator.CountKaons(topHadronic, topLeptonic);
		//std::cout << "Kaons 1 :\n";
		for (unsigned int i = 0; i < kaons1.size(); i++) 
		{
			_stats._Top1KaonCharges[i] = kaons1[i]->getCharge();
			_stats._Top1KaonMomentum[i] = MathOperator::getModule(kaons1[i]->getMomentum());
			//std::cout << "\tq: " <<  kaons1[i]->getCharge() << " p: " << MathOperator::getModule(kaons1[i]->getMomentum()) <<"\n";
		}
		//std::cout << "Kaons 2 :\n";
		for (unsigned int i = 0; i < kaons2.size(); i++) 
		{
			_stats._Top2KaonCharges[i] = kaons2[i]->getCharge();
			_stats._Top2KaonMomentum[i] = MathOperator::getModule(kaons2[i]->getMomentum());
			//std::cout << "\tq: " <<  kaons2[i]->getCharge() << " p: " << MathOperator::getModule(kaons2[i]->getMomentum()) <<"\n";
		}
	}
	void TTBarProcessor:: check( LCEvent * evt ) 
	{
	}
	void TTBarProcessor::test(TopQuark * top, TopQuark * top2, vector< RecoJet * > * bjets, vector< RecoJet * > * wjets, ReconstructedParticle * lepton)
	{	
		_stats._totalEnergy = top->getEnergy() + top2->getEnergy();
		_stats._missedEnergy = 2*_EBeamparameter - _stats._totalEnergy;
		float momentum[3];
		float wmomentum[3];
		float sum = 0;
		float sum2 = 0;
		float sum3 = 0;
		for (unsigned int i = 0; i < 3; i++) 
		{
			momentum[i] = (-1)*(bjets->at(0)->getMomentum()[i] + bjets->at(1)->getMomentum()[i] + wjets->at(1)->getMomentum()[i] + wjets->at(0)->getMomentum()[i] + lepton->getMomentum()[i]);
			sum += std::pow(momentum[i] + lepton->getMomentum()[i],2);
			wmomentum[i] = momentum[i] + lepton->getMomentum()[i];
			sum2 +=  std::pow(wmomentum[i] + top2->GetB()->getMomentum()[i],2);
			sum3 += std::pow(top->GetB()->getMomentum()[i] + top2->GetB()->getMomentum()[i] ,2);
		}
		float energy = MathOperator::getModule(momentum);
		std::cout << "Np: " << MathOperator::getModule(momentum) << " diff: " << _stats._missedEnergy - MathOperator::getModule(momentum) << "\n";
		_stats._W2mass = std::sqrt( pow(energy + lepton->getEnergy(),2) - sum);
		float wenergy =std::sqrt( _stats._W2mass * _stats._W2mass + std::pow(MathOperator::getModule(wmomentum),2));
 		_stats._Top2mass = std::sqrt( pow(wenergy + top2->GetB()->getEnergy(),2) - sum2);
		_stats._Top2energy = std::sqrt( _stats._Top2mass * _stats._Top2mass  +  sum2);
		_stats._Top2gamma = _stats._Top2energy / _stats._Top2mass;
		std::cout << "W2m: " << _stats._W2mass << "\n";
		_stats._bProduct = std::sqrt(sum3);
		_stats._chiTop2Mass = std::pow(_stats._Top2mass - _TopMassparameter, 2) / std::pow( _TopMassSigmaparameter, 2) ;
		_stats._chiTop2E = std::pow(_stats._Top2energy - _EBeamparameter , 2) / std::pow( _EBeamSigmaparameter, 2); 
		_stats._chiGammaT2 = std::pow( _stats._Top2gamma - _GammaTparameter, 2) / std::pow( _GammaTSigmaparameter, 2); 
	}
	void TTBarProcessor::ComputeCharge(TopQuark * top1, TopQuark * top2)
	{
		//Top hadronic
		_stats._qCostheta[0] = -2.;
		_stats._Top1bcharge = top1->GetHadronCharge();
		_stats._Top1bmomentum = top1->GetHadronMomentum();
		_stats._Top1bdistance = top1->GetMinHadronDistance();
		vector<float> direction = MathOperator::getDirection(top1->getMomentum());
		_stats._Top1costheta =  std::cos( MathOperator::getAngles(direction)[1] );
		_stats._Top1bntracks = top1->GetNumberOfVertexParticles();
		_stats._Top1btag = top1->GetBTag();
		_stats._Top1gamma = top1->getEnergy()/top1->getMass();
		_stats._W1gamma = top1->GetW()->getEnergy()/top1->GetW()->getMass();
		float Top1nvtx = top1->GetNumberOfVertices();
		vector<float> bdirection = MathOperator::getDirection(top1->GetB()->getMomentum());
		_stats._Top1bcostheta =std::abs( std::cos( MathOperator::getAngles(bdirection)[1] ));
		std::cout << "Top charge: " << _stats._Top1bcharge
			  << " top pB: " << _stats._Top1bmomentum
			  << " top W: " <<  MathOperator::getModule(top1->GetW()->getMomentum())//_stats._Top1bmomentum
			  << " top bntracks: " << _stats._Top1bntracks
			  << " top btag: " << _stats._Top1btag
			  << " top bcostheta: " << _stats._Top1bcostheta
			  << "\n";
		//Top Leptonic
		_stats._Top2bmomentum = top2->GetHadronMomentum();
		_stats._Top2bdistance = top2->GetMinHadronDistance();
		_stats._Top2bcharge = top2->GetHadronCharge();
		vector<float> direction2 = MathOperator::getDirection(top2->getMomentum());
		_stats._Top2costheta =  std::cos( MathOperator::getAngles(direction2)[1] );
		_stats._Top2bntracks = top2->GetNumberOfVertexParticles();
		_stats._Top2leptonCharge = top2->GetW()->getCharge();
		vector<float> ldirection = MathOperator::getDirection(top2->GetW()->getMomentum());
		_stats._Top2leptonE = top2->GetW()->getEnergy();
		_stats._Top2leptonCos = std::abs( std::cos( MathOperator::getAngles(ldirection)[1] ));
		_stats._Top2btag = top2->GetBTag();
		vector<float> bdirection2 = MathOperator::getDirection(top2->GetB()->getMomentum());
		_stats._Top2bcostheta =std::abs( std::cos( MathOperator::getAngles(bdirection2)[1] ));
		float Top2nvtx = top2->GetNumberOfVertices();
		std::cout << "Top charge: " << _stats._Top2bcharge
			  << " top pB: " << _stats._Top2bmomentum
			  << " top W: " << MathOperator::getModule(top2->GetW()->getMomentum())
			  << " top bntracks: " << _stats._Top2bntracks
			  << " top btag: " << _stats._Top2btag
			  << " top bcostheta: " << _stats._Top2bcostheta
			  << "\n";
		float pcut = 0;
		bool trustTop1 = false;
		bool trustTop2 = false;
		if (//_stats._Top1bmomentum > pcut && 
		    abs( _stats._Top1bcharge) > 0 &&
		    //_stats._Top1bmomentum > pcut
		    _stats._Top1bntracks > 0 //&& 
		    //_stats._Top1bntracks < 9 //&& 


		    )  //  _stats._Top1bntracks < 7 
		{
			_stats._Top1Vtx = 1;
			JetCharge & topCharge = top1->GetComputedCharge();
			topCharge.ByTrackCount = new int(top1->GetHadronCharge());
		}
		if (//_stats._Top2bmomentum > pcut2 && 
		    abs( _stats._Top2bcharge) > 0 && 
		    //_stats._Top2bmomentum > pcut
		    _stats._Top2bntracks > 0 //&& 
		    //_stats._Top2bntracks < 9 //&& 
		    ) 
		{
			_stats._Top2Vtx = 1;
			JetCharge & topCharge = top2->GetComputedCharge();
			topCharge.ByTrackCount = new int(top2->GetHadronCharge());
		}
		if ( _stats._Top2bcharge != 0 ||  _stats._Top1bcharge != 0) 
		{
			_summary._nChargedB++;
		}
		if (top2->GetW()) 
		{
			top2->GetComputedCharge().ByLepton = new int(top2->GetW()->getCharge());
		}
		bool useLepton = (_ePolarization > 0)? false : true; //Lepton by default
	}
	void TTBarProcessor::DecideOnAsymmetry(TopQuark * top1, TopQuark * top2)
	{
		//Print
		
		std::cout << "\t\tTracks\tTVCM\tLepton\tp\n";
		std::cout << "Top1:\t" << intToStr(top1->GetComputedCharge().ByTrackCount) <<"\t"
					<< intToStr(top1->GetComputedCharge().ByTVCM) <<"\t"
					<< intToStr(top1->GetComputedCharge().ByLepton) <<"\t"
					<< _stats._Top1bmomentum <<"\n";
		std::cout << "Top2:\t" << intToStr(top2->GetComputedCharge().ByTrackCount) <<"\t"
					<< intToStr(top2->GetComputedCharge().ByTVCM) <<"\t"
					<< intToStr(top2->GetComputedCharge().ByLepton) <<"\t"
					<< _stats._Top2bmomentum <<"\n";
		_stats._qCostheta[0] = -2.;

		vector<float> direction = MathOperator::getDirection(top1->getMomentum());
		float costheta =  std::cos( MathOperator::getAngles(direction)[1] );
		//float m = top1->getMass();
		//float e = top1->getEnergy();
		float btagcut = 0.8;
		float pcut = 25.;
		vector<int> samecharge;
		vector<int> goodcharge;
		vector<int> chargevalue;
		//if (e/m > _GammaTparameter -0.1 && top2->GetComputedCharge().ByLepton) 
		//Track charge * Track charge
		float gammacut1 = 1.23;
		if (top2->GetComputedCharge().ByTrackCount && top1->GetComputedCharge().ByTrackCount) 
		{
			int top1charge = *(top1->GetComputedCharge().ByTrackCount );
			int top2charge = *(top2->GetComputedCharge().ByTrackCount );
			if (top1charge * top2charge < 0 && ((_stats._Top1btag > btagcut && _stats._Top1bmomentum > pcut) ||( _stats._Top2btag > btagcut && _stats._Top2bmomentum > pcut))) 
			{
				//_stats._qCostheta[0] = (top1charge < 0)? costheta: - costheta;
				std::cout << "Two vertices are used!\n";
				//_stats._methodUsed = 1;
				chargevalue.push_back(top1charge);
				goodcharge.push_back(1);
				_stats._methodCorrect = top1->__GetMCCharge() * top1charge < 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				//_summary._nAfterKinematicCuts++;
				//return;
			}
			if (top1charge * top2charge > 0 && ((_stats._Top1btag > btagcut && _stats._Top1bmomentum > pcut) || (_stats._Top2btag > btagcut && _stats._Top2bmomentum > pcut)))
			{
				samecharge.push_back(1);
			}
		}
		
		//Kaon * Kaon
		if (top2->GetComputedCharge().ByTVCM && top1->GetComputedCharge().ByTVCM) 
		{
			int top1charge = *(top1->GetComputedCharge().ByTVCM );
			int top2charge = *(top2->GetComputedCharge().ByTVCM );
			if (top1charge * top2charge < 0) 
			{
				//_stats._qCostheta[0] = (top1charge < 0)? costheta: - costheta;
				std::cout << "Two kaons are used!\n";
				chargevalue.push_back(top1charge);
				//_stats._methodUsed = 2;
				goodcharge.push_back(2);
				_stats._methodCorrect = top1->__GetMCCharge() * top1charge < 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				///_summary._nAfterKinematicCuts++;
				//return;
			}
			if (top1charge * top2charge > 0)
			{
				samecharge.push_back(2);
			}
		}//*/
		//Track charge + Kaon
		if (top1->GetComputedCharge().ByTrackCount && top1->GetComputedCharge().ByTVCM) 
		{
			int top1charge = *(top1->GetComputedCharge().ByTrackCount );
			int top1kaon = *(top1->GetComputedCharge().ByTVCM );
			if (top1charge * top1kaon > 0 && _stats._Top1btag > btagcut && _stats._Top1bmomentum > pcut) 
			{
				//_stats._qCostheta[0] = (top1charge < 0)? costheta: -costheta; 
				std::cout << "Vertex + kaon for top1 is used!\n";
				//_stats._methodUsed = 3;
				chargevalue.push_back(top1charge);
				goodcharge.push_back(3);
				_stats._methodCorrect = top1->__GetMCCharge() * top1charge < 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				//_summary._nAfterKinematicCuts++;
				//return;
			}
			if (top1charge * top1kaon < 0 && _stats._Top1btag > btagcut && _stats._Top1bmomentum > pcut) 
			{
				samecharge.push_back(3);
			}
		}
		if (top2->GetComputedCharge().ByTrackCount && top2->GetComputedCharge().ByTVCM) 
		{
			int top2charge = *(top2->GetComputedCharge().ByTrackCount );
			int top2kaon = *(top2->GetComputedCharge().ByTVCM );
			if (top2charge * top2kaon > 0 && _stats._Top2btag > btagcut && _stats._Top2bmomentum > pcut) 
			{
				//_stats._qCostheta[0] = (top2charge > 0)? costheta: - costheta; 
				std::cout << "Vertex + kaon for top2 is used!\n";
				//_stats._methodUsed = 4;
				goodcharge.push_back(3);
				chargevalue.push_back(-top2charge);
				_stats._methodCorrect = top1->__GetMCCharge() * top2charge > 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				//_summary._nAfterKinematicCuts++;
				//return;
			}
			if (top2charge * top2kaon < 0 && _stats._Top2btag > btagcut && _stats._Top2bmomentum > pcut) 
			{
				samecharge.push_back(3);
			}
		}
		if (top1->GetComputedCharge().ByTrackCount && top2->GetComputedCharge().ByTVCM) 
		{
			int top1charge = *(top1->GetComputedCharge().ByTrackCount );
			int top2kaon = *(top2->GetComputedCharge().ByTVCM );
			if (top1charge * top2kaon < 0 &&  _stats._Top1btag > btagcut && _stats._Top1bmomentum > pcut) 
			{
				//_stats._qCostheta[0] = (top1charge < 0)? costheta: - costheta; 
				std::cout << "Vertex1 + kaon2 is used!\n";
				//_stats._methodUsed = 5;
				chargevalue.push_back(top1charge);
				goodcharge.push_back(4);
				_stats._methodCorrect = top1->__GetMCCharge() *  top1charge < 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				//_summary._nAfterKinematicCuts++;
				//return;
			}
			if (top1charge * top2kaon > 0 &&  _stats._Top1btag > btagcut && _stats._Top1bmomentum > pcut)
			{
				samecharge.push_back(4);
			}
		}
		if (top2->GetComputedCharge().ByTrackCount && top1->GetComputedCharge().ByTVCM) 
		{
			int top2charge = *(top2->GetComputedCharge().ByTrackCount );
			int top1kaon = *(top1->GetComputedCharge().ByTVCM );
			if (top2charge * top1kaon < 0 &&  _stats._Top2btag > btagcut && _stats._Top2bmomentum > pcut) 
			{
				//_stats._qCostheta[0] = (top2charge > 0)? costheta: - costheta; 
				std::cout << "Vertex2 + kaon1 is used!\n";
				//_stats._methodUsed = 6;
				goodcharge.push_back(4);
				chargevalue.push_back(-top2charge);
				_stats._methodCorrect = top1->__GetMCCharge() * top1kaon< 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				//_summary._nAfterKinematicCuts++;
				//return;
			}
			if (top2charge * top1kaon > 0 &&  _stats._Top2btag > btagcut && _stats._Top2bmomentum > pcut) 
			{
				samecharge.push_back(4);
			}
		}//*/
		//LEPTON
		if (top2->GetComputedCharge().ByLepton && top2->GetComputedCharge().ByTrackCount)
		{
			int top2lepton = *(top2->GetComputedCharge().ByLepton );
			int top2charge = *(top2->GetComputedCharge().ByTrackCount );
			if (top2charge * top2lepton < 0  &&  ((_stats._Top2btag > btagcut && _stats._Top2bmomentum > pcut) || _stats._Top1gamma > gammacut1))
			{
				std::cout << "Vertex + lepton for top2 is used!\n";
				goodcharge.push_back(5);
				chargevalue.push_back(-top2charge);
				_stats._methodCorrect = top1->__GetMCCharge() * top2charge < 0;
				//_summary._nAfterKinematicCuts++;
			}
			if (top2charge * top2lepton > 0 &&  _stats._Top2btag > btagcut && _stats._Top2bmomentum > pcut) 
			{
				samecharge.push_back(5);
			}
		}
		if (top2->GetComputedCharge().ByLepton && top1->GetComputedCharge().ByTrackCount)
		{
			int top2lepton = *(top2->GetComputedCharge().ByLepton );
			int top1charge = *(top1->GetComputedCharge().ByTrackCount );
			if (top1charge * top2lepton > 0 &&  ((_stats._Top1btag > btagcut && _stats._Top1bmomentum > pcut)  || _stats._Top1gamma > gammacut1))
			{
				std::cout << "Vertex + lepton for top1 is used!\n";
				goodcharge.push_back(5);
				chargevalue.push_back(top1charge);
				_stats._methodCorrect = top1->__GetMCCharge() * top1charge < 0;
			}
			if (top1charge * top2lepton < 0 &&  ((_stats._Top1btag > btagcut && _stats._Top1bmomentum > pcut) ||_stats._Top1gamma > gammacut1)) 
			{
				samecharge.push_back(5);
			}
		}//
		if (top2->GetComputedCharge().ByLepton && top1->GetComputedCharge().ByTVCM)
		{
			int top2lepton = *(top2->GetComputedCharge().ByLepton );
			int top1charge = *(top1->GetComputedCharge().ByTVCM );
			if (top1charge * top2lepton > 0  && _stats._Top1gamma > gammacut1)
			{
				std::cout << "Vertex + lepton for top1 is used!\n";
				goodcharge.push_back(6);
				chargevalue.push_back(top1charge);
				_stats._methodCorrect = top1->__GetMCCharge() * top1charge < 0;
			}
			if (top1charge * top2lepton < 0 && _stats._Top1gamma > gammacut1) 
			{
				samecharge.push_back(6);
			}
		}
		if (top2->GetComputedCharge().ByLepton && top2->GetComputedCharge().ByTVCM)
		{
			int top2lepton = *(top2->GetComputedCharge().ByLepton );
			int top2charge = *(top2->GetComputedCharge().ByTVCM );
			if (top2charge * top2lepton < 0 &&  _stats._Top1gamma > gammacut1)
			{
				std::cout << "Vertex + lepton for top1 is used!\n";
				goodcharge.push_back(6);
				chargevalue.push_back(-top2charge);
				_stats._methodCorrect = top1->__GetMCCharge() * top2charge < 0;
			}
			if (top2charge * top2lepton > 0  && _stats._Top1gamma > gammacut1) 
			{
				samecharge.push_back(6);
			}
		}//
		//float chi2 = _stats._chiTopMass + _stats._chiTopE + _stats._chiGammaT + _stats._chiCosWb + _stats._chiPbstar;
		/*float chi2 =  _stats._chiGammaT + _stats._chiCosWb + _stats._chiPbstar;
		//if (top2->GetComputedCharge().ByLepton &&  _stats._Top1gamma > gammacut1+0.1  && goodcharge.size() == 0) 
		if (top2->GetComputedCharge().ByLepton &&  chi2 < 15) 
		{
			int top2lepton = *(top2->GetComputedCharge().ByLepton );
			//_stats._qCostheta[0] = (top2lepton < 0)? _stats._Top1costheta: -_stats._Top1costheta;
			goodcharge.push_back(7);
			chargevalue.push_back(top2lepton);
			_stats._methodCorrect = _stats._MCBWcorrect == 1;
			//_stats._methodUsed = 1;
			//_stats._methodTaken[0] = 7;
			//_summary._nAfterKinematicCuts++;
			//return;
		}//*/
		_stats._methodRefused = samecharge.size();
		_stats._methodUsed = goodcharge.size();
		if (samecharge.size() > 0) 
		{
			std::cout << "REFUSED BY CHARGE: ";
			for (unsigned int i = 0; i < samecharge.size(); i++) 
			{
				std::cout << " " << samecharge[i];
				_stats._methodSameCharge[i] = samecharge[i];
			}
			std::cout << "\n";
		}
		if (goodcharge.size() > 0)//(_stats._methodUsed > 0 && _stats._MCMass > _massCutparameter) //CRUNCH
		{
			std::cout << "ACCEPTED BY CHARGE: ";
			for (unsigned int i = 0; i < goodcharge.size(); i++) 
			{
				std::cout << " " << goodcharge[i];
				_stats._methodTaken[i] = goodcharge[i];
			}
			std::cout << "\n";
		}
		if (chargevalue.size() > 0) 
		{
			std::cout << "CHARGE VALUE: ";
			int sum = 0;
			for (unsigned int i = 0; i < chargevalue.size(); i++) 
			{
				std::cout << " " << chargevalue[i];
				sum += chargevalue[i];
			}
			std::cout << "\n";
			if (sum == 0) 
			{
				std::cout << "CONTRADICTING RESULT\n";
				_stats._qCostheta[0] = -2.;
				_stats._qCostheta[1] = -2.;
			}
			else 
			{
				//if ( _stats._MCMass > _massCutparameter) 
				{
					_summary._nAfterKinematicCuts++;
				}
				_stats._qCostheta[0] = (sum < 0)? _stats._Top1costheta: - _stats._Top1costheta;
				_stats._qCostheta[0] = (sum > 0)? -_stats._Top1costheta:  _stats._Top1costheta;
				_stats._qCostheta1 = (sum < 0)? _stats._Top1costheta: - _stats._Top1costheta;
				_stats._qCostheta1 = (sum > 0)? -_stats._Top1costheta:  _stats._Top1costheta;
			}
		}
		// Lepton only
		//_stats._gammaT = e/m;
		//if (_stats._methodUsed < 1 && samecharge > 0) 
		//{
			//_stats._methodSameCharge = 1;
		//}
		//_stats._methodUsed = 0;
		
	}
	void TTBarProcessor::AnalyseTTBarHadronic( LCEvent * evt )
	{
		try
		{
			std::cout << "***************Analysis*" <<_summary._nEvt++<<"*****************\n";
			LCCollection * jetcol = evt->getCollection(_JetsColName);
			LCCollection * jetrelcol = evt->getCollection(_JetsRelColName);
			LCCollection * mccol = evt->getCollection(_MCColName);
			LCCollection * mcvtxcol = evt->getCollection(_MCVtxColName);
			QQBarMCOperator opera(mccol);
			VertexChargeOperator vtxOperator(evt->getCollection(_colName),evt->getCollection(_colRelName));
			vector < MCParticle * > mctops = AnalyseGenerator(opera);
			vector< RecoJet * > * jets = getJets(jetcol, jetrelcol);
			//std::sort(jets->begin(), jets->end(), sortByBtag);
			vector< RecoJet * > * wjets = new vector< RecoJet * >();
			vector< RecoJet * > * bjets = getBTagJets(jets, wjets);
			if ( bjets->at(0)->GetBTag() < _highBTagCutparameter ||  bjets->at(1)->GetBTag() < _lowBTagCutparameter) 
			{
				return;	
			}
			_summary._nAfterBtagCuts++;
			vector< RecoJet * > * wbosons = formW(wjets);
			vector< TopQuark * > * tops = composeTops(bjets,wbosons);
			_hTree->Fill();
			ClearVariables();
		}
		catch(DataNotAvailableException &e)
		{
			std::cout << e.what() <<"\n";
		}
	}
	vector< TopQuark * > * TTBarProcessor::composeTops(vector< RecoJet * > * bjets, vector< RecoJet * > * wjets)
	{
		vector< TopQuark * > * result = new vector< TopQuark * > ();
		TopQuark * candidateb0w0 = new TopQuark(bjets->at(0), wjets->at(0));
		TopQuark * candidateb1w1 = new TopQuark(bjets->at(1), wjets->at(1));
		float chi00 = getChi2(candidateb0w0);
		float chi11 = getChi2(candidateb1w1);
		std::cout << "Chi2: " << chi00 << " " << chi11 << "\n";
		//vs
		TopQuark * candidateb0w1 = new TopQuark(bjets->at(0), wjets->at(1));
		TopQuark * candidateb1w0 = new TopQuark(bjets->at(1), wjets->at(0));
		float chi01 = getChi2(candidateb0w1);
		float chi10 = getChi2(candidateb1w0);
		if (chi00 + chi11 < chi01 + chi10) 
		{
			result->push_back(candidateb0w0);
			result->push_back(candidateb1w1);
			_stats._Top1mass = candidateb0w0->getMass();
			_stats._Top2mass = candidateb1w1->getMass();
		}
		else 
		{
			result->push_back(candidateb0w1);
			result->push_back(candidateb1w0);
			_stats._Top1mass = candidateb0w1->getMass();
			_stats._Top2mass = candidateb1w0->getMass();
		}
		std::cout << "Chi2: " << chi01 << " " << chi10 << "\n";
		return result;


	}
	vector< RecoJet * > * TTBarProcessor::formW(vector< RecoJet * > * wjets)
	{
		vector< RecoJet * > * result = new vector< RecoJet * > ();
		float bestdifference = 1000.;
		int icandidate = -1;
		int jcandidate = -1;
		int kcandidate = -1;
		int lcandidate = -1;
		for (unsigned int i = 0; i < wjets->size(); i++) 
		{
			for (unsigned int j = 0; j < i; j++) 
			{
				TopQuark * candidate1 = new TopQuark(wjets->at(i), wjets->at(j));
				float mass1 = candidate1->getMass();
				vector<int> ncandidates = getOpposite(i,j);
				TopQuark * candidate2 = new TopQuark(wjets->at(ncandidates[0]), wjets->at(ncandidates[1]));
				float mass2 = candidate2->getMass();
				std::cout << "W candidate masses: " << mass1 << " " << mass2 << "\n";
				//float ctag1 = wjets->at(i)->GetCTag();
				//float ctag2 = wjets->at(j)->GetCTag();
				//float ctag3 = wjets->at(ncandidates[0])->GetCTag();
				//float ctag4 = wjets->at(ncandidates[1])->GetCTag();
				std::cout << "W candidate c-tags: " << wjets->at(i)->GetCTag() << " " << wjets->at(j)->GetCTag() << "\n";
				if (std::abs(mass1 + mass2 - 2*_WMassparameter) < bestdifference)// && 
						//(!(ctag1 > _lowBTagCutparameter &&  ctag2 > _lowBTagCutparameter) &&
						//!(ctag3 > _lowBTagCutparameter && ctag4 > _lowBTagCutparameter) 
						//|| ctag1+ctag2+ctag3 > 3*_lowBTagCutparameter)) 
				{
					bestdifference = std::abs(mass1+mass2 - 2*_WMassparameter);
					icandidate = i;
					jcandidate = j;
					kcandidate = ncandidates[0];
					lcandidate = ncandidates[1];
				}

			}
		}
		result->push_back(new TopQuark(wjets->at(icandidate), wjets->at(jcandidate)));
		_stats._W1mass = result->at(0)->getMass();
		std::cout << "i: " << icandidate << " j: " << jcandidate  << " k: " << kcandidate  << " l: " << lcandidate <<"\n";
		result->push_back(new TopQuark(wjets->at(kcandidate), wjets->at(lcandidate)));
		_stats._W2mass = result->at(1)->getMass();
		return result;
	}
	void TTBarProcessor::AnalyseBBBar(LCEvent * evt)
	{
		LCCollection * mcvtxcol = NULL;
		try
		{
			mcvtxcol = evt->getCollection(_MCVtxColName);
		}
		catch(DataNotAvailableException &e)
		{
			std::cout << e.what() <<"\n";
		}
		try
		{
			std::cout << "***************Analysis*" <<_summary._nEvt++<<"*****************\n";
			LCCollection * jetcol = evt->getCollection(_JetsColName);
			LCCollection * jetrelcol = evt->getCollection(_JetsRelColName);
			LCCollection * mccol = evt->getCollection(_MCColName);
			//LCCollection * mcvtxcol = evt->getCollection(_MCVtxColName);
			LCCollection * pfocol = evt->getCollection(_colName);
			QQBarMCOperator opera(mccol);
			vector < MCParticle * > mcbs = AnalyseGeneratorBBBar(opera);
			VertexChargeOperator vtxOperator(evt->getCollection(_colName),evt->getCollection(_colRelName));
			vector< RecoJet * > * jets = getJets(jetcol, jetrelcol);
			std::sort(jets->begin(), jets->end(), sortByBtag);
			if (jets->size() < 2 || (jets->at(0)->GetBTag() < _highBTagCutparameter || jets->at(1)->GetBTag() < _lowBTagCutparameter)) 
			//if (jets->size() < 2) 
			{
				return;
			}
			std::cout << "MCPDG: " << _stats._MCPDG << "\n";
			if (_stats._MCMass > _massCutparameter && _stats._MCPDG == 5) 
			{
				_summary._nAfterBtagCuts++;
			}
			std::cout << "--EVENT ACCEPTED--\n";
			_stats._B1momentum = jets->at(0)->GetHadronMomentum(); //MathOperator::getModule(jets->at(0)->getMomentum());
			_stats._B2momentum = jets->at(1)->GetHadronMomentum();//MathOperator::getModule(jets->at(1)->getMomentum());
			_stats._B1Jetmomentum = MathOperator::getModule(jets->at(0)->getMomentum());
			_stats._B2Jetmomentum = MathOperator::getModule(jets->at(1)->getMomentum());
			_stats._B1btag = jets->at(0)->GetBTag();
			_stats._B2btag = jets->at(1)->GetBTag();
			_stats._B1mass = jets->at(0)->getMass();
			_stats._B2mass = jets->at(1)->getMass();
			_stats._B1charge =  jets->at(0)->GetHadronCharge();
			_stats._B2charge =  jets->at(1)->GetHadronCharge();
			_stats._bbbarAngle = MathOperator::getAngleBtw(jets->at(0)->getMomentum(), jets->at(1)->getMomentum());
			RecoJet * Zboson = new TopQuark(jets->at(0), jets->at(1));
			_stats._InvMass = Zboson->getMass();
			findPhoton(pfocol, jets);
			if (_stats._InvMass > _massCutparameter-20 && _stats._maxPhotonEnergy < 40 && _stats._MCPDG == 5) //CRUNCH
			{
				_summary._nAfterMassCuts++;
			}
			vector< const double * > momentums;// = {jets->at(0)->getMomentum(), jets->at(1)->getMomentum()};
			momentums.push_back(jets->at(0)->getMomentum());
			momentums.push_back(jets->at(1)->getMomentum());
			vector< ReconstructedParticle * > kaons1 = vtxOperator.GetKaons(jets->at(0));
			if (kaons1.size() > 0) 
			{
				_stats._kaonMomentum = MathOperator::getModule(kaons1[0]->getMomentum());
			}
			 _stats._bbbarPt = 0.0;
			 _stats._bbbarP = 0.0;
			for (unsigned int i = 0; i < 3; i++) 
			{
				 _stats._bbbarP += std::pow(jets->at(0)->getMomentum()[i] + jets->at(1)->getMomentum()[i],2);
				 _stats._bbbarPt += (i == 2)? 0.0:std::pow(jets->at(0)->getMomentum()[i] + jets->at(1)->getMomentum()[i],2);
			}
			_stats._bbbarPt = std::sqrt(_stats._bbbarPt);
			_stats._bbbarP = std::sqrt(_stats._bbbarP);
			std::cout << "Missing pt: " << _stats._bbbarPt << "\n";
			vector<float> thrust;
			pfocol->getParameters().getFloatVals ("majorThrustAxis", thrust);
			//_stats._ThrustCos = std::cos(MathOperator::getAngles(thrust)[1]); 
			_stats._SumCos = -2; 
			MatchB(jets, mcbs, mcvtxcol);
			ComputeCharge(jets, vtxOperator);
			ClusteringOperator cloperator;
			if (jets->at(0)->GetRecoVertices() && jets->at(0)->GetRecoVertices()->size() > 0) 
			{
				_stats._B1Y = MathOperator::getAngleBtw(jets->at(0)->getMomentum(), jets->at(0)->GetRecoVertices()->at(0)->getAssociatedParticle()->getMomentum()); //cloperator.ClusterizeJet(jets->at(0));
			}
			//_stats._B2Y = cloperator.ClusterizeJet(jets->at(1));
			_stats._B1chargeBalance = getChargeBalance(jets->at(0));
			_stats._B2chargeBalance = getChargeBalance(jets->at(1));
			_stats._Sphericity = pfocol->getParameters().getFloatVal("sphericity");
			_stats._Thrust = pfocol->getParameters().getFloatVal("majorThrustValue");
			//getZZ(evt->getCollection("FourJets"));
			//_stats._bbbar4JetMass = pfocol->getNumberOfElements();
			vector<float> axis;
			//pfocol->getParameters().getFloatVals("principleThrustAxis", axis);
			pfocol->getParameters().getFloatVals("majorThrustAxis", axis);
			//getThrust(axis, pfocol);
			_hTree->Fill();
			ClearVariables();
			delete Zboson;
			for (unsigned int i = jets->size(); i > -1; i--) 
			{
				RecoJet * jet = jets->at(i);
				jets->pop_back();
				delete jet;
			}
			delete jets;
		}
		catch(DataNotAvailableException &e)
		{
			std::cout << e.what() <<"\n";
		}
		
	}
	void TTBarProcessor::ComputeCharge(std::vector< RecoJet * > *jets, VertexChargeOperator & vtxOperator)
	{
		//std::sort(jets->begin(), jets->end(), sortByCostheta);
		PrintJets(jets);
		if (jets->at(1)->__GetMCNtracks() != jets->at(1)->GetNumberOfVertexParticles()) 
		{
			//return;	
		}
		/*if (jets->at(0)->GetNumberOfVertexParticles() > 2) 
		{
			
		}
		JetCharge & topCharge = top2->GetComputedCharge();
		topCharge.ByTrackCount = new int(top2->GetHadronCharge());
*/	
		vector<float> bdirection1 = MathOperator::getDirection(jets->at(0)->getMomentum());
		vector<float> bdirection2 = MathOperator::getDirection(jets->at(1)->getMomentum());
		_stats._B1costheta = jets->at(0)->GetCostheta();//( std::cos( MathOperator::getAngles(bdirection1)[1] ));
		_stats._B2costheta = jets->at(1)->GetCostheta();//( std::cos( MathOperator::getAngles(bdirection2)[1] ));
		//_stats._B1costheta = ( std::cos( MathOperator::getAngles(bdirection1)[1] ));
		//_stats._B2costheta = ( std::cos( MathOperator::getAngles(bdirection2)[1] ));
		for (unsigned int i = 0; i < jets->size(); i++) 
		{
			if (jets->at(i)->GetNumberOfVertexParticles() > 2)// && abs(jets->at(i)->GetHadronCharge()) > 0 ) 
			{
				JetCharge & bCharge = jets->at(i)->GetComputedCharge();
				int chargeValue = (jets->at(i)->GetHadronCharge() > -5)? jets->at(i)->GetHadronCharge() : 0;
				bCharge.ByTrackCount = (chargeValue == 0)? new int(0): new int((float)chargeValue / abs(chargeValue));
				
			}
			if (jets->at(i)->GetNumberOfVertexParticles() > 0) 
			{
				vtxOperator.ComputeCharge(jets->at(i));
			}
			vector<float> bdirection = MathOperator::getDirection(jets->at(i)->getMomentum());
			float costheta =  std::cos( MathOperator::getAngles(bdirection)[1] );
		}
		float btagcut = .80;
		float pcut = 20.;
		RecoJet * top1 = jets->at(0);
		RecoJet * top2 = jets->at(1);
		std::cout << "Jet \t VTX \t KAON\n";
		std::cout << "B1:\t" << intToStr(top1->GetComputedCharge().ByTrackCount) <<"\t"
					<< intToStr(top1->GetComputedCharge().ByTVCM) <<"\n";
					//<< intToStr(top1->GetComputedCharge().ByLepton) <<"\n";
		std::cout << "B2:\t" << intToStr(top2->GetComputedCharge().ByTrackCount) <<"\t"
					<< intToStr(top2->GetComputedCharge().ByTVCM) <<"\n";
					//<< intToStr(top2->GetComputedCharge().ByLepton) <<"\n";
		if ((top1->GetComputedCharge().ByTVCM || top2->GetComputedCharge().ByTVCM) && _stats._MCMass > _massCutparameter && _stats._MCPDG ==5) 
		{
			_summary._nKaons++;
		}
		if ((top1->GetComputedCharge().ByTrackCount || top2->GetComputedCharge().ByTrackCount) && _stats._MCMass > _massCutparameter && _stats._MCPDG == 5) 
		{
			_summary._nChargedB++;
		}
		_stats._qCostheta[0] = -2.;
		_stats._qCostheta[1] = -2.;
		_stats._B1VtxTag = (top1->GetComputedCharge().ByTrackCount)? *(top1->GetComputedCharge().ByTrackCount) : -5;
		_stats._B2VtxTag = (top2->GetComputedCharge().ByTrackCount)? *(top2->GetComputedCharge().ByTrackCount) : -5;
		_stats._B1KaonTag = (top1->GetComputedCharge().ByTVCM)? *(top1->GetComputedCharge().ByTVCM) : -5;
		_stats._B2KaonTag = (top2->GetComputedCharge().ByTVCM)? *(top2->GetComputedCharge().ByTVCM) : -5;
		vector<int> samecharge;
		vector<int> goodcharge;
		vector<int> zerocharge;
		vector<int> chargevalue;
		///BBB
		//Track charge * Track charge
		if (top2->GetComputedCharge().ByTrackCount && top1->GetComputedCharge().ByTrackCount) 
		{
			int top1charge = *(top1->GetComputedCharge().ByTrackCount );
			int top2charge = *(top2->GetComputedCharge().ByTrackCount );
			if (top1charge * top2charge < 0 && _stats._B2btag > btagcut && _stats._B2momentum > pcut)//  && _stats._methodUsed < 1) 
			{
				//_stats._qCostheta[0] = (top1charge < 0)? _stats._B1costheta: - _stats._B1costheta;
				//_stats._qCostheta[1] = (top1charge > 0)? _stats._B2costheta: - _stats._B2costheta;
				chargevalue.push_back(top1charge);
				std::cout << "Two vertices are used!\n";
				//_stats._methodUsed = 1;
				goodcharge.push_back(1);
				_stats._methodCorrect = top1->__GetMCCharge() * top1charge > 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				//_summary._nAfterKinematicCuts++;
				//return;
			}
			if ( _stats._B2btag > btagcut && _stats._B2momentum > pcut && (top1charge == 0 || top2charge == 0)) 
			{
				zerocharge.push_back(1);
			}
			if (top1charge * top2charge > 0 && _stats._B2btag > btagcut && _stats._B2momentum > pcut) 
			{
				samecharge.push_back(1);
			}
		}
		
		//Kaon * Kaon
		if (top2->GetComputedCharge().ByTVCM && top1->GetComputedCharge().ByTVCM) 
		{
			int top1charge = *(top1->GetComputedCharge().ByTVCM );
			int top2charge = *(top2->GetComputedCharge().ByTVCM );
			bool zero = false;
			if (top1charge * top2charge < 0)//  && _stats._methodUsed < 1) 
			{
				//_stats._qCostheta[0] = (top1charge < 0)? _stats._B1costheta: - _stats._B1costheta;
				//_stats._qCostheta[1] = (top1charge > 0)? _stats._B2costheta: - _stats._B2costheta;
				std::cout << "Two kaons are used!\n";
				//_stats._methodUsed = 2;
				goodcharge.push_back(2);
				chargevalue.push_back(top1charge);
				_stats._methodCorrect = top1->__GetMCCharge() * top1charge > 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				//_summary._nAfterKinematicCuts++;
				//return;
			}
			if (top1->GetComputedCharge().ByTrackCount && *(top1->GetComputedCharge().ByTrackCount ) == 0) 
			{
				zerocharge.push_back(2);
				zero = true;
			}
			if (top2->GetComputedCharge().ByTrackCount && *(top2->GetComputedCharge().ByTrackCount ) == 0 && !zero) 
			{
				zerocharge.push_back(2);
			}
			if (top1charge * top2charge > 0) 
			{
				samecharge.push_back(2);
			}
		}////
		//Track charge + Kaon
		/*if (top1->GetComputedCharge().ByTrackCount && top1->GetComputedCharge().ByTVCM) 
		{
			int top1charge = *(top1->GetComputedCharge().ByTrackCount );
			int top1kaon = *(top1->GetComputedCharge().ByTVCM );
			if (top1charge * top1kaon > 0 && _stats._B1btag > btagcut && _stats._B1momentum > pcut) 
			{
				//_stats._qCostheta[0] = (top1charge < 0)? _stats._B1costheta: -_stats._B1costheta; 
				//_stats._qCostheta[1] = (top1charge > 0)? _stats._B2costheta: - _stats._B2costheta;
				std::cout << "Vertex + kaon for B1 is used!\n";
				//_stats._methodUsed = 3;
				goodcharge.push_back(3);
				chargevalue.push_back(top1charge);
				_stats._methodCorrect = top1->__GetMCCharge() * top1charge > 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
					 std::cout << "True VTX: " <<  top1->__GetMCCharge() << "\n";
				}
				//_summary._nAfterKinematicCuts++;
				//return;
			}
			if ( _stats._B1btag > btagcut && _stats._B1momentum > pcut && (top1charge == 0)) 
			{
				zerocharge.push_back(3);
			}
			if (top1charge * top1kaon < 0 &&_stats._B1btag > btagcut && _stats._B1momentum > pcut) 
			{
				samecharge.push_back(3);
			}

		}
		if (top2->GetComputedCharge().ByTrackCount && top2->GetComputedCharge().ByTVCM) 
		{
			int top2charge = *(top2->GetComputedCharge().ByTrackCount );
			int top2kaon = *(top2->GetComputedCharge().ByTVCM );
			if (top2charge * top2kaon > 0 && _stats._B2btag > btagcut && _stats._B2momentum > pcut ) 
			{
				//_stats._qCostheta[0] = (top2charge > 0)? _stats._B1costheta: - _stats._B1costheta; 
				//_stats._qCostheta[1] = (top2charge < 0)? _stats._B2costheta: - _stats._B2costheta;
				std::cout << "Vertex + kaon for B2 is used!\n";
				//_stats._methodUsed = 4;
				chargevalue.push_back(-top2charge);
				goodcharge.push_back(4);
				_stats._methodCorrect = top1->__GetMCCharge() * top2charge < 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
					 std::cout << "True VTX: " <<  top2->__GetMCCharge() << "\n";
				}
				//_summary._nAfterKinematicCuts++;
				//return;
			}
			if ( _stats._B2btag > btagcut && _stats._B2momentum > pcut && (top2charge == 0)) 
			{
				zerocharge.push_back(4);
			}
			if (top2charge * top2kaon < 0 && _stats._B2btag > btagcut && _stats._B2momentum > pcut) 
			{
				samecharge.push_back(4);
			}
		}//
		if (top2->GetComputedCharge().ByTrackCount && top1->GetComputedCharge().ByTVCM) 
		{
			int top2charge = *(top2->GetComputedCharge().ByTrackCount );
			int top1kaon = *(top1->GetComputedCharge().ByTVCM );
			if (top2charge * top1kaon < 0 &&  _stats._B2btag > btagcut && _stats._B2momentum > pcut) 
			{
				//_stats._qCostheta[0] = (top2charge > 0)? _stats._B1costheta: - _stats._B1costheta; 
				//_stats._qCostheta[1] = (top2charge < 0)? _stats._B2costheta: - _stats._B2costheta;
				std::cout << "Vertex2 + kaon1 is used!\n";
				//_stats._methodUsed = 6;
				chargevalue.push_back(-top2charge);
				goodcharge.push_back(6);
				_stats._methodCorrect = top1->__GetMCCharge() * top1kaon > 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				//_summary._nAfterKinematicCuts++;
				//return;
			}
			if ( _stats._B2btag > btagcut && _stats._B2momentum > pcut && top1->GetComputedCharge().ByTrackCount && *(top1->GetComputedCharge().ByTrackCount ) == 0) 
			{
				zerocharge.push_back(6);
			}
			if (top2charge * top1kaon > 0 &&  _stats._B2btag > btagcut && _stats._B2momentum > pcut)
			{
				samecharge.push_back(6);
			}
		}///
		if (top1->GetComputedCharge().ByTrackCount && top2->GetComputedCharge().ByTVCM) 
		{
			int top1charge = *(top1->GetComputedCharge().ByTrackCount );
			int top2kaon = *(top2->GetComputedCharge().ByTVCM );
			if (top1charge * top2kaon < 0 &&  _stats._B1btag > btagcut && _stats._B1momentum > pcut) 
			{
				//_stats._qCostheta[0] = (top1charge < 0)? _stats._B1costheta: - _stats._B1costheta; 
				//_stats._qCostheta[1] = (top1charge > 0)? _stats._B2costheta: - _stats._B2costheta;
				std::cout << "Vertex1 + kaon2 is used!\n";
				//_stats._methodUsed = 5;
				chargevalue.push_back(top1charge);
				goodcharge.push_back(5);
				_stats._methodCorrect = top1->__GetMCCharge() *  top1charge > 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				//_summary._nAfterKinematicCuts++;
				//return;
			}
			if ( _stats._B1btag > btagcut && _stats._B1momentum > pcut && top2->GetComputedCharge().ByTrackCount &&  *(top2->GetComputedCharge().ByTrackCount ) == 0) 
			{
				zerocharge.push_back(5);
			}
			if (top1charge * top2kaon > 0 &&_stats._B1btag > btagcut && _stats._B1momentum > pcut) 
			{
				samecharge.push_back(5);
			}

		}//*/
		_stats._methodRefused = samecharge.size();
		_stats._methodUsed = goodcharge.size();
		_stats._methodZero = zerocharge.size();
		if (samecharge.size() > 0) 
		{
			std::cout << "REFUSED BY CHARGE: ";
			for (unsigned int i = 0; i < samecharge.size(); i++) 
			{
				std::cout << " " << samecharge[i];
				_stats._methodSameCharge[i] = samecharge[i];
			}
			std::cout << "\n";
		}
		if (zerocharge.size() > 0) 
		{
			std::cout << "ZERO CHARGE: ";
			for (unsigned int i = 0; i < zerocharge.size(); i++) 
			{
				std::cout << " " << zerocharge[i];
				_stats._methodZeroCharge[i] = zerocharge[i];
			}
			std::cout << "\n";
		}
		if (goodcharge.size() > 0)//(_stats._methodUsed > 0 && _stats._MCMass > _massCutparameter) //CRUNCH
		{
			std::cout << "ACCEPTED BY CHARGE: ";
			for (unsigned int i = 0; i < goodcharge.size(); i++) 
			{
				std::cout << " " << goodcharge[i];
				_stats._methodTaken[i] = goodcharge[i];
			}
			std::cout << "\n";
		}
		if (chargevalue.size() > 0) 
		{
			std::cout << "CHARGE VALUE: ";
			int sum = 0;
			for (unsigned int i = 0; i < chargevalue.size(); i++) 
			{
				std::cout << " " << chargevalue[i];
				sum += chargevalue[i];
			}
			std::cout << "\n";
			if (sum == 0) 
			{
				std::cout << "CONTRADICTING RESULT\n";
				_stats._qCostheta[0] = -2.;
				_stats._qCostheta[1] = -2.;
				_stats._methodUsed = 0;
			}
			else 
			{
				if ( _stats._MCMass > _massCutparameter && _stats._MCPDG == 5) 
				{
					_summary._nAfterKinematicCuts++;
				}
				double momentum[3]; 
				for (unsigned int i = 0; i < 3; i++) 
				{
					momentum[i] =  jets->at(0)->getMomentum()[i] - jets->at(1)->getMomentum()[i];
				}
				vector <float> direction = MathOperator::getDirection(momentum);
				float costheta = std::cos(MathOperator::getAngles(direction)[1]);
				_stats._ThrustCos = (sum < 0)?costheta:-costheta;
				_stats._SumCos = (sum < 0)?costheta:-costheta;
				_stats._qCostheta[0] = (sum < 0)? _stats._B1costheta: - _stats._B1costheta;
				_stats._qCostheta[1] = (sum > 0)? _stats._B2costheta: - _stats._B2costheta;
				_stats._qCostheta1 = (sum < 0)? _stats._B1costheta: - _stats._B1costheta;
				_stats._qCostheta2 = (sum > 0)? _stats._B2costheta: - _stats._B2costheta;
			}
		}
		return;//*/
		/*if (_stats._B1charge * _stats._B2charge < 0 && _stats._B1momentum > 0 && _stats._B2momentum > 0) 
		{
			int b1charge = _stats._B1charge  / std::abs(_stats._B1charge);
			int b2charge = _stats._B2charge  / std::abs(_stats._B2charge);
			_stats._qCostheta[0] = -_stats._B1costheta * b1charge;
			_stats._qCostheta[1] = -_stats._B2costheta * b2charge;
		}*/
		
	}
	ReconstructedParticle * TTBarProcessor::findPhoton(LCCollection * pfocol, vector< RecoJet * > *jets)
	{
		int pfonumber = pfocol->getNumberOfElements();
		ReconstructedParticle * winner = NULL;
		float maxenergy = 0.0;
		for (int i = 0; i < pfonumber; i++) 
		{
			ReconstructedParticle * particle = dynamic_cast< ReconstructedParticle * >(pfocol->getElementAt(i));
			if ((particle->getType() == 22 || particle->getType() ==2112) && particle->getEnergy() > maxenergy) 
			{
				winner = particle;
				maxenergy = particle->getEnergy();
			}
		}
		if (winner) 
		{
			vector<float> gdirection = MathOperator::getDirection(winner->getMomentum());
			_stats._maxPhotonCostheta = (std::cos( MathOperator::getAngles(gdirection)[1] ));
			_stats._maxPhotonEnergy = winner->getEnergy();
			if(!jets)
			{
				return winner;
			}
			float angle1 = MathOperator::getAngleBtw(winner->getMomentum(), jets->at(0)->getMomentum());
			float angle2 = MathOperator::getAngleBtw(winner->getMomentum(), jets->at(1)->getMomentum());
			_stats._maxPhotonAngle = (angle1 < angle2)? angle1 : angle2;
		}
		return winner;
	}

	float TTBarProcessor::getChi2(TopQuark * candidate)
	{
		float pT = MathOperator::getModule(candidate->getMomentum());
		float mT = candidate->getMass();
		float ET = candidate->getEnergy();
		float pB = MathOperator::getModule(candidate->GetB()->getMomentum());
		float beta = pT / ET;
		float gamma =1.0 / std::sqrt( 1.0 - beta * beta);
		float cosbT = std::cos( MathOperator::getAngle(candidate->getMomentum(), candidate->GetB()->getMomentum()) );
		float bpstar = gamma * pB * ( 1.0 - beta * cosbT);
		float cosbW = std::cos( MathOperator::getAngle(candidate->GetB()->getMomentum(), candidate->GetW()->getMomentum()) );
		_stats._Top1pstarb = bpstar;

		_stats._chiTopMass = std::pow(mT - _TopMassparameter, 2) / std::pow( _TopMassSigmaparameter, 2) ;
		_stats._chiTopE = std::pow(ET - _EBeamparameter , 2) / std::pow( _EBeamSigmaparameter, 2); 
		_stats._chiPbstar = std::pow( bpstar - _PStarparameter, 2) / std::pow( _PStarSigmaparameter, 2);
		_stats._chiGammaT = std::pow( gamma - _GammaTparameter, 2) / std::pow( _GammaTSigmaparameter, 2); 
		_stats._chiCosWb = std::pow( cosbW - _CosbWparameter, 2) / std::pow( _CosbWSigmaparameter, 2);
		float chi2 = _stats._chiTopMass + _stats._chiTopE + _stats._chiPbstar  + _stats._chiCosWb + _stats._chiGammaT;
		//float chi2 = std::pow( mT - _stats._TopMassparameter, 2) / std::pow( _stats._TopMassSigmaparameter, 2) + 
		//	     std::pow(ET - _EBeamparameter , 2) / std::pow( _EBeamSigmaparameter, 2) +
		//	     std::pow( bpstar - _PStarparameter, 2) / std::pow( _PStarSigmaparameter, 2) +
		//	     std::pow( gamma - _GammaTparameter, 2) / std::pow( _GammaTSigmaparameter, 2) +
		//	     std::pow( cosbW - _CosbWparameter, 2) / std::pow( _CosbWSigmaparameter, 2);
		std::cout << " chi2: " << chi2
			  << " bpstar: " << bpstar
			  << " cosw: " << cosbW 
			  << "\n";
		return chi2;
	}
	vector< RecoJet * > * TTBarProcessor::getBTagJets(std::vector< RecoJet * > * alljets, std::vector< RecoJet * > * wjets)
	{
		vector< RecoJet * > * result = new vector< RecoJet * >();
		if (!wjets) 
		{
			std::cout << "Wjets are NULL!\n";
			return result;
		}
		if (alljets->size() < 4) 
		{
			return result;
		}
		//vector<RecoJet>
		std::cout << "Before sorting:\n";
		for (unsigned int i = 0; i < alljets->size(); i++) 
		{
			RecoJet * jet = alljets->at(i);
			std::cout << "\tB tag: " << jet->GetBTag()  << "\n";
			/*if (jet->GetBTag() > _lowBTagCutparameter) 
			{
				result->push_back(jet);
			}
			else 
			{
				wjets->push_back(jet);
			}*/
		}
		std::sort(alljets->begin(), alljets->end(), sortByBtag);
		std::cout << "After sorting:\n";
		for (unsigned int i = 0; i < alljets->size(); i++) 
		{
			RecoJet * jet = alljets->at(i);
			std::cout << "\tB tag: " << jet->GetBTag()  << "\n";

		}
		result->push_back(alljets->at(0));
		result->push_back(alljets->at(1));
		for (unsigned int i = 2; i < alljets->size(); i++) 
		{
			wjets->push_back(alljets->at(i));
		}
		return result;
	}
	vector< RecoJet * > * TTBarProcessor::getJets(LCCollection * jetcol, LCCollection *jetrelcol)
	{
		int jetnumber = jetcol->getNumberOfElements();
		vector< RecoJet * > * result = new vector< RecoJet * >();
		LCRelationNavigator navigator(jetrelcol);
		PIDHandler pidh(jetcol);
		int alid = -1;
		try
		{
			alid = pidh.getAlgorithmID("vtxrec");
		}
		catch(UTIL::UnknownAlgorithm &e)
		{
			std::cout << "No algorithm vtxrec!\n";
			alid = -1;
		}
		if (alid < 0) 
		{
			try
			{
				alid = pidh.getAlgorithmID("lcfiplus");
			}
			catch(UTIL::UnknownAlgorithm &e)
			{
				std::cout << "No algorithm lcfiplus!\n";
				alid = -1;
			}
			
		}
		std::cout << "Algorithm id: " << jetnumber << "\n";
		for (int j = 0; j < jetnumber; j++) 
		{
			ReconstructedParticle * jetpart = dynamic_cast< ReconstructedParticle * >(jetcol->getElementAt(j));
			vector< Vertex * > * vertices = convert(navigator.getRelatedToObjects(jetpart));
			const vector< ReconstructedParticle * > components = jetpart->getParticles();
			int nvtx = vertices->size();
			float btag = 0.0;
			float ctag = 0.0;
			if (alid > -1) 
			{
				const ParticleID& pid = pidh.getParticleID(jetpart,alid);
				vector<float> params = pid.getParameters();
				btag = params[pidh.getParameterIndex(alid,"BTag")];
				ctag = params[pidh.getParameterIndex(alid,"CTag")];
			}
			RecoJet * jet = new RecoJet(jetpart, btag, ctag, nvtx);
			jet->SetRecoVertices(vertices);
			PrintJet(jet);
			result->push_back(jet);
		}
		return result;
		
	}
	vector< Vertex * > * TTBarProcessor::convert(const std::vector< LCObject * > & objs)
	{
		std::vector< Vertex * > * result = new std::vector< Vertex * >();
		for (unsigned int i = 0; i < objs.size(); i++) 
		{
			result->push_back(dynamic_cast< Vertex * >(objs[i]));
		}
		return result;
	}
	void TTBarProcessor::PrintJet(RecoJet * jet)
	{
		std::cout << "Jet E: " << jet->getEnergy()
			  << " m: " << jet->getMass()
			  << " btag: " << jet->GetBTag()
			  << " costheta: " << jet->GetCostheta()
			  << " p_B: " << jet->GetHadronMomentum()
			  << " ntracks: " << jet->GetNumberOfVertexParticles()
			  << "\n";
	}
	void TTBarProcessor::PrintParticle(ReconstructedParticle * jet)
	{
		 std::cout << "E: " << jet->getEnergy()
		 	   <<" m: " << jet->getMass()
			   <<" PDG: " << jet->getType()
			   << "\n";
	}
	void TTBarProcessor::PrintParticle(MCParticle * jet)
	{
		 std::cout << "E: " << jet->getEnergy()
		 	   <<" m: " << jet->getMass()
		 	   <<" q: " << jet->getCharge()
			   <<" PDG: " << jet->getPDG()
			   << "\n";
	}
	void TTBarProcessor::PrintJets(std::vector< RecoJet * > *jets)
	{
		for (unsigned int i = 0; i < jets->size(); i++) 
		{
			PrintJet(jets->at(i));
		}
	}
	void TTBarProcessor::__ComputeChargeCheat(TopQuark * top1, TopQuark * top2)
	{
		//_Top1bcharge = (top1->GetB()->__GetMCCharge() > 0.0)? -1:1;
		if (top1->GetB()->__GetMCCharge() > 0.0) 
		{
			_stats._Top1bcharge = -1.0;
		}
		if (top1->GetB()->__GetMCCharge() < 0.0) 
		{
			_stats._Top1bcharge = 1.0;
		}
		//_Top1bcharge = top1->GetHadronCharge();
		//_Top2bcharge = top2->GetHadronCharge();

		_stats._Top1bmomentum = top1->GetHadronMomentum();
		vector<float> direction = MathOperator::getDirection(top1->getMomentum());
		_stats._Top1costheta =  std::cos( MathOperator::getAngles(direction)[1] );
		_stats._Top1bntracks = top1->GetNumberOfVertexParticles();
		int Top1Genbntracks = top1->GetB()->__GetMCNtracks();
		_stats._Top1btag = top1->GetBTag();
		std::cout << "Top charge: " << _stats._Top1bcharge
			  << " top pB: " << _stats._Top1bmomentum
			  << " top bntracks: " << _stats._Top1bntracks
			  << " top Gen bntracks: " << Top1Genbntracks 
			  << " top btag: " << _stats._Top1btag
			  << " top costheta: " << _stats._Top1costheta
			  << "\n";
		//Top Leptonic
		_stats._Top2bcharge = top2->__GetMCCharge();
		_stats._Top2bmomentum = top2->GetHadronMomentum();
		vector<float> direction2 = MathOperator::getDirection(top2->getMomentum());
		_stats._Top2costheta =  std::cos( MathOperator::getAngles(direction2)[1] );
		_stats._Top2bntracks = top2->GetNumberOfVertexParticles();
		int Top2Genbntracks = top2->GetB()->__GetMCNtracks();
		_stats._Top2btag = top2->GetBTag();
		std::cout << "Top charge: " << _stats._Top2bcharge
			  << " top pB: " << _stats._Top2bmomentum
			  << " top bntracks: " << _stats._Top2bntracks
			  << " top Gen bntracks: " << Top2Genbntracks 
			  << " top btag: " << _stats._Top2btag
			  << " top costheta: " << _stats._Top2costheta
			  << "\n";
		/*if (Top1Genbntracks == _stats._Top1bntracks && abs( _stats._Top1bcharge) > 0) 
		{
			_stats._qCostheta[0] = (_stats._Top1bcharge < 0)? _stats._Top1costheta : -_stats._Top1costheta;
		}
		if (Top2Genbntracks == _stats._Top2bntracks  && abs( _stats._Top2bcharge) > 0) 
		{
			_stats._qCostheta[0] = (_stats._Top2bcharge > 0)? _stats._Top1costheta : -_stats._Top1costheta;
		}*/
		_stats._qCostheta[0] = _stats._Top1costheta * _stats._Top1bcharge;
		
	}
	void TTBarProcessor::Match(vector< MCParticle * > & mctops,vector< MCParticle * > &  mcbs, vector< MCParticle * > & mcws, TopQuark * topHadronic,  TopQuark * top2)
	{
		_stats._Top1truthAngle = 4.0;
		float charge = 0.0;
		float minwangle = 4.0;
		float minbangle = 4.0;
		MCParticle * whad = NULL;
		MCParticle * bhad = NULL;
		for (unsigned int i = 0; i < mcbs.size(); i++) 
		{
			float angleB =  MathOperator::getAngle(mcbs[i]->getMomentum(), topHadronic->GetB()->getMomentum());
			float angleW =  MathOperator::getAngle(mcws[i]->getMomentum(), topHadronic->GetW()->getMomentum());
				std::cout << "\tB angle: " << angleB << "\n";
				std::cout << "\tW angle: " << angleW << "\n";
			if (angleW < minwangle) 
			{
				whad = mcws[i];
				minwangle = angleW;
			}
			if (angleB < minbangle) 
			{
				bhad = mcbs[i];
				minbangle = angleB;
			}
		}
		std::cout << "Min w angle: " << minwangle << " min b angle: "  << minbangle << "\n";
		MCParticle * tophad = NULL;
		for (unsigned int i = 0; i < mctops.size(); i++) 
		{
			MCParticle * mctop = mctops[i];
			float angle =  MathOperator::getAngle(mctop->getMomentum(), topHadronic->getMomentum());
			if (angle < _stats._Top1truthAngle) 
			{
				_stats._Top1truthAngle = angle;
				charge = mctop->getCharge();
				topHadronic->__SetMCCharge(charge);
				tophad = mctop;
			}
		}
		if (whad && bhad) 
		{
			std::cout << "MC charge: " << whad->getCharge() + bhad->getCharge() << "\n";
			std::cout << "W charge: " << whad->getCharge() <<  " w had: " << whad->getGeneratorStatus() << "\n";
			_stats._MCBWcorrect = (std::abs(whad->getCharge() + bhad->getCharge()) < 1 );
		}
		else 
		{
			_stats._MCBWcorrect = -1;
		}
		std::cout << "Truth Angle: " << _stats._Top1truthAngle << " charge: " << charge << "\n";
		if (top2) 
		{
			top2->__SetMCCharge(-charge);
		}
	}
	void TTBarProcessor::MatchB(vector<RecoJet*> * bjets, vector< MCParticle * > & mcbs, LCCollection * mcvtxcol)
	{
		float charge = 0.0;
		float angle = 4.0;
		int bgenntracks = 0;
		int bbargenntracks = 0;
		_stats._MCBOscillation = 0;
		_stats._MCBBarOscillation = 0;
		if (mcvtxcol && mcvtxcol->getNumberOfElements() > 0) 
		{
			for (int i = 0; i < mcvtxcol->getNumberOfElements(); i++) 
			{
				Vertex * vertex = dynamic_cast< Vertex * >(mcvtxcol->getElementAt(i));
				if (vertex->getParameters()[1] > 0) 
				{
					bgenntracks += vertex->getAssociatedParticle()->getParticles().size();
				}
				if (vertex->getParameters()[1] < 0) 
				{
					bbargenntracks += vertex->getAssociatedParticle()->getParticles().size();
				}
			}
			_stats._MCBBarOscillation = mcvtxcol->getParameters().getIntVal("BBarOscillation");
			_stats._MCBOscillation = mcvtxcol->getParameters().getIntVal("BOscillation");
		}
		for (unsigned int i = 0; i < mcbs.size(); i++) 
		{
			MCParticle * mcb = mcbs[i];
			float angleCurrent =  MathOperator::getAngle(mcb->getMomentum(), bjets->at(0)->getMomentum());
			if (angleCurrent < angle) 
			{ 
				angle = angleCurrent;
				charge = mcb->getCharge();
				bjets->at(0)->__SetMCCharge(charge);
			}
		}
		std::cout << "Truth Angle: " << angle << " charge: " << charge << "\n";
		_stats._B1truthAngle = angle;
		bjets->at(1)->__SetMCCharge(-charge);
		//*/
		/*float angle00 = MathOperator::getAngle(mcbs[0]->getMomentum(), bjets->at(0)->getMomentum());
		float angle11 = MathOperator::getAngle(mcbs[1]->getMomentum(), bjets->at(1)->getMomentum());
		//float energy0 = mcbs[0]->getEnergy() - bjets->at(0)->getEnergy() + mcbs[1]->getEnergy() -  bjets->at(1)->getEnergy();
		float angle01 = MathOperator::getAngle(mcbs[0]->getMomentum(), bjets->at(1)->getMomentum());
		float angle10 = MathOperator::getAngle(mcbs[1]->getMomentum(), bjets->at(0)->getMomentum());
		//float energy1 = mcbs[1]->getEnergy() - bjets->at(0)->getEnergy() + mcbs[0]->getEnergy() -  bjets->at(1)->getEnergy();
		std::cout << "Angle btw quarks: " << _stats._bbbarAngle << " a1+a2: " << angle00+angle11 << " b1+b2: " << angle01+angle10 <<"\n"; 
		if (angle00 + angle11 < angle01 + angle10) 
		{
			bjets->at(0)->__SetMCCharge(mcbs[0]->getCharge());
			bjets->at(1)->__SetMCCharge(mcbs[1]->getCharge());
			_stats._B1truthAngle = angle00 + angle11;
		}
		else 
		{
			_stats._B1truthAngle = angle01 + angle10;
			bjets->at(1)->__SetMCCharge(mcbs[0]->getCharge());
			bjets->at(0)->__SetMCCharge(mcbs[1]->getCharge());
			
		}*/
		if ( bjets->at(0)->__GetMCCharge() < 0) 
		{
			bjets->at(0)->__SetMCNtracks(bgenntracks);
			bjets->at(0)->__SetMCOscillation(_stats._MCBOscillation);
			bjets->at(1)->__SetMCNtracks(bbargenntracks);
			bjets->at(1)->__SetMCOscillation(_stats._MCBBarOscillation);
		}
		else 
		{
			bjets->at(0)->__SetMCNtracks(bbargenntracks);
			bjets->at(0)->__SetMCOscillation(_stats._MCBBarOscillation);
			bjets->at(1)->__SetMCNtracks(bgenntracks);
			bjets->at(1)->__SetMCOscillation(_stats._MCBOscillation);
		}
	}
	void TTBarProcessor::MatchB(std::vector< EVENT::MCParticle * > & mcbs, TopQuark * topHadronic, TopQuark * top2, LCCollection * mcvtxcol)
	{
		float charge = 0.0;
		float angle = 4.0;
		int bgenntracks = 0;
		int bbargenntracks = 0;
		_stats._MCBOscillation = 0;
		_stats._MCBBarOscillation = 0;
		if (mcvtxcol && mcvtxcol->getNumberOfElements() > 0) 
		{
			for (int i = 0; i < mcvtxcol->getNumberOfElements(); i++) 
			{
				Vertex * vertex = dynamic_cast< Vertex * >(mcvtxcol->getElementAt(i));
				if (vertex->getParameters()[1] > 0) 
				{
					bgenntracks += vertex->getAssociatedParticle()->getParticles().size();
				}
				if (vertex->getParameters()[1] < 0) 
				{
					bbargenntracks += vertex->getAssociatedParticle()->getParticles().size();
				}
			}
			_stats._MCBBarOscillation = mcvtxcol->getParameters().getIntVal("BBarOscillation");
			_stats._MCBOscillation = mcvtxcol->getParameters().getIntVal("BOscillation");
		}
		for (unsigned int i = 0; i < mcbs.size(); i++) 
		{
			MCParticle * mcb = mcbs[i];
			float angleCurrent =  MathOperator::getAngle(mcb->getMomentum(), topHadronic->GetB()->getMomentum());
			if (angleCurrent < angle) 
			{ 
				angle = angleCurrent;
				charge = mcb->getCharge();
				topHadronic->GetB()->__SetMCCharge(charge);
			}

		}
		std::cout << "Truth Angle: " << angle << " charge: " << charge << "\n";
		if (top2) 
		{
			if (topHadronic->GetB()->__GetMCCharge() < 0) 
			{
				topHadronic->GetB()->__SetMCNtracks(bgenntracks);
				topHadronic->GetB()->__SetMCOscillation(_stats._MCBOscillation);
				top2->GetB()->__SetMCNtracks(bbargenntracks);
				top2->GetB()->__SetMCOscillation(_stats._MCBBarOscillation);
			}
			else 
			{
				topHadronic->GetB()->__SetMCNtracks(bbargenntracks);
				topHadronic->GetB()->__SetMCOscillation(_stats._MCBBarOscillation);
				top2->GetB()->__SetMCNtracks(bgenntracks);
				top2->GetB()->__SetMCOscillation(_stats._MCBOscillation);
			}
			std::cout << "bgenntracks: " << topHadronic->GetB()->__GetMCNtracks()<< " bbargenntracks: " << top2->GetB()->__GetMCNtracks() << "\n";
			top2->GetB()->__SetMCCharge(-charge);
			if ((-charge) * top2->GetW()->getCharge() < 0) 
			{
				_stats._Top2leptonCorrect = 1;
				std::cout << "Charge correct! \n";
			}
		}
	}
	void TTBarProcessor::ClearVariables()
	{
		_stats._qMCcostheta[0] = -2;
		_stats._qMCcostheta[1] = -2;
		_stats._MCMass = -1;
		_stats.Clear();
	}

	void TTBarProcessor::end()
	{   
		_hSumTree->Fill();
		_hfile->Write();
		_hfile->Close();
	}
	string TTBarProcessor::intToStr(int * number)
	{
		std::stringstream ss;
		if (number) 
		{
			ss << *number;
		}
		else 
		{
			ss << "NAN";
		}
		string str = ss.str();
		return str;
	}
	vector<int>  TTBarProcessor::getOpposite(int icandidate, int jcandidate)
	{
		vector<int> result;
		int kcandidate = -1;
		int lcandidate = -1;
		for ( int i = 0; i < 4; i++) 
		{
			if (icandidate != i && jcandidate != i && kcandidate < 0) 
			{
				kcandidate = i;
			}
			if (jcandidate != i && icandidate != i &&  lcandidate < 0 && kcandidate != i)
			{
				lcandidate = i;
			}
		}
		result.push_back(kcandidate);
		result.push_back(lcandidate);
		return result;
	}

	float TTBarProcessor::getChargeBalance(RecoJet * jet)
	{
		int nparticles = jet->getParticles().size();
		float e_c = 0.0;
		float e = 0.0;
		for (unsigned int i = 0; i < nparticles; i++) 
		{
			ReconstructedParticle * particle = jet->getParticles()[i];
			if (std::abs(particle->getCharge()) > 0.0) 
			{
				e_c += particle->getEnergy();
			}
			e += particle->getEnergy();
		}
		return e_c / e;
	}
	void TTBarProcessor::getThrust(vector<float> & thrust, LCCollection * pfos)
	{
		std::cout << "Size: " << thrust.size() << "\n";
		double taxis[3];
		taxis[0] = thrust[0]; taxis[1] = thrust[1]; taxis[2] = thrust[2];
		vector<float> direction = MathOperator::getDirection(thrust);
		float axisAngle = MathOperator::getAngles(direction)[1];
		std::cout << "THRUST COS0: " << std::cos(axisAngle) << "\n";
		int nparticles = pfos->getNumberOfElements();
		int counter = 0;
		vector< ReconstructedParticle * > jetPositive;
		vector< ReconstructedParticle * > jetNegative;
		for (unsigned int i = 0; i < nparticles; i++) 
		{
			ReconstructedParticle * particle = dynamic_cast< ReconstructedParticle * >(pfos->getElementAt(i));
			float angle = MathOperator::getAngleBtw(taxis, particle->getMomentum());
			//std::cout << "Angle: " << angle << "\n";
			if (angle > 3.1416/2) 
			{
				jetPositive.push_back(particle);
				counter++;
			}
			else 
			{
				jetNegative.push_back(particle);
			}
		}
		std::cout << "In positive direction " << counter << " particles, in negative - " << nparticles - counter << "\n";
		_stats._ZZMass1 = getMass(jetPositive);
		_stats._ZZMass2 = getMass(jetNegative);
		std::cout << "Mass_+: " << _stats._ZZMass1 << " Mass_-: " << _stats._ZZMass2 << "\n";
	}
	float TTBarProcessor::getMass(std::vector< EVENT::ReconstructedParticle * > & hemisphere)
	{
		float mass = 0.0;
		float energy  = 0.0;
		float momentum[3] = { 0., 0., 0.};
		for (unsigned int i = 0; i < hemisphere.size(); i++) 
		{
			energy += hemisphere[i]->getEnergy();
			momentum[0] += hemisphere[i]->getMomentum()[0];
			momentum[1] += hemisphere[i]->getMomentum()[1];
			momentum[2] += hemisphere[i]->getMomentum()[2];
		}
		mass = std::sqrt(energy * energy - momentum[0] * momentum[0] - momentum[1]*momentum[1] - momentum[2] * momentum[2] );
		
		return mass;
	}
	void TTBarProcessor::getZZ(LCCollection * fourjetcol)
	{
		if (!fourjetcol || fourjetcol->getNumberOfElements() < 4) 
		{
			return;
		}
		int njets = fourjetcol->getNumberOfElements();

		vector<ReconstructedParticle * > jets;
		vector<float> masses;
		float average = 0.0;
		for (unsigned int i = 0; i < njets; i++) 
		{
			ReconstructedParticle * jetpart = dynamic_cast< ReconstructedParticle * >(fourjetcol->getElementAt(i));
			jets.push_back(jetpart);
			masses.push_back(jetpart->getMass());
			//average += std::sqrt(std::pow(jetpart->getEnergy(),2) -std::pow( MathOperator::getModule(jetpart->getMomentum()),2));
			//std::cout << "m_" << i << ": " << std::sqrt(std::pow(jetpart->getEnergy(),2) -std::pow( MathOperator::getModule(jetpart->getMomentum()),2)) << "\n";
		}
		std::sort(jets.begin(), jets.end(), sortByEnergy);

		int n = 0;
		vector< vector<float> > differences;
		/*vector<ReconstructedParticle * > jetsZ1;
		jetsZ1.push_back(jets[0]);
		jetsZ1.push_back(jets[1]);
		vector<ReconstructedParticle * > jetsZ2;
		jetsZ2.push_back(jets[2]);
		jetsZ2.push_back(jets[3]);
		float mass1 = getMass(jetsZ1);
		float mass2 = getMass(jetsZ2);
		vector<float> pair1;
		pair1.push_back(mass1);
		pair1.push_back(mass2);
		differences.push_back(pair1);*/

		vector<ReconstructedParticle * > jetsZ3;
		jetsZ3.push_back(jets[0]);
		jetsZ3.push_back(jets[2]);
		vector<ReconstructedParticle * > jetsZ4;
		jetsZ4.push_back(jets[1]);
		jetsZ4.push_back(jets[3]);
		float mass3 = getMass(jetsZ3);
		float mass4 = getMass(jetsZ4);
		vector<float> pair2;
		pair2.push_back(mass3);
		pair2.push_back(mass4);
		differences.push_back(pair2);
		
		vector<ReconstructedParticle * > jetsZ5;
		jetsZ5.push_back(jets[0]);
		jetsZ5.push_back(jets[3]);
		vector<ReconstructedParticle * > jetsZ6;
		jetsZ6.push_back(jets[1]);
		jetsZ6.push_back(jets[2]);
		float mass5 = getMass(jetsZ5);
		float mass6 = getMass(jetsZ6);
		vector<float> pair3;
		pair3.push_back(mass5);
		pair3.push_back(mass6);
		differences.push_back(pair3);
		
		std::cout << "Masses: " //<< mass1 << " & " << mass2 << "; "
					<< mass3 << " & " << mass4 << "; "
					<< mass5 << " & " << mass6 << ";\n";
		float mindiff = 1000.;
		int chosen = -1;
		for (unsigned int i = 0; i < differences.size(); i++) 
		{
			float diff = abs(differences[i][0] - differences[i][1]);
			if (diff < mindiff) 
			{
				mindiff = diff;
				chosen = i;
			}
		}
		_stats._ZZMass1 = differences[chosen][0];
		_stats._ZZMass2 = differences[chosen][1];

		for (unsigned int i = 0; i < njets; i++) 
		{
			for (unsigned int j = 1; j < i; j++) 
			{
				
			}
		}
		/*for (unsigned int i = 1; i < njets; i++) 
		{
			for (unsigned int j = 0; j < i; j++) 
			{
				float ei = jets[i]->getEnergy();
				float ej = jets[j]->getEnergy();
				float momentum = 0.;
				for (unsigned int k = 0; k < 3; k++) 
				{
					momentum+= jets[i]->getMomentum()[k]*jets[j]->getMomentum()[k];
				}
				float mass = std::sqrt(ei*ej - momentum);
				average += mass;
				massesZ.push_back(mass);
				n++;
				std::cout << "Z_" << i << "|" << j << ": " << mass << "\n";
			}
		}
		float deviation = 0.0;
		for (unsigned int i = 0; i < masses.size(); i++) 
		{
			deviation += std::pow(average / n - massesZ[i],2);
		}
		_stats._bbbar4JetMass = average / n;//average / n;*/

	}

}
