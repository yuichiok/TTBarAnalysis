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
namespace TTbarAnalysis
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
		_type = 2;
		registerProcessorParameter( "AnalysisType",
				      "Analysis Type",
				      _type,
				      _type);
		
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
				break;
		}
		//AnalyseTTBarSemiLeptonic(evt);

		
		//_hTree->Branch("detectedTotal", &_stats._detectedTotal, "detectedTotal/I");
	}
	void TTBarProcessor::processRunHeader( LCRunHeader* run) 
	{ 
		_nRun++ ;
	} 
	vector< MCParticle * > TTBarProcessor::AnalyseGeneratorBBBar(MCOperator & opera)
	{
		vector <MCParticle *> genfinalstate = opera.GetPairParticles(5);
		for (unsigned int i = 0; i < genfinalstate.size(); i++) 
		{
			PrintParticle(genfinalstate[i]);
			vector<float> direction = MathOperator::getDirection(genfinalstate[i]->getMomentum());
			int charge = genfinalstate[i]->getCharge() / std::abs(genfinalstate[i]->getCharge());
			_stats._qMCcostheta[i] = - std::cos( MathOperator::getAngles(direction)[1]) * charge;
			if (_stats._qMCcostheta[i] < -0.6 && _stats._qMCcostheta[i] > -2) 
			{
				std::cout << "Bkg process: " << _stats._qMCcostheta[i] << "\n";
			}
		}
		if (genfinalstate.size() == 2) 
		{
			MCParticle * Z = opera.CombineParticles(genfinalstate[0],genfinalstate[1]);
			_stats._MCMass = Z->getMass();
		}
		_hGenTree->Fill();
		return genfinalstate;
		
	}
	vector< MCParticle * > TTBarProcessor::AnalyseGenerator(MCOperator & opera)
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
				break;
		}
		//AnalyseTTBarSemiLeptonic(evt);
	}

	void TTBarProcessor::AnalyseTTBarSemiLeptonic( LCEvent * evt )
	{
		try
		{
			std::cout << "***************Analysis*" <<_summary._nEvt++<<"*****************\n";

			LCCollection * isoleptoncol = evt->getCollection(_IsoLeptonColName);
			LCCollection * jetcol = evt->getCollection(_JetsColName);
			LCCollection * jetrelcol = evt->getCollection(_JetsRelColName);
			LCCollection * mccol = evt->getCollection(_MCColName);
			LCCollection * mcvtxcol = evt->getCollection(_MCVtxColName);
			MCOperator opera(mccol);
			VertexChargeOperator vtxOperator(evt->getCollection(_colName),evt->getCollection(_colRelName));
			vector < MCParticle * > mctops = AnalyseGenerator(opera);
			vector< RecoJet * > * jets = getJets(jetcol, jetrelcol);
			vector< RecoJet * > * wjets = new vector< RecoJet * >();
			vector< RecoJet * > * bjets = getBTagJets(jets, wjets);
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
			Match(mctops, topHadronic);
			vector< EVENT::MCParticle * > mcbquarks = opera.GetBquarkPair();
			MatchB(mcbquarks, topHadronic, topLeptonic, mcvtxcol);
			ComputeCharge(topHadronic, topLeptonic);
			ComputeChargeTVCM(topHadronic, topLeptonic, vtxOperator);

			DecideOnAsymmetry(topHadronic, topLeptonic);
			if (_stats._methodUsed && wHadronic->getMass() < 110 && topHadronic->getMass() < 200) 
			{
				_summary._nAfterMassCuts++;
			}
			//test(topHadronic, topLeptonic);
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
		//_stats._qCostheta[0] =  costheta;
		//_stats._Top1bTVCM = top->GetResultTVCM();
		//_stats._Top2bTVCM = top2->GetResultTVCM();
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
		/*vector<float> directionTop = MathOperator::getDirection(top->getMomentum());
		_stats._Top1costheta =  std::cos( MathOperator::getAngles(directionTop)[1] );
		_stats._Top1bdistance = top->GetMaxHadronDistance();
		_stats._Top2bdistance = top2->GetMaxHadronDistance();
		_stats._Top1bntracks = top->GetNumberOfVertexParticles();
		_stats._Top2bntracks = top2->GetNumberOfVertexParticles();
		*/
		_stats._qCostheta[0] = -2;
		std::cout << "\t Costheta: " << _stats._Top1costheta << "\n";
		if (top->GetComputedCharge().ByTVCM) 
		{
			int top1charge = *(top->GetComputedCharge().ByTVCM );
			 _stats._qCostheta[0] = (top1charge < 0)? _stats._Top1costheta : -_stats._Top1costheta;
			_stats._methodUsed = 1;
			_stats._methodCorrect = top->__GetMCCharge() * top1charge < 0;
			_summary._nAfterKinematicCuts++;
		}
		if (top2->GetComputedCharge().ByTVCM) 
		{
			int top2charge = *(top2->GetComputedCharge().ByTVCM );
			 _stats._qCostheta[0] = (top2charge < 0)? -_stats._Top1costheta : _stats._Top1costheta;
			_stats._methodUsed = 2;
			_stats._methodCorrect = top2->__GetMCCharge() * top2charge < 0;
			_summary._nAfterKinematicCuts++;
		}
		if (top2->GetComputedCharge().ByTVCM && top->GetComputedCharge().ByTVCM) 
		{
			int top1charge = *(top->GetComputedCharge().ByTVCM );
			int top2charge = *(top2->GetComputedCharge().ByTVCM );
			if (top1charge * top2charge < 0) 
			{
				_stats._qCostheta[0] = (top1charge < 0)? _stats._Top1costheta : -_stats._Top1costheta;
			}
			_stats._methodUsed = 3;
			_stats._methodCorrect = top->__GetMCCharge() * top1charge < 0;
			_summary._nAfterKinematicCuts++;
		}
		std::cout << "\t Costheta final: " << _stats._qCostheta[0]  << "\n";
	}
	void TTBarProcessor:: check( LCEvent * evt ) 
	{
	}
	void TTBarProcessor::test(TopQuark * top, TopQuark * top2)
	{	
		_totalEnergy = top->getEnergy() + top2->getEnergy();
		_missedEnergy = 2*_EBeamparameter - _totalEnergy;
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
		
		bool trustTop1 = false;
		bool trustTop2 = false;
		//float pcut = 15.0;
		//float pcut2 = 15.0;
		//float btagcut = 0.8;
		//float btagcut2 = 0.8;
		int nvtxcut = 0;
		int nvtxcut2 = 0;
		if (_stats._Top1bntracks %2 == 1 && _stats._Top1bntracks != 3) 
		{
			nvtxcut = 1;
			//pcut = 20.;
			//btagcut = 0.99;
		}
		if (_stats._Top2bntracks %2 == 1 && _stats._Top2bntracks != 3) 
		{
			nvtxcut2 = 1;
			//pcut2 = 20.;
			//btagcut2 = 0.99;
		}
		if (//_stats._Top1bmomentum > pcut && 
		    abs( _stats._Top1bcharge) > 0 && 
		    _stats._Top1bntracks > 2 //&& 

		    //top1->GetB()->__GetMCNtracks() > 0 &&
		    //Top1nvtx > nvtxcut &&
		    //_stats._Top1bntracks < 9 &&
		    //_stats._Top1bntracks != 7 &&
		    ///_stats._Top1bcostheta < 0.9 &&
		    //_stats._Top1bdistance > 0.5 &&
		    //(_stats._Top1bntracks %2 == 0 || _stats._Top1bntracks == 3) && 
		    //_stats._Top1btag > btagcut 
		    )  //  _stats._Top1bntracks < 7 
		{
			std::cout << "Use top 1!\n";
			trustTop1 = true;
			_stats._Top1Vtx = 1;
			TopCharge & topCharge = top1->GetComputedCharge();
			topCharge.ByTrackCount = new int(top1->GetHadronCharge());
			
		}
		if (//_stats._Top2bmomentum > pcut2 && 
		    abs( _stats._Top2bcharge) > 0 && 
		    _stats._Top2bntracks > 2 //&& 
		    //top2->GetB()->__GetMCNtracks() > 0 &&
		    //Top2nvtx > nvtxcut2 &&
		    //_stats._Top2bntracks < 9 &&
		    //_stats._Top2bntracks != 7 &&
		    //_stats._Top2bcostheta < 0.9 &&
		    //_stats._Top2bdistance > 0.5 &&
		    //(_stats._Top2bntracks %2 == 0 || _stats._Top2bntracks == 3) && 
		    //_stats._Top2btag > btagcut2 
		    ) 
		{
			std::cout << "Use top 2!\n";
			trustTop2 = true;
			_stats._Top2Vtx = 1;
			TopCharge & topCharge = top2->GetComputedCharge();
			topCharge.ByTrackCount = new int(top2->GetHadronCharge());
		}
		std::cout << "\tTop1 weighted charge: " << top1->GetHadronCharge(true) << "\n";
		std::cout << "\tTop2 weighted charge: " << top2->GetHadronCharge(true) << "\n";
		if ( _stats._Top2bcharge != 0 ||  _stats._Top1bcharge != 0) 
		{
			_summary._nChargedB++;
		}
		if (trustTop1 || trustTop2) 
		{
			//_summary._nAfterKinematicCuts++;
		}
		if (top2->GetW()) 
		{
			top2->GetComputedCharge().ByLepton = new int(top2->GetW()->getCharge());
		}
		bool useLepton = (_ePolarization > 0)? false : true; //Lepton by default
		/*if (trustTop2 && _stats._Top2bcharge * top2->GetW()->getCharge() < 0.0) 
		{
			_stats._qCostheta[0] = (_stats._Top2bcharge > 0)? _stats._Top1costheta : -_stats._Top1costheta;
		}
		if (trustTop1 && _stats._Top1bcharge * top2->GetW()->getCharge() > 0.0) 
		{
			_stats._qCostheta[0] = (_stats._Top1bcharge < 0)? _stats._Top1costheta : -_stats._Top1costheta;
		}
		
		//*/
		if (trustTop1 && !trustTop2 && !useLepton) 
		{
			_stats._qCostheta[0] = (_stats._Top1bcharge < 0)? _stats._Top1costheta : -_stats._Top1costheta;
		}
		if (trustTop2 && !trustTop1 && !useLepton) 
		{
			_stats._qCostheta[0] = (_stats._Top2bcharge > 0)? _stats._Top1costheta : -_stats._Top1costheta;
		}
		if (trustTop1 && trustTop2 && !useLepton) 
		{
			if (_stats._Top1bntracks < _stats._Top2bntracks || (_stats._Top1bntracks % 2 == 0 && _stats._Top2bntracks % 2 == 1)) 
			{
				_stats._qCostheta[0] = (_stats._Top1bcharge < 0)? _stats._Top1costheta : -_stats._Top1costheta;
			}
			if (_stats._Top1bntracks > _stats._Top2bntracks || ( _stats._Top1bntracks % 2 == 1 && _stats._Top2bntracks % 2 == 0)) 
			{
				_stats._qCostheta[0] = (_stats._Top2bcharge > 0)? _stats._Top1costheta : -_stats._Top1costheta;
			}
			if (_stats._Top2bntracks == _stats._Top1bntracks) 
			{
				if (Top1nvtx == 2 && Top2nvtx == 1) 
				{
					_stats._qCostheta[0] = (_stats._Top1bcharge < 0)? _stats._Top1costheta : -_stats._Top1costheta;
				}
				if (Top2nvtx == 2 && Top1nvtx == 1) 
				{
					_stats._qCostheta[0] = (_stats._Top2bcharge > 0)? _stats._Top1costheta : -_stats._Top1costheta;
				}
				if (_stats._Top1bmomentum > _stats._Top2bmomentum && Top1nvtx == Top2nvtx) 
				{
					_stats._qCostheta[0] = (_stats._Top1bcharge < 0)? _stats._Top1costheta : -_stats._Top1costheta;
				}
				if (_stats._Top2bmomentum > _stats._Top1bmomentum && Top2nvtx == Top1nvtx) 
				{
					_stats._qCostheta[0] = (_stats._Top2bcharge > 0)? _stats._Top1costheta : -_stats._Top1costheta;
				}

			}
		}//*/
		if ((trustTop1 ||  trustTop2) && _stats._qCostheta[0] < -1.0) 
		{
			std::cout << "ERROR!\n";
		}
		if (!trustTop1 && !trustTop2 && !useLepton) 
		{
			if (_stats._Top2bcharge * top2->GetW()->getCharge() < 0.0) 
			{
				//_stats._qCostheta[0] = (_stats._Top2bcharge > 0)? _stats._Top1costheta : -_stats._Top1costheta;
			}
		}
		if (useLepton) // || _stats._chiHad < 20) 
		{
			//_stats._qCostheta[0] = (top2->GetW()->getCharge() > 0)? -_stats._Top1costheta : _stats._Top1costheta;
		}
	}
	void TTBarProcessor::DecideOnAsymmetry(TopQuark * top1, TopQuark * top2)
	{
		//Print
		
		std::cout << "\t\tTracks\tTVCM\tLepton\n";
		//std::cout << "Top1:\t" << (top1->GetComputedCharge().ByTrackCount)?top1->GetComputedCharge().ByTrackCount:"NAN\t"
		//			<< (top1->GetComputedCharge().ByTVCM)?top1->GetComputedCharge().ByTVCM:"NAN\t"
		//			<< (top1->GetComputedCharge().ByLepton)?top1->GetComputedCharge().ByLepton:"NAN\t"
		std::cout << "Top1:\t" << intToStr(top1->GetComputedCharge().ByTrackCount) <<"\t"
					<< intToStr(top1->GetComputedCharge().ByTVCM) <<"\t"
					<< intToStr(top1->GetComputedCharge().ByLepton) <<"\n";
		std::cout << "Top2:\t" << intToStr(top2->GetComputedCharge().ByTrackCount) <<"\t"
					<< intToStr(top2->GetComputedCharge().ByTVCM) <<"\t"
					<< intToStr(top2->GetComputedCharge().ByLepton) <<"\n";
		_stats._qCostheta[0] = -2.;

		vector<float> direction = MathOperator::getDirection(top1->getMomentum());
		float costheta =  std::cos( MathOperator::getAngles(direction)[1] );
		float m = top1->getMass();
		float e = top1->getEnergy();
		float btagcut = 0.8;
		float pcut = 20.;
		//if (e/m > _GammaTparameter -0.1 && top2->GetComputedCharge().ByLepton) 
		//Track charge * Track charge
		if (top2->GetComputedCharge().ByTrackCount && top1->GetComputedCharge().ByTrackCount) 
		{
			int top1charge = *(top1->GetComputedCharge().ByTrackCount );
			int top2charge = *(top2->GetComputedCharge().ByTrackCount );
			if (top1charge * top2charge < 0 && _stats._Top1btag > btagcut && _stats._Top1bmomentum > pcut) 
			{
				_stats._qCostheta[0] = (top1charge < 0)? costheta: - costheta;
				std::cout << "Two vertices are used!\n";
				_stats._methodUsed = 1;
				_stats._methodCorrect = top1->__GetMCCharge() * top1charge < 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				_summary._nAfterKinematicCuts++;
				return;
			}
		}
		
		//Kaon * Kaon
		if (top2->GetComputedCharge().ByTVCM && top1->GetComputedCharge().ByTVCM) 
		{
			int top1charge = *(top1->GetComputedCharge().ByTVCM );
			int top2charge = *(top2->GetComputedCharge().ByTVCM );
			if (top1charge * top2charge < 0) 
			{
				_stats._qCostheta[0] = (top1charge < 0)? costheta: - costheta;
				std::cout << "Two kaons are used!\n";
				_stats._methodUsed = 2;
				_stats._methodCorrect = top1->__GetMCCharge() * top1charge < 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				_summary._nAfterKinematicCuts++;
				return;
			}
		}//*/
		//Track charge + Kaon
		if (top1->GetComputedCharge().ByTrackCount && top1->GetComputedCharge().ByTVCM) 
		{
			int top1charge = *(top1->GetComputedCharge().ByTrackCount );
			int top1kaon = *(top1->GetComputedCharge().ByTVCM );
			if (top1charge * top1kaon > 0 && _stats._Top1btag > btagcut && _stats._Top1bmomentum > pcut) 
			{
				_stats._qCostheta[0] = (top1charge < 0)? costheta: -costheta; 
				std::cout << "Vertex + kaon for top1 is used!\n";
				_stats._methodUsed = 3;
				_stats._methodCorrect = top1->__GetMCCharge() * top1charge < 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				_summary._nAfterKinematicCuts++;
				return;
			}
		}
		if (top2->GetComputedCharge().ByTrackCount && top2->GetComputedCharge().ByTVCM) 
		{
			int top2charge = *(top2->GetComputedCharge().ByTrackCount );
			int top2kaon = *(top2->GetComputedCharge().ByTVCM );
			if (top2charge * top2kaon > 0 && _stats._Top2btag > btagcut && _stats._Top2bmomentum > pcut) 
			{
				_stats._qCostheta[0] = (top2charge > 0)? costheta: - costheta; 
				std::cout << "Vertex + kaon for top2 is used!\n";
				_stats._methodUsed = 4;
				_stats._methodCorrect = top1->__GetMCCharge() * top2charge > 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				_summary._nAfterKinematicCuts++;
				return;
			}
		}
		if (top1->GetComputedCharge().ByTrackCount && top2->GetComputedCharge().ByTVCM) 
		{
			int top1charge = *(top1->GetComputedCharge().ByTrackCount );
			int top2kaon = *(top2->GetComputedCharge().ByTVCM );
			if (top1charge * top2kaon < 0 &&  _stats._Top1btag > btagcut && _stats._Top1bmomentum > pcut) 
			{
				_stats._qCostheta[0] = (top1charge < 0)? costheta: - costheta; 
				std::cout << "Vertex1 + kaon2 is used!\n";
				_stats._methodUsed = 5;
				_stats._methodCorrect = top1->__GetMCCharge() *  top1charge < 0;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				_summary._nAfterKinematicCuts++;
				return;
			}
		}
		if (top2->GetComputedCharge().ByTrackCount && top1->GetComputedCharge().ByTVCM) 
		{
			int top2charge = *(top2->GetComputedCharge().ByTrackCount );
			int top1kaon = *(top1->GetComputedCharge().ByTVCM );
			if (top2charge * top1kaon < 0 &&  _stats._Top2btag > btagcut && _stats._Top2bmomentum > pcut) 
			{
				_stats._qCostheta[0] = (top2charge > 0)? costheta: - costheta; 
				std::cout << "Vertex2 + kaon1 is used!\n";
				_stats._methodUsed = 6;
				if (!_stats._methodCorrect) 
				{
					 std::cout << "Not Correct!\n";
				}
				_stats._methodCorrect = top1->__GetMCCharge() * top1kaon< 0;
				_summary._nAfterKinematicCuts++;
				return;
			}
		}//*/
		//LEPTON
		/*if (top2->GetComputedCharge().ByLepton && top2->GetComputedCharge().ByTrackCount)
		{
			int top2lepton = *(top2->GetComputedCharge().ByLepton );
			int top2charge = *(top2->GetComputedCharge().ByTrackCount );
			if (top2charge * top2lepton < 0)
			{
				_stats._qCostheta[0] = (top2charge > 0)? costheta: - costheta; 
				std::cout << "Vertex + kaon for top2 is used!\n";
				_stats._methodUsed = 7;
				_stats._methodCorrect = top1->__GetMCCharge() * top2charge < 0;
				_summary._nAfterKinematicCuts++;
				return;
			}
		}
		if (top2->GetComputedCharge().ByLepton && top1->GetComputedCharge().ByTrackCount)
		{
			int top2lepton = *(top2->GetComputedCharge().ByLepton );
			int top1charge = *(top1->GetComputedCharge().ByTrackCount );
			if (top1charge * top2lepton > 0)
			{
				_stats._qCostheta[0] = (top1charge < 0)? costheta: - costheta; 
				std::cout << "Vertex + kaon for top2 is used!\n";
				_stats._methodUsed = 7;
				_stats._methodCorrect = top1->__GetMCCharge() * top1charge < 0;
				_summary._nAfterKinematicCuts++;
				return;
			}
		}//
		// Lepton only
		_stats._gammaT = e/m;
		if (top2->GetComputedCharge().ByLepton && _stats._chiHad < 15) 
		{
			int top2lepton = *(top2->GetComputedCharge().ByLepton );
			_stats._qCostheta[0] = (top2lepton < 0)? costheta: - costheta;
			_stats._methodUsed = 8;
			_summary._nAfterKinematicCuts++;
			return;
		}//*/
		_stats._methodUsed = 0;
		
	}
	void TTBarProcessor::AnalyseBBBar(LCEvent * evt)
	{
		try
		{
			std::cout << "***************Analysis*" <<_summary._nEvt++<<"*****************\n";
			LCCollection * jetcol = evt->getCollection(_JetsColName);
			LCCollection * jetrelcol = evt->getCollection(_JetsRelColName);
			LCCollection * mccol = evt->getCollection(_MCColName);
			LCCollection * mcvtxcol = evt->getCollection(_MCVtxColName);
			MCOperator opera(mccol);
			VertexChargeOperator vtxOperator(evt->getCollection(_colName),evt->getCollection(_colRelName));
			vector < MCParticle * > mcbs = AnalyseGeneratorBBBar(opera);
			vector< RecoJet * > * jets = getJets(jetcol, jetrelcol);
			std::sort(jets->begin(), jets->end(), sortByBtag);
			if (mcbs.size() != 2 || jets->at(0)->GetBTag() < 0.8 ) 
			{
				return;
			}
			std::cout << "--EVENT ACCEPTED--\n";
			_stats._B1momentum = jets->at(0)->GetHadronMomentum(); //MathOperator::getModule(jets->at(0)->getMomentum());
			_stats._B2momentum = jets->at(1)->GetHadronMomentum();//MathOperator::getModule(jets->at(1)->getMomentum());
			_stats._B1btag = jets->at(0)->GetBTag();
			_stats._B2btag = jets->at(1)->GetBTag();
			_stats._B1charge =  jets->at(0)->GetHadronCharge();
			_stats._B2charge =  jets->at(1)->GetHadronCharge();
			vector<float> bdirection1 = MathOperator::getDirection(jets->at(0)->getMomentum());
			vector<float> bdirection2 = MathOperator::getDirection(jets->at(1)->getMomentum());
			_stats._B1costheta =( std::cos( MathOperator::getAngles(bdirection1)[1] ));
			_stats._B2costheta =( std::cos( MathOperator::getAngles(bdirection2)[1] ));
			MatchB(jets, mcbs, mcvtxcol);
			ComputeCharge(jets);
			_hTree->Fill();
			ClearVariables();


		}
		catch(DataNotAvailableException &e)
		{
			std::cout << e.what() <<"\n";
		}
		
	}
	void TTBarProcessor::ComputeCharge(std::vector< RecoJet * > *jets)
	{
		//std::sort(jets->begin(), jets->end(), sortByCostheta);
		PrintJets(jets);
		if (jets->at(0)->__GetMCNtracks() != jets->at(0)->GetNumberOfVertexParticles()) 
		{
			return;	
		}
		/*if (jets->at(0)->__GetMCCharge() < 0) 
		{
			_stats._qCostheta[0] = _stats._B1costheta;
			_stats._qCostheta[1] =-_stats._B2costheta;
		}
		else 
		{
			_stats._qCostheta[0] = -_stats._B1costheta;
			_stats._qCostheta[1] = _stats._B2costheta;

		}
		return;*/
		if (_stats._B1charge * _stats._B2charge < 0 && _stats._B1momentum > 0 && _stats._B2momentum > 0) 
		{
			int b1charge = _stats._B1charge  / std::abs(_stats._B1charge);
			int b2charge = _stats._B2charge  / std::abs(_stats._B2charge);
			_stats._qCostheta[0] = -_stats._B1costheta * b1charge;
			_stats._qCostheta[1] = -_stats._B2costheta * b2charge;
		}
		
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
		_stats._chiTopMass = std::pow(mT - _TopMassparameter, 2) / std::pow( _TopMassSigmaparameter, 2) ;
		_stats._chiTopE = std::pow(ET - _EBeamparameter , 2) / std::pow( _EBeamSigmaparameter, 2); 
		_stats._chiPbstar = std::pow( bpstar - _PStarparameter, 2) / std::pow( _PStarSigmaparameter, 2);
		_stats._chiGammaT = std::pow( gamma - _GammaTparameter, 2) / std::pow( _GammaTSigmaparameter, 2); 
		_stats._chiCosWb = std::pow( cosbW - _CosbWparameter, 2) / std::pow( _CosbWSigmaparameter, 2);
		float chi2 = _stats._chiTopMass + _stats._chiTopE + _stats._chiPbstar  + _stats._chiCosWb;
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
		if (alljets->size() != 4) 
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
		for (int i = 0; i < alljets->size(); i++) 
		{
			RecoJet * jet = alljets->at(i);
			std::cout << "\tB tag: " << jet->GetBTag()  << "\n";

		}
		result->push_back(alljets->at(0));
		result->push_back(alljets->at(1));
		wjets->push_back(alljets->at(2));
		wjets->push_back(alljets->at(3));
		return result;
	}
	vector< RecoJet * > * TTBarProcessor::getJets(LCCollection * jetcol, LCCollection *jetrelcol)
	{
		int jetnumber = jetcol->getNumberOfElements();
		vector< RecoJet * > * result = new vector< RecoJet * >();
		LCRelationNavigator navigator(jetrelcol);
		PIDHandler pidh(jetcol);
		int alid = 0;
		try
		{
			alid = pidh.getAlgorithmID("vtxrec");
		}
		catch(UTIL::UnknownAlgorithm &e)
		{
			std::cout << "No algorithm vtxrec!\n";
			alid = pidh.getAlgorithmID("lcfiplus");
		}
		for (int j = 0; j < jetnumber; j++) 
		{
			ReconstructedParticle * jetpart = dynamic_cast< ReconstructedParticle * >(jetcol->getElementAt(j));
			vector< Vertex * > * vertices = convert(navigator.getRelatedToObjects(jetpart));
			const vector< ReconstructedParticle * > components = jetpart->getParticles();
			int nvtx = vertices->size();
			const ParticleID& pid = pidh.getParticleID(jetpart,alid);
			vector<float> params = pid.getParameters();
			float btag = params[pidh.getParameterIndex(alid,"BTag")];
			float ctag = params[pidh.getParameterIndex(alid,"CTag")];
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
	void TTBarProcessor::Match(std::vector< EVENT::MCParticle * > & mctops, TopQuark * topHadronic,  TopQuark * top2)
	{
		_stats._Top1truthAngle = 4.0;
		float charge = 0.0;

		for (int i = 0; i < mctops.size(); i++) 
		{
			MCParticle * mctop = mctops[i];
			float angle =  MathOperator::getAngle(mctop->getMomentum(), topHadronic->getMomentum());
			if (angle < _stats._Top1truthAngle) 
			{
				_stats._Top1truthAngle = angle;
				charge = mctop->getCharge();
				topHadronic->__SetMCCharge(charge);
			}
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
		for (int i = 0; i < mcbs.size(); i++) 
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
		for (int i = 0; i < mcbs.size(); i++) 
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
}
