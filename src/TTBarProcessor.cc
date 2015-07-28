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

	}
	
	void TTBarProcessor::init() 
	{ 
		// usually a good idea to
		printParameters() ;

		_nRun = 0 ;
		_nEvt = 0 ;   

		_hfile = new TFile( _hfilename.c_str(), "RECREATE", _hfilename.c_str() ) ;

		_hTree = new TTree( "Stats", "tree" );

		_hTree->Branch("MCTopmass", &_MCTopmass, "MCTopmass/F");
		_hTree->Branch("MCTopmomentum", &_MCTopmomentum, "MCTopmomentum/F");
		_hTree->Branch("MCTopcostheta", &_MCTopcostheta, "MCTopcostheta/F");
		_hTree->Branch("MCTopBarmass", &_MCTopBarmass, "MCTopBarmass/F");
		_hTree->Branch("MCTopBarmomentum", &_MCTopBarmomentum, "MCTopBarmomentum/F");
		_hTree->Branch("MCTopBarcostheta", &_MCTopBarcostheta, "MCTopBarcostheta/F");
		_hTree->Branch("qMCcostheta", _qMCcostheta, "qMCcostheta[2]/F");
		_hTree->Branch("W1mass", &_W1mass, "W1mass/F");
		_hTree->Branch("Top1mass", &_Top1mass, "Top1mass/F");
		_hTree->Branch("Top1charge", &_Top1charge, "Top1charge/I");
		_hTree->Branch("Top1bmomentum", &_Top1bmomentum, "Top1bmomentum/F");
		_hTree->Branch("Top1costheta", &_Top1costheta, "Top1costheta/F");
		_hTree->Branch("Top1bntracks", &_Top1bntracks, "Top1bntracks/I");
		_hTree->Branch("Top2bmomentum", &_Top2bmomentum, "Top2bmomentum/F");
		_hTree->Branch("Top2bntracks", &_Top2bntracks, "Top2bntracks/I");
		_hTree->Branch("qCostheta", _qCostheta, "qCostheta[2]/F");
		_hTree->Branch("chiHad", &_chiHad, "chiHad/F");
		
		//_hTree->Branch("detectedTotal", &_detectedTotal, "detectedTotal/I");
	}
	void TTBarProcessor::processRunHeader( LCRunHeader* run) 
	{ 
		_nRun++ ;
	} 
	void TTBarProcessor::AnalyseGenerator(MCOperator & opera)
	{
		std::vector< EVENT::MCParticle * > mctops = opera.GetTopPairParticles();

		if (mctops.size() < 2) 
		{
			_mctag = 0;
			return;
		}
		MCParticle * top = mctops[0];
		MCParticle * topbar = mctops[1];
		vector<float> direction = MathOperator::getDirection(top->getMomentum());
		vector<float> directionbar = MathOperator::getDirection(topbar->getMomentum());
		_MCTopmass = top->getMass();
		_MCTopBarmass = topbar->getMass();
		_MCTopmomentum = MathOperator::getModule(top->getMomentum());
		_MCTopBarmomentum = MathOperator::getModule(topbar->getMomentum());
		_MCTopcostheta = std::cos( MathOperator::getAngles(direction)[1] );
		_MCTopBarcostheta = std::cos( MathOperator::getAngles(directionbar)[1] );
		std::cout << "Top costheta: " << _MCTopcostheta << "\n";
		std::cout << "TopBar costheta: " << _MCTopBarcostheta << "\n";
		_qMCcostheta[0] = _MCTopcostheta;
		_qMCcostheta[1] = -_MCTopBarcostheta;
		if (_MCTopmass > 170.0 &&_MCTopmass < 181.0 && _MCTopBarmass > 170.0 &&_MCTopBarmass < 181.0) 
		{
			_mctag = 1;
		}
	}
	void TTBarProcessor::processEvent( LCEvent * evt )
	{
		try
		{
			std::cout << "***************Analysis*" <<_nEvt++<<"*****************\n";

		//	LCCollection * isoleptoncol = evt->getCollection(_IsoLeptonColName);
			LCCollection * jetcol = evt->getCollection(_JetsColName);
			LCCollection * jetrelcol = evt->getCollection(_JetsRelColName);
			LCCollection * mccol = evt->getCollection(_MCColName);
			MCOperator opera(mccol);
			AnalyseGenerator(opera);
			vector< RecoJet * > * jets = getJets(jetcol, jetrelcol);
			vector< RecoJet * > * wjets = new vector< RecoJet * >();
			vector< RecoJet * > * bjets = getBTagJets(jets, wjets);
			std::cout << "B jets: \n";
			PrintJets(bjets);
			std::cout << "W jets: \n";
			PrintJets(wjets);
			if (bjets->size() != 2 ) 
			{
				std::cout << "No B jets!!! \n";
				return;
			}
			if (bjets->at(0)->GetBTag() < _highBTagCutparameter && bjets->at(1)->GetBTag() < _highBTagCutparameter) 
			{
				std::cout << "No high B jets!!! \n";
				return;
			}
			float chimin = 100000.0;
			RecoJet * wHadronic = new TopQuark(wjets->at(0), wjets->at(1));
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
				}
			}
			if (chosen) 
			{
				topLeptonic = new TopQuark(bjets->at(0));
			}
			else 
			{
				topLeptonic = new TopQuark(bjets->at(1));
			}
			_W1mass = wHadronic->getMass();
			_Top1mass = topHadronic->getMass();
			_chiHad = chimin;
			ComputeCharge(topHadronic, topLeptonic);
			_hTree->Fill();
			ClearVariables();
		}
		catch(DataNotAvailableException &e)
		{
			std::cout << e.what() <<"\n";
		}
	}
	void TTBarProcessor:: check( LCEvent * evt ) 
	{
	}
	void TTBarProcessor::ComputeCharge(TopQuark * top1, TopQuark * top2)
	{
		//Top hadronic
		_Top1charge = top1->GetHadronCharge();
		_Top1bmomentum = top1->GetHadronMomentum();
		vector<float> direction = MathOperator::getDirection(top1->getMomentum());
		_Top1costheta =  std::cos( MathOperator::getAngles(direction)[1] );
		_Top1bntracks = top1->GetNumberOfVertexParticles();
		_Top1btag = top1->GetBTag();
		std::cout << "Top charge: " << _Top1charge
			  << " top pB: " << _Top1bmomentum
			  << " top bntracks: " << _Top1bntracks
			  << " top btag: " << _Top1btag
			  << " top costheta: " << _Top1costheta
			  << "\n";
		//Top Leptonic
		_Top2charge = top2->GetHadronCharge();
		_Top2bmomentum = top2->GetHadronMomentum();
		vector<float> direction2 = MathOperator::getDirection(top2->getMomentum());
		_Top2costheta =  std::cos( MathOperator::getAngles(direction2)[1] );
		_Top2bntracks = top2->GetNumberOfVertexParticles();
		_Top2btag = top2->GetBTag();
		std::cout << "Top charge: " << _Top2charge
			  << " top pB: " << _Top2bmomentum
			  << " top bntracks: " << _Top2bntracks
			  << " top btag: " << _Top2btag
			  << " top costheta: " << _Top2costheta
			  << "\n";
		
		bool trustTop1 = false;
		bool trustTop2 = false;

		if (_Top1bmomentum > 20.0 && abs( _Top1charge) > 0 && _Top1bntracks > 2 && _Top1bntracks < 7 && _Top1btag > 0.9 ) 
		{
			std::cout << "Use top 1!\n";
			trustTop1 = true;
		}
		if (_Top2bmomentum > 20.0 && abs( _Top2charge) > 0 && _Top2bntracks > 2 && _Top2bntracks < 7 && _Top2btag > 0.9 ) 
		{
			std::cout << "Use top 2!\n";
			trustTop2 = true;
		}
		if (trustTop1 && !trustTop2) 
		{
			_qCostheta[0] = (_Top1charge < 0)? _Top1costheta : -_Top1costheta;
		}
		if (trustTop2 && !trustTop1) 
		{
			_qCostheta[0] = (_Top2charge > 0)? _Top1costheta : -_Top1costheta;
		}
		if (trustTop1 && trustTop2) 
		{
			if (_Top1bntracks < _Top2bntracks) 
			{
				//_qCostheta[0] = (_Top1charge < 0)? _Top1costheta : -_Top1costheta;
			}
			if (_Top2bntracks > _Top1bntracks) 
			{
				//_qCostheta[0] = (_Top2charge > 0)? _Top1costheta : -_Top1costheta;
			}
			if (_Top2bntracks == _Top1bntracks) 
			{
				if (_Top1bmomentum > _Top2bmomentum) 
				{
					//_qCostheta[0] = (_Top1charge < 0)? _Top1costheta : -_Top1costheta;
				}
				else 
				{
					//_qCostheta[0] = (_Top2charge > 0)? _Top1costheta : -_Top1costheta;
				}
			}
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
		float chi2 = std::pow( mT - _TopMassparameter, 2) / std::pow( _TopMassSigmaparameter, 2) + 
			     std::pow(ET - _EBeamparameter , 2) / std::pow( _EBeamSigmaparameter, 2) +
			     std::pow( bpstar - _PStarparameter, 2) / std::pow( _PStarSigmaparameter, 2) +
			     std::pow( cosbW - _CosbWparameter, 2) / std::pow( _CosbWSigmaparameter, 2);
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
		for (unsigned int i = 0; i < alljets->size(); i++) 
		{
			RecoJet * jet = alljets->at(i);
			if (jet->GetBTag() > _lowBTagCutparameter) 
			{
				result->push_back(jet);
			}
			else 
			{
				wjets->push_back(jet);
			}
		}
		if (result->size() < 3) 
		{
			return result;
		}
		RecoJet * bjet1 = NULL;
		RecoJet * bjet2 = NULL;
		if (result->size() > 2) 
		{
			float maxbtag = 0.0;
			for (unsigned int i = 0; i < result->size(); i++) 
			{
				RecoJet * jet = result->at(i);
				if (jet->GetBTag() > maxbtag) 
				{
					bjet1 = jet;
					maxbtag = jet->GetBTag();
				}
			}
			maxbtag = 0.0;
			for (unsigned int i = 0; i < result->size(); i++) 
			{
				RecoJet * jet = result->at(i);
				if (jet != bjet1 && jet->GetBTag() > maxbtag) 
				{
					bjet2 = jet;
					maxbtag = jet->GetBTag();
				}
			}  
			for (unsigned int i = 0; i < result->size(); i++) 
			{
				RecoJet * jet = result->at(i);
				if (jet != bjet1 && jet != bjet2) 
				{
					wjets->push_back(jet);
				}
			}  
		}
		result->clear();
		result->push_back(bjet1);
		result->push_back(bjet2);
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
			  << " ctag: " << jet->GetCTag()
			  << " ntracks: " << jet->GetNumberOfVertexParticles()
			  << "\n";
	}
	void TTBarProcessor::PrintJets(std::vector< RecoJet * > *jets)
	{
		for (unsigned int i = 0; i < jets->size(); i++) 
		{
			PrintJet(jets->at(i));
		}
	}

	void TTBarProcessor::ClearVariables()
	{
		_MCTopmass = -1.0;
		_MCTopBarmass = -1.0;
		_MCTopmomentum = -1.0;
		_MCTopBarmomentum = -1.0;
		_MCTopcostheta = -2.0;
		_MCTopBarcostheta = -2.0;
		_qMCcostheta[0] = -2.0;
		_qMCcostheta[1] = -2.0;

		_Top1mass = -1.0;
		_W1mass = -1.0;
		_chiHad = -1.0;
		_qCostheta[0] = -2.0;
		_qCostheta[1] = -2.0;
	}

	void TTBarProcessor::end()
	{   
		_hfile->Write();
		_hfile->Close();
	}
}
	
