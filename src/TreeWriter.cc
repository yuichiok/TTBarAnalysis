#include "TreeWriter.hh"
using std::string;
using std::vector;

namespace TTbarAnalysis 
{
	TreeWriter:: TreeWriter() {}
	void TreeWriter::InitializeSummaryTree(TTree * _hSumTree, SummaryData & data)
	{	
		_hSumTree->Branch("nEvents", &data._nEvt, "nEvents/I");
		_hSumTree->Branch("nGenUsed", &data._nGenUsed, "nGenUsed/I");
		_hSumTree->Branch("nAfterLeptonCuts", &data._nAfterLeptonCuts, "nAfterLeptonCuts/I");
		_hSumTree->Branch("nAfterBtagCuts", &data._nAfterBtagCuts, "nAfterBtagCuts/I");
		_hSumTree->Branch("nKaons", &data._nKaons, "nKaons/I");
		_hSumTree->Branch("nChargedB", &data._nChargedB, "nChargedB/I");
		_hSumTree->Branch("nAfterKinematicCuts", &data._nAfterKinematicCuts, "nAfterKinematicCuts/I");
		
	}
	void TreeWriter::InitializeStatsTree(TTree * _hTree, StatsData & data)
	{
		//Generated
		_hTree->Branch("mctag", &data._mctag, "mctag/I");
		_hTree->Branch("MCTopmass", &data._MCTopmass, "MCTopmass/F");
		_hTree->Branch("MCTopmomentum", &data._MCTopmomentum, "MCTopmomentum/F");
		_hTree->Branch("MCTopcostheta", &data._MCTopcostheta, "MCTopcostheta/F");
		_hTree->Branch("MCTopBarmass", &data._MCTopBarmass, "MCTopBarmass/F");
		_hTree->Branch("MCTopBarmomentum", &data._MCTopBarmomentum, "MCTopBarmomentum/F");
		_hTree->Branch("MCTopBarcostheta", &data._MCTopBarcostheta, "MCTopBarcostheta/F");
		_hTree->Branch("MCTopBangle", &data._MCTopBangle, "MCTopBangle/F");
		_hTree->Branch("MCTopcosWb", &data._MCTopcosWb, "MCTopcosWb/F");
		_hTree->Branch("MCNeutrinoEnergy", &data._MCNeutrinoEnergy, "MCNeutrinoEnergy/F");
		_hTree->Branch("MCBOscillation", &data._MCBOscillation, "MCBOscillation/I");
		_hTree->Branch("MCBBarOscillation", &data._MCBBarOscillation, "MCBBarOscillation/I");
		_hTree->Branch("qMCBcostheta", data._qMCBcostheta, "qMCBcostheta[2]/F");
		_hTree->Branch("qMCcostheta", data._qMCcostheta, "qMCcostheta[2]/F");

		//Reconstructed
		_hTree->Branch("totalEnergy", &data._totalEnergy, "totalEnergy/F");
		_hTree->Branch("missedEnergy", &data._missedEnergy, "missedEnergy/F");
		_hTree->Branch("W1mass", &data._W1mass, "W1mass/F");
		_hTree->Branch("W1momentum", &data._W1momentum, "W1momentum/F");
		_hTree->Branch("W1costheta", &data._W1costheta, "W1costheta/F");
		_hTree->Branch("W2momentum", &data._W2momentum, "W2momentum/F");
		_hTree->Branch("Top1mass", &data._Top1mass, "Top1mass/F");
		_hTree->Branch("Top1energy", &data._Top1energy, "Top1energy/F");
		_hTree->Branch("Top1bcharge", &data._Top1bcharge, "Top1bcharge/I");
		_hTree->Branch("Top1bmomentum", &data._Top1bmomentum, "Top1bmomentum/F");
		_hTree->Branch("Top1bdistance", &data._Top1bdistance, "Top1bdistance/F");
		_hTree->Branch("Top1costheta", &data._Top1costheta, "Top1costheta/F");
		_hTree->Branch("Top1bcostheta", &data._Top1bcostheta, "Top1bcostheta/F");
		_hTree->Branch("Top1truthAngle", &data._Top1truthAngle, "Top1truthAngle/F");
		_hTree->Branch("Top1bntracks", &data._Top1bntracks, "Top1bntracks/I");
		_hTree->Branch("Top1bTVCM", &data._Top1bTVCM, "Top1bTVCM/I");
		_hTree->Branch("Top1cosWb", &data._Top1cosWb, "Top1cosWb/F");
		_hTree->Branch("Top2bmomentum", &data._Top2bmomentum, "Top2bmomentum/F");
		_hTree->Branch("Top2bdistance", &data._Top2bdistance, "Top2bdistance/F");
		_hTree->Branch("Top2bcharge", &data._Top2bcharge, "Top2bcharge/I");
		_hTree->Branch("Top2bcostheta", &data._Top2bcostheta, "Top2bcostheta/F");
		_hTree->Branch("Top2bTVCM", &data._Top2bTVCM, "Top2bTVCM/I");
		_hTree->Branch("Top2bntracks", &data._Top2bntracks, "Top2bntracks/I");
		_hTree->Branch("Top2leptonCharge", &data._Top2leptonCharge, "Top2leptonCharge/I");
		_hTree->Branch("Top2leptonCorrect", &data._Top2leptonCorrect, "Top2leptonCorrect/I");
		_hTree->Branch("UsedBTVCM", &data._UsedBTVCM, "UsedBTVCM/I");
		_hTree->Branch("methodUsed", &data._methodUsed, "methodUsed/I");
		_hTree->Branch("methodCorrect", &data._methodCorrect, "methodCorrect/I");
		_hTree->Branch("qCostheta", data._qCostheta, "qCostheta[2]/F");
		_hTree->Branch("chiHad", &data._chiHad, "chiHad/F");
		_hTree->Branch("chiTopMass", &data._chiTopMass, "chiTopMass/F");
		_hTree->Branch("chiTopE", &data._chiTopE, "chiTopE/F");
		_hTree->Branch("chiPbstar", &data._chiPbstar, "chiPbstar/F");
		_hTree->Branch("chiCosWb", &data._chiCosWb, "chiCosWb/F");
		_hTree->Branch("chiGammaT", &data._chiGammaT, "chiGammaT/F");
		
	}
} /* TTBarAnalysis */
