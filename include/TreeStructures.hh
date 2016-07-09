#ifndef _SummaryData_hh
#define _SummaryData_hh
namespace TTbarAnalysis 
{
	struct SummaryData 
	{
		public: 
		int _nEvt ;
		int _nGenUsed;
		int _nAfterBtagCuts;
		int _nAfterKinematicCuts;
		int _nAfterLeptonCuts;
		int _nAfterMassCuts;
		int _nChargedB;
		int _nKaons;
		void Clear()
		{
			_nEvt  = 0;
			_nGenUsed = 0;
			_nAfterBtagCuts = 0;
			_nAfterKinematicCuts = 0;
			_nAfterLeptonCuts = 0;
			_nAfterMassCuts = 0;
			_nChargedB = 0;
			_nKaons = 0;

		}
	};
	struct StatsData 
	{
		float _B1momentum;
		float _B2momentum;
		float _B2btag;
		float _B1btag;
		float _B1costheta;
		float _B2costheta;
		float _B1truthAngle;
		int _B1charge;
		int _B2charge;


		int _mctag;
		int _recotag;
		float _W1mass;
		float _W1momentum;
		float _W1costheta;
		float _W2mass;
		float _W2momentum;
		float _W2costheta;
		int _Top1bcharge;
		float _Top1btag;
		float _Top1bmomentum;
		int _Top1bntracks;
		float _Top1mass;
		float _Top1momentum;
		float _Top1energy;
		float _Top1bdistance;
		float _Top1costheta;
		float _Top1bcostheta;
		float _Top1cosWb;
		float _Top1truthAngle;
		int _Top1Vtx;
		int _Top1Kaon;
		int _Top1KaonNumber;
		int _Top1KaonCharges[10];
		float _Top1KaonMomentum[10];
		int _Top1bTVCM;
		int _Top2bntracks;
		int _Top2bTVCM;
		int _UsedBTVCM;
		float _Top2bmomentum;
		float _Top2bdistance;
		float _Top2btag;
		float _Top2mass;
		int _Top2bcharge;
		float _Top2momentum;
		float _Top2costheta;
		float _Top2bcostheta;
		int _Top2Vtx;
		int _Top2Kaon;
		int _Top2KaonNumber;
		int _Top2KaonCharges[10];
		float _Top2KaonMomentum[10];
		int _Top2leptonCharge;
		int _Top2leptonCorrect;
		float _qCostheta[2];
		int _methodUsed;
		int _methodCorrect;
		float _chiHad;
		float _chiTopMass;
		float _chiTopE;
		float _chiPbstar;
		float _chiCosWb;
		float _chiGammaT;

		int _MCBOscillation;
		int _MCBBarOscillation;
		float _MCTopBangle;
		//float _MCTopBarBangle;
		int _MCLeptonPDG;
		float _MCLeptonMomentum;
		float _MCLeptonCostheta;
		float _MCTopmomentum;
		float _MCTopmass;
		float _MCTopcostheta;
		float _MCTopBarmomentum;
		float _MCTopBarmass;
		float _MCTopBarcostheta;
		float _MCTopcosWb;
		float _MCTopcosWt;
		float _qMCcostheta[2];
		float _MCMass;
		float _qMCBcostheta[2];
		float _MCNeutrinoEnergy;
		float _totalEnergy;
		float _missedEnergy;
		float _gammaT;
		
		void Clear()
		{
			_B1momentum = -1.0;
			_B2momentum = -1.0;
			_B2btag = -1.;
			_B1btag = -1.0;
			_B1costheta = -2.0;
			_B2costheta = -2.0;
			_B1charge = -5.0;
			_B2charge = -5.0;
			_B1truthAngle = -1.0;


			_MCLeptonMomentum = -1.0;
			_MCLeptonCostheta = -2.0;
			_MCLeptonPDG = 0;
			_MCTopmass = -1.0;
			_MCMass = -1.0;
			_MCTopBarmass = -1.0;
			_MCTopmomentum = -1.0;
			_MCTopBarmomentum = -1.0;
			_MCTopcostheta = -2.0;
			_MCTopBarcostheta = -2.0;
			_qMCcostheta[0] = -2.0;
			_qMCcostheta[1] = -2.0;
			_qMCBcostheta[0] = -2.0;
			_qMCBcostheta[1] = -2.0;
			_MCTopBangle = -1.0;
			_MCTopcosWb = -2.0;
			_Top1mass = -1.0;
			_W1mass = -1.0;
			_W2costheta = -2.0;
			_chiHad = -1.0;
			_methodUsed = -1;
			_methodCorrect = -1;
			_qCostheta[0] = -2.0;
			_qCostheta[1] = -2.0;
		}
	};
}
#endif
