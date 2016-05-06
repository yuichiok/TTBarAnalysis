
#include "MathOperator.hh"
#include "RecoJet.hh"
#include <stdio.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#ifndef Top_Quark_h
#define Top_Quark_h 1
namespace TTbarAnalysis 
{
	struct TopCharge
	{
		TopCharge()
		{
			ByTVCM = NULL;
			ByTrackCount = NULL;
			ByLepton = NULL;
		}
		~TopCharge()
		{
			delete ByTVCM;
			delete ByTrackCount;
			delete ByLepton;
		}

		int * ByTVCM;
		int * ByTrackCount;
		int * ByLepton;
	};
	class TopQuark : public RecoJet
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			TopQuark (RecoJet * b, EVENT::ReconstructedParticle * w);
			TopQuark (RecoJet * b);
			virtual ~TopQuark () 
			{
			};
		//
		//	Methods
		//
			RecoJet * GetB();
			EVENT::ReconstructedParticle * GetW();
			float getMass();
			int GetResultTVCM();
			void SetResultTVCM(int used = 0);
			bool IsHadronic();
			TopCharge & GetComputedCharge();
		private:
		//
		//	Data
		//
			RecoJet * myB;
			EVENT::ReconstructedParticle * myW;
			int myResultTVCM;
			TopCharge myComputedCharge;
		//
		//	Private methods
		//
	};
} /* TTbarAnalisys */
#endif
