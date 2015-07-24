
#include "MathOperator.hh"
#include "RecoJet.hh"
#include <stdio.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#ifndef Top_Quark_h
#define Top_Quark_h 1
namespace TTbarAnalysis 
{
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
			virtual ~TopQuark () {};
		//
		//	Methods
		//
			EVENT::ReconstructedParticle * GetB();
			EVENT::ReconstructedParticle * GetW();
			float getMass();
			bool IsHadronic();
		private:
		//
		//	Data
		//
			EVENT::ReconstructedParticle * myB;
			EVENT::ReconstructedParticle * myW;
		//
		//	Private methods
		//
	};
} /* TTbarAnalisys */
#endif
