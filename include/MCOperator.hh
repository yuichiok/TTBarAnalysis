#include <stdlib.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <EVENT/MCParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <EVENT/LCCollection.h>
#include "MathOperator.hh"

#include "marlin/VerbosityLevels.h"
#ifndef _MCOperator_hh
#define _MCOperator_hh
namespace TTbarAnalysis 
{
	class MCOperator 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			MCOperator (EVENT::LCCollection * col);
			virtual ~MCOperator () {};
		//
		//	Methods
		//
			//DO NOT USE THAT ON T-QUARKS!!!
			std::vector< EVENT::MCParticle * > GetPairParticles(int pdg);
			std::vector< EVENT::MCParticle * > GetTopPairParticles(float & topBangle, float & topBarBangle);
			std::vector< EVENT::MCParticle * > GetFinalState();
			EVENT::MCParticle * CombineParticles(EVENT::MCParticle * b, EVENT::MCParticle * w);
			EVENT::MCParticle * FindParticle(int pdg);
			std::vector< EVENT::MCParticle * > GetBquarkPair();
			EVENT::MCParticle * GetNeutrino();

		private:
		//
		//	Data
		//
			EVENT::LCCollection * myCollection;
			std::vector< EVENT::MCParticle * > myBquarkPair;
			EVENT::MCParticle * myNeutrino;
		//
		//	Private methods
		//
	};
}
#endif
