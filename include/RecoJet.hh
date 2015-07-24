#include "MathOperator.hh"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/Vertex.h>

#ifndef _RecoJet_hh
#define _RecoJet_hh
namespace TTbarAnalysis 
{
	class RecoJet : public IMPL::ReconstructedParticleImpl
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			RecoJet (EVENT::ReconstructedParticle * rawjet, float btag, float ctag, int number);// : IMPL::ReconstructedParticleImpl(rawjet);
			RecoJet ();
			virtual ~RecoJet () {};
		//
		//	Methods
		//
			float GetBTag();
			float GetCTag();
			void SetBTag(float value);
			void SetCTag(float value);
			void SetRecoVertices(std::vector<  EVENT::Vertex * > * vertices);
			std::vector<  EVENT::Vertex * > * GetRecoVertices();
			int GetNumberOfVertices();
			int GetNumberOfVertexParticles();
			float GetHadronCharge();
			float GetHadronMomentum();
			float GetHadronMass();
			const double * GetMomentum();
			int GetMCPDG();
			void SetMCPDG(int pdg);
			
		protected:
		//
		//	Data
		//
			 float myBTag;
			 float myCTag;
			 int myNumber;
			 int myMCPDG;
			 EVENT::ReconstructedParticle * myRawRecoJet;
			 std::vector<  EVENT::Vertex * > * myRecoVertices;
		//
		//	Private methods
		//
	};
}
#endif
