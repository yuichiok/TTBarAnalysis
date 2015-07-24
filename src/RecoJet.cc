
#include "RecoJet.hh"
using EVENT::ReconstructedParticle;
using IMPL::ReconstructedParticleImpl;
using std::vector;
namespace TTbarAnalysis
{
	RecoJet::RecoJet (ReconstructedParticle * rawjet, float btag,float ctag, int number)
	//: IMPL::ReconstructedParticleImpl(rawjet)
	{
		myBTag = btag;
		myCTag = ctag;
		myNumber = number;
		myMCPDG = 0;
		myRecoVertices = NULL;
		myRawRecoJet = rawjet;
		setMomentum(rawjet->getMomentum());

		setMass(rawjet->getMass());
		setEnergy(rawjet->getEnergy());
		_particles.reserve(rawjet->getParticles().size());
		_particles.insert(_particles.begin(),rawjet->getParticles().begin(),rawjet->getParticles().end());

	}
	RecoJet::RecoJet ()
	{
		myRecoVertices = NULL;
		myMCPDG = 0;
		
	}
	void RecoJet::SetRecoVertices(std::vector<  EVENT::Vertex * > * vertices)
	{
		myRecoVertices = vertices;
	}
	std::vector<  EVENT::Vertex * > * RecoJet::GetRecoVertices()
	{
		return myRecoVertices;
	}
	void RecoJet::SetBTag(float value)
	{
		myBTag = value;
	}
	void RecoJet::SetCTag(float value)
	{
		myCTag = value;
	}
	float RecoJet::GetBTag()
	{
		return myBTag;
	}
	float RecoJet::GetCTag()
	{
		return myCTag;
	}
	int RecoJet::GetNumberOfVertices()
	{
		if (myRecoVertices) 
		{
			return myRecoVertices->size();
		}
		return myNumber;
	}
	int RecoJet::GetNumberOfVertexParticles()
	{
		int sum = -1;
		if (myRecoVertices) 
		{
			sum = 0;
			for (unsigned int i = 0; i < myRecoVertices->size(); i++) 
			{
				sum += myRecoVertices->at(i)->getAssociatedParticle()->getParticles().size();
			}
		}
		return sum;
	}
	float RecoJet::GetHadronCharge()
	{
		float charge = -5.0;
		if (myRecoVertices) 
		{
			charge = 0.0;
			for (unsigned int i = 0; i < myRecoVertices->size(); i++) 
			{
				charge += myRecoVertices->at(i)->getAssociatedParticle()->getCharge();
			}
		}
		return charge;
	}
	float RecoJet::GetHadronMomentum()
	{
		float momentum = -1.0;
		if (myRecoVertices) 
		{
			momentum = 0.0;
			for (unsigned int i = 0; i < myRecoVertices->size(); i++)
			{
				momentum += MathOperator::getModule(myRecoVertices->at(i)->getAssociatedParticle()->getMomentum()); // CRUNCH!!!
			}
		}
		return momentum;
	}
	float RecoJet::GetHadronMass()
	{
		float mass = -1.0;
		if (myRecoVertices) 
		{
			 mass = 0.0;
			 for (int i = 0; i < myRecoVertices->size(); i++) 
			 {
			 	mass += myRecoVertices->at(i)->getAssociatedParticle()->getMass();
			 }
		}
		return mass;
	}
	
	int RecoJet::GetMCPDG()
	{
		return myMCPDG;
	}
	void RecoJet::SetMCPDG(int pdg)
	{
		myMCPDG = pdg;
	}

}
