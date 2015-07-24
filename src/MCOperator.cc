#include "MCOperator.hh"
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::MCParticle;
using IMPL::MCParticleImpl;
namespace TTbarAnalysis
{
	MCOperator:: MCOperator (LCCollection * col)
	{
		myCollection = col;
	}
	//DO NOT USE THAT ON T-QUARKS!!!
	vector< MCParticle * > MCOperator::GetPairParticles(int pdg)
	{
		pdg = abs(pdg);
		vector< MCParticle * > pair;
		if (pdg < 1) 
		{
			return pair;
		}
		int number = myCollection->getNumberOfElements();
		MCParticle * b = NULL;
		MCParticle * bbar = NULL;
		for (int i = 0; i < number; i++) 
		{
			MCParticle * particle = dynamic_cast<MCParticle*>( myCollection->getElementAt(i) );
			if (particle->getPDG() == pdg)// && countParticle == 0) 
			{
				b = particle;
			}
			if (particle->getPDG() == -pdg)// && countAntiparticle == 0) 
			{
				bbar =  particle;
			}
		}
		if (b) 
		{
			pair.push_back(b);
		}
		if (bbar) 
		{
			pair.push_back(bbar);
		}
		return pair;
	}
	vector< MCParticle * > MCOperator::GetTopPairParticles()
	{
		vector< MCParticle * > pair;
		MCParticle * b = FindParticle(5);
		MCParticle * bbar = FindParticle(-5);
		MCParticle * wplus = FindParticle(24);
		MCParticle * wminus = FindParticle(-24);
		if (!b || !bbar ) 
		{
			return pair;
		}
		if (!wplus || !wminus) 
		{
			vector< MCParticle * > final = GetFinalState();
			std::cout << "Nfinal - bbbar: " << final.size() << "\n";
			if (final.size() < 4) 
			{
				std::cout << "CRUNCH\n";
				return pair; // CRUNCH
			}
			for (int i = 0; i < final.size(); i++) 
			{
				for (int j = 0; j < i; j++) 
				{
					if (i == j) 
					{
						continue;
					}
					MCParticle * candidate = CombineParticles(final[i],final[j]);
					if (std::abs(candidate->getCharge()) > 0.99 && std::abs(candidate->getCharge()) < 1.01 ) 
					{
						if (candidate->getCharge() > 0) 
						{
							wplus = candidate;
							std::cout << "W+ found: q " << candidate->getCharge() << " m " << candidate->getMass() << "\n";
						}
						else 
						{
							wminus = candidate;
							std::cout << "W- found: q " << candidate->getCharge() << " m " << candidate->getMass() << "\n";
						}
					}
				}
			}
		}
		MCParticle * top = CombineParticles(b, wplus);
		MCParticle * topbar = CombineParticles(bbar, wminus);
		pair.push_back(top);
		pair.push_back(topbar);
		return pair;
	}
	MCParticle * MCOperator::CombineParticles(EVENT::MCParticle * b, EVENT::MCParticle * w)
	{
		MCParticleImpl * result = new MCParticleImpl();
		double energy = b->getEnergy() + w->getEnergy();
		double momentum[3];
		for (int i = 0; i < 3; i++) 
		{
			momentum[i] = b->getMomentum()[i] + w->getMomentum()[i];
		}
		double charge = b->getCharge() + w->getCharge();
		double mass = std::sqrt(energy * energy - momentum[0] * momentum[0] - momentum[1]*momentum[1] - momentum[2] * momentum[2] );
		int pdg = charge / std::abs(charge) * 6;
		//result->setEnergy(energy);
		result->setMass(mass);
		result->setMomentum(momentum);
		result->setCharge(charge);
		result->setPDG(pdg);
		return result;
	}

	vector <MCParticle *> MCOperator::GetFinalState()
	{
		vector <MCParticle *> result;
		MCParticle * particle = dynamic_cast<MCParticle*>( myCollection->getElementAt(2) );
		if (abs(particle->getPDG()) == 11) 
		{
			const vector< MCParticle * > daughters = particle->getDaughters();
			for (unsigned int i = 0; i < daughters.size(); i++) 
			{
				if (abs(daughters[i]->getPDG()) != 5) 
				{
					result.push_back(daughters[i]);
				}
			}
		}
		else 
		{
			std::cout << "ERROR: Particle has " << particle->getPDG() << ", e- not found!\n";
		}
		return result;
	}
	MCParticle * MCOperator::FindParticle(int pdg)
	{
		int number = myCollection->getNumberOfElements();
		MCParticle * result = NULL;
		for (int i = 0; i < number; i++) 
		{
			MCParticle * particle = dynamic_cast<MCParticle*>( myCollection->getElementAt(i) );
			if (particle->getPDG() == pdg) 
			{
				result = particle;
				break;
			}
		}
		if (!result) 
		{
			std::cout << "Particle " << pdg << " not found!\n";
		}
		return result;
	}
}
