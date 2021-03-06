#include "QQBarMCOperator.hh"
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::MCParticle;
using IMPL::MCParticleImpl;
namespace TTBarProcessor
{
	QQBarMCOperator:: QQBarMCOperator (LCCollection * col)
	{
		myCollection = col;
		myNeutrino = NULL;
	}
	//DO NOT USE THAT ON T-QUARKS!!!
	vector< MCParticle * > QQBarMCOperator::GetPairParticles(int pdg)
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
			if (particle->getPDG() == pdg&& !b)// &&  countParticle == 0) 
			{
				b = particle;
			}
			if (particle->getPDG() == -pdg  && !bbar)//&& countAntiparticle == 0) 
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
	vector< MCParticle * > QQBarMCOperator::GetTopPairParticles(float & topBangle, float & topcosWb)
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
			int finalsize = final.size();
			for (int i = 0; i < finalsize; i++) 
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
		myBquarkPair.push_back(b);
		myBquarkPair.push_back(bbar);
		myWPair.push_back(wplus);
		myWPair.push_back(wminus);
		MCParticle * top = CombineParticles(b, wplus);
		topBangle = MathOperator::getAngle(top->getMomentum(), b->getMomentum());
		MCParticle * topbar = CombineParticles(bbar, wminus);
		topcosWb = std::cos( MathOperator::getAngle(wplus->getMomentum(), b->getMomentum()));
		pair.push_back(top);
		pair.push_back(topbar);
		return pair;
	}
	vector< MCParticle * > QQBarMCOperator::GetWPair()
	{
		return myWPair;
	}
	vector< MCParticle * > QQBarMCOperator::GetBquarkPair()
	{
		return myBquarkPair;
	}
	MCParticle * QQBarMCOperator::CombineParticles(EVENT::MCParticle * b, EVENT::MCParticle * w)
	{
		MCParticleImpl * result = new MCParticleImpl();
		double energy = b->getEnergy() + w->getEnergy();
		double momentum[3];
		for (int i = 0; i < 3; i++) 
		{
			momentum[i] = b->getMomentum()[i] + w->getMomentum()[i];
		}
		double charge = b->getCharge() + w->getCharge();
		if (abs(b->getPDG())< 5 && abs(w->getPDG()) < 5) 
		{
			result->setGeneratorStatus(1);
		}
		else 
		{
			result->setGeneratorStatus(0);
		}
		double mass = std::sqrt(energy * energy - momentum[0] * momentum[0] - momentum[1]*momentum[1] - momentum[2] * momentum[2] );
		int pdg = charge / std::abs(charge) * 6;
		//result->setEnergy(energy);
		result->setMass(mass);
		result->setMomentum(momentum);
		result->setCharge(charge);
		result->setPDG(pdg);
		return result;
	}
	vector <MCParticle *> QQBarMCOperator::GetFinalStateBkg()
	{
		vector <MCParticle *> result;
		/*for (unsigned int i = 0; i < 4; i++) 
		{
			MCParticle * particle = dynamic_cast<MCParticle*>( myCollection->getElementAt(i+2) );
			result.push_back(particle);
		}*/
		result.push_back(dynamic_cast<MCParticle*>( myCollection->getElementAt( 2)));
		result.push_back(dynamic_cast<MCParticle*>( myCollection->getElementAt( 3)));

		MCParticle * higgs = FindParticle(25);
		if (higgs->getDaughters().size() == 2) 
		{
			result.push_back(higgs->getDaughters()[0]);
			result.push_back(higgs->getDaughters()[1]);
		}
		if (higgs->getDaughters().size() == 1) 
		{
			if (higgs->getDaughters()[0]->getDaughters().size() == 2) 
			{
				result.push_back(higgs->getDaughters()[0]->getDaughters()[0]);
				result.push_back(higgs->getDaughters()[0]->getDaughters()[1]);

			}
		}
		return result;
	}
	vector <MCParticle *> QQBarMCOperator::GetFinalState()
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
				if (abs(daughters[i]->getPDG()) == 12 ||
				    abs(daughters[i]->getPDG()) == 14 ||	
				    abs(daughters[i]->getPDG()) == 16) 
				{
					myNeutrino = daughters[i];
				}
			}
		}
		else 
		{
			std::cout << "ERROR: Particle has " << particle->getPDG() << ", e- not found!\n";
		}
		return result;
	}
	MCParticle * QQBarMCOperator::FindParticle(int pdg)
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
	MCParticle * QQBarMCOperator::GetNeutrino()
	{
		return myNeutrino;
	}
	MCParticle * QQBarMCOperator::GetTauLepton()
	{
		MCParticle * tau = FindParticle(15);
		MCParticle * taubar = FindParticle(-15);
		MCParticle * particle = (tau)? tau : taubar;
		if (!particle) 
		{
			return NULL;
		}
		MCParticle * result = NULL;
		if (particle->getDaughters().size() == 1 && abs(particle->getDaughters()[0]->getPDG()) == 15) 
		{
			particle = particle->getDaughters()[0];
		}
		if (particle->getDaughters().size() == 3) 
		{
			for (unsigned int i = 0; i < particle->getDaughters().size(); i++) 
			{
				if (abs(particle->getDaughters()[i]->getPDG()) == 13) 
				{
					std::cout << "Particle tau found!\n";
					result = particle->getDaughters()[i];
					break;
				}
			}
		}
		return result;
	}
}
