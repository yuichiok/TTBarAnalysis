#include "VertexChargeOperator.hh"
using EVENT::ReconstructedParticle;
using UTIL::LCRelationNavigator;
using EVENT::MCParticle;
using EVENT::Vertex;
using EVENT::LCObject;
using std::vector;
using EVENT::LCCollection;
namespace TTbarAnalysis 
{
	const float VertexChargeOperator::BMASS = 5.279;
	const float VertexChargeOperator::CTAU = 0.4554;
	VertexChargeOperator:: VertexChargeOperator(LCCollection * rel)
	{
		myResultingB = 0;
		myRelCollection = rel;
	}
	int VertexChargeOperator::GetResultingB()
	{
		return myResultingB;
	}
	float VertexChargeOperator::checkAsymmetry(TopQuark * top, Vertex * vertex1,ReconstructedParticle * kaon1, float & top1costheta)
	{
		float __mccharge = top->GetB()->__GetMCCharge();
		float bmomentum = top->GetB()->GetHadronMomentum() * 2;
		float gamma = std::sqrt(bmomentum + BMASS) / BMASS;
		float distance = MathOperator::getModule(vertex1->getPosition());
		int ntracks = top->GetB()->GetNumberOfVertexParticles();
		float coefficient = distance / CTAU;
		std::cout << "We have a kaon " << kaon1->getType() 
			  << " with charge " << kaon1->getCharge() 
			  << " from b with charge of " << __mccharge
			  << " coefficient: " << coefficient
			  << "\n";
		float result = (kaon1->getCharge() < 0.0) ? top1costheta : -top1costheta;
		if (__mccharge * kaon1->getCharge() < 0) 
		{
			top->SetResultTVCM(-1);
			std::cout << "FATAL: Charges are different!!!\n";
			std::cout << "N particles: " << top->GetB()->GetNumberOfVertexParticles()
				  << " Ternary charge: " << vertex1->getAssociatedParticle()->getCharge()
				  << " distance: " << MathOperator::getModule(vertex1->getPosition())
				  << "\n";
		}
		else 
		{
			top->SetResultTVCM(1);
		}
		return result;
	}

	int VertexChargeOperator::CountKaons(TopQuark * top1, TopQuark * top2)
	{
		Vertex * vertex1 =  getTernaryVertex(top1);
		Vertex * vertex2 =  getTernaryVertex(top2);
		ReconstructedParticle * kaon1 = __getKaonCheat(vertex1);
		ReconstructedParticle * kaon2 = __getKaonCheat(vertex2);
		if (kaon1 || kaon2) 
		{
			return 1;
		}
		return 0;
	}
	float VertexChargeOperator::GetAsymmetryTVCM(TopQuark * top1, TopQuark * top2)
	{
		std::cout << "\n\n----TVCM method is called----\n";
		Vertex * vertex1 =  getTernaryVertex(top1);
		Vertex * vertex2 =  getTernaryVertex(top2);
		ReconstructedParticle * kaon1 = __getKaonCheat(vertex1);
		ReconstructedParticle * kaon2 = __getKaonCheat(vertex2);

		if (!kaon1 && !kaon2) 
		{
			std::cout << "We have no kaons!\n";
			return -2.0;
		}
		vector<float> direction = MathOperator::getDirection(top1->getMomentum());
		float top1costheta =  std::cos( MathOperator::getAngles(direction)[1] );
		float result = -2.0;
		if (kaon1 && !kaon2) 
		{
			checkAsymmetry(top1, vertex1, kaon1, top1costheta);
			result = (kaon1->getCharge() < 0.0) ? top1costheta : -top1costheta;
			myResultingB = kaon1->getCharge();
			TopCharge & topCharge = top1->GetComputedCharge();
			topCharge.ByTVCM = new int(kaon1->getCharge());
		}
		if (kaon2 && !kaon1) 
		{
			checkAsymmetry(top2, vertex2, kaon2, top1costheta);
			result = (kaon2->getCharge() < 0.0) ? -top1costheta : top1costheta;
			myResultingB = -kaon2->getCharge();
			TopCharge & topCharge = top2->GetComputedCharge();
			topCharge.ByTVCM = new int(kaon2->getCharge());
		}
		if (kaon2 && kaon1) 
		{
			std::cout << "We have twwo kaons!\n";
			if (kaon1->getCharge() * kaon2->getCharge() < 0) 
			{
				result = (kaon1->getCharge() < 0.0) ? top1costheta : -top1costheta;
			}
			TopCharge & topCharge1 = top1->GetComputedCharge();
			topCharge1.ByTVCM = new int(kaon1->getCharge());
			TopCharge & topCharge2 = top2->GetComputedCharge();
			topCharge2.ByTVCM = new int(kaon2->getCharge());
			/*if ( MathOperator::getModule(vertex1->getPosition()) < MathOperator::getModule(vertex2->getPosition())) 
			{
				result = (kaon1->getCharge() < 0.0) ? top1costheta : -top1costheta;
				float __mccharge = top->GetB()->__GetMCCharge();
				if (__mccharge * kaon1->getCharge() < 0.0) 
				{
					top->SetResultTVCM(-1);
				}
				else 
				{
					top->SetResultTVCM(1);
				}
			}
			else 
			{
				result = (kaon2->getCharge() < 0.0) ? -top1costheta : top1costheta;
				float __mccharge = top2->GetB()->__GetMCCharge();
				if (__mccharge * kaon2->getCharge() < 0.0) 
				{
					top2->SetResultTVCM(-1);
				}
				else 
				{
					top2->SetResultTVCM(1);
				}
			}*/
		}
		/*if (leptonkaon1) 
		{
			std::cout << "We have a lepton " << leptonkaon1->getType() << " with charge " << leptonkaon1->getCharge() << "\n";
			result = (leptonkaon1->getCharge() < 0.0) ? -top1costheta : top1costheta;  
		}
		if (leptonkaon2) 
		{
			std::cout << "We have a lepton " << leptonkaon2->getType() << " with charge " << leptonkaon2->getCharge() << "\n";
			result = (leptonkaon2->getCharge() < 0.0) ? top1costheta : -top1costheta;  
		}*/
		
		return result;
	}
	Vertex * VertexChargeOperator::getTernaryVertex(TopQuark * top)
	{
		vector< Vertex * > * vertices = top->GetRecoVertices();
		if (!vertices || vertices->size() == 0) 
		{
			return NULL;
		}
		if (vertices->size() == 1) 
		{
			std::cout << "Only one vertex is reconstructed!\n";
			return NULL;
		}
		float maxdistance = 0.0;
		Vertex * winner = NULL;
		for (int i = 0; i < vertices->size(); i++) 
		{
			Vertex * vertex  = vertices->at(i);
			float distance = MathOperator::getModule(vertex->getPosition());
			if (distance > maxdistance) 
			{
				winner = vertex;
				maxdistance = distance;
			}
		}
		return winner;
	}
	ReconstructedParticle * VertexChargeOperator::getFlavourParticle(Vertex * ternary)
	{
		if (!ternary) 
		{
			return NULL;
		}
		const vector<ReconstructedParticle *> particles = ternary->getAssociatedParticle()->getParticles();
		ReconstructedParticle * result = NULL;
		int count = 0;
		if ( particles.size() > 1) 
		{
			//return NULL;
		}
		for (int i = 0; i < particles.size(); i++) 
		{
			ReconstructedParticle * particle = particles[i];
			if (abs(particle->getType()) == 11 ||      //electron
			    abs(particle->getType()) == 13) 	   //muon
			{
				result = particle;
				count++;
			}
		}
		if (count>1) 
		{
			std::cout << "ERROR: Multiple leptons found!\n";
			return NULL;
		}
		return result;
	}
	ReconstructedParticle * VertexChargeOperator::__getKaonCheat(Vertex * ternary)
	{
		if (!ternary) 
		{
			return NULL;
		}
		LCRelationNavigator navigator(myRelCollection);
		const vector<ReconstructedParticle *> particles = ternary->getAssociatedParticle()->getParticles();
		ReconstructedParticle * result = NULL;
		int count = 0;
		for (int i = 0; i < particles.size(); i++) 
		{
			ReconstructedParticle * particle = particles[i];
			std::cout << "Charge: " << particle->getCharge() << "\n";
			vector<float> direction = MathOperator::getDirection(particle->getMomentum());
			int tpchits = particle->getTracks()[0]->getSubdetectorHitNumbers()[6];
			float costheta =  std::abs(std::cos( MathOperator::getAngles(direction)[1] ));
			if (costheta > 0.95 || tpchits < 60) 
			{
				continue;
			}
			vector< LCObject * > obj = navigator.getRelatedToObjects(particle);
			vector< float > weights = navigator.getRelatedToWeights(particle);
			MCParticle * winner = NULL;
			float maxweight = 0.50;
			for (unsigned int i = 0; i < obj.size(); i++) 
			{
				if (weights[i] > maxweight) 
				{
					winner = dynamic_cast< MCParticle * >(obj[i]);
					maxweight = weights[i];
				}
			}
			if (!winner) 
			{
				std::cout << "ERROR: no genparticle!\n";
			}
			if (winner && abs(winner->getPDG()) == 321) 
			{
				result = particle;
				count++;
			}
		}
		if (count > 1) 
		{
			std::cout << "ERROR: Multiple kaons found!\n";
			return NULL;
		}
		return result;
	}

} /* TTbarAnalysis */
