#include <vector>
#include "MathOperator.hh"
#include "TopQuark.hh"
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Vertex.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCObject.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/LCCollection.h>
namespace TTbarAnalysis 
{
	class VertexChargeOperator 
	{
		public:
		//
		//	Constants
		//
			static const float BMASS;// = 5.279; //GeV
			static const float CTAU;// = 0.4554; 
		//
		//	Constructors
		//
			VertexChargeOperator (EVENT::LCCollection * rel = NULL) ;
			virtual ~VertexChargeOperator () {};
		//
		//	Methods
		//
			float GetAsymmetryTVCM(TopQuark * top, TopQuark * top2); //TernaryVertexChargeMeasurement
			int GetResultingB();
			int CountKaons(TopQuark * top, TopQuark * top2);
		private:
		//
		//	Data
		//
			EVENT::LCCollection * myRelCollection;
			int myResultingB;
			
		//
		//	Private methods
		//
			EVENT::Vertex * getTernaryVertex(TopQuark * top);
			EVENT::ReconstructedParticle * getFlavourParticle(EVENT::Vertex * ternary);
			EVENT::ReconstructedParticle * __getKaonCheat(EVENT::Vertex * ternary);
			float checkAsymmetry(TopQuark * top, EVENT::Vertex * vertex1, EVENT::ReconstructedParticle * kaon1, float & top1costheta);
	};
} /* TTbarAnalysis */
