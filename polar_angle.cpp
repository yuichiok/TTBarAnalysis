#include <iostream>
#include <string>
#include <TFile.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1.h>
#include <TCut.h>
#include <TEventList.h>

using namespace std ;


#ifdef __CINT__
int polar_angle(){
#else
int main( int argc, char **argv ){
	TApplication app( "app", &argc, argv ) ;
#endif


	string rootfilename = "TTBarProcessorLeftGamma.root" ;
	string treename = "Stats" ;

	TFile file( rootfilename.c_str() ) ;
	if( file.IsZombie() ){
		cout << "cannot open the file" << endl ;
		return 1 ;
	}

	int Top1bcharge, Top2bcharge ;
	float Top1costheta, Top1bcostheta, Top2costheta, Top2bcostheta ;

	int samesignnum = 0, both0num = 0, usednum = 0 ;

	TCut btag = " ( Top1btag > 0.80 ) && ( Top2btag > 0.30 ) " ;
	TCut chi2_1 = " chiTopMass1 + chiTopE1 + chiPbstar1 < 30 " ;
	TCut chi2_2 = " chiTopMass2 + chiTopE2 + chiPbstar2 < 30 " ;
	TCut chi2 = chi2_1 + chi2_2 ;
	TCut kinematic = " ( Top1mass > 140 ) && ( Top1mass < 210 ) && ( Top2mass > 140 ) && ( Top2mass < 210 ) " ;

	TTree* tree = (TTree*)file.Get( treename.c_str() ) ;
	int eventnum = tree->GetEntries() ;
	cout << "eventnum            = " << eventnum << " (100%)" << endl ;
	int afterbtag = tree->GetEntries( btag ) ;
	cout << "after b-tag cut     = " << afterbtag << " (" << (float)100*afterbtag/eventnum << "%)" << endl ;
	int afterkinematic = tree->GetEntries( btag && kinematic ) ;
	cout << "atfer kinematic cut = " << afterkinematic << " (" << (float)100*afterkinematic/eventnum << "%)" << endl;
	int afterchi2 = tree->GetEntries( btag && kinematic && chi2 ) ;
	cout << "after chi2 cut      = " << afterchi2 << " (" << (float)100*afterchi2/eventnum << "%)" << endl;

	TEventList *elist = new TEventList( "elist", "Event List Title" ) ;
	tree->Draw( ">>elist", btag && kinematic && chi2 ) ;
	//afterchi2 = elist->GetN() ;

	tree->SetBranchAddress( "Top1bcharge", &Top1bcharge );
	tree->SetBranchAddress( "Top1costheta", &Top1costheta );
	tree->SetBranchAddress( "Top1bcostheta", &Top1bcostheta );
	tree->SetBranchAddress( "Top2bcharge", &Top2bcharge );
	tree->SetBranchAddress( "Top2costheta", &Top2costheta );
	tree->SetBranchAddress( "Top2bcostheta", &Top2bcostheta );

	TCanvas* c1 = new TCanvas( "c1", "c1" ) ;
	TH1F* htop = new TH1F( "htop", "top polar angle", 40, -1, 1 ) ;
	htop->SetLineColor(2) ;
	TH1F* hbottom = new TH1F( "hbottom", "bottom polar angle", 40, -1, 1 ) ;
	htop->SetLineColor(4) ;

	for( int i=0 ; i<afterchi2 << i++ ){

		tree->GetEntry( elist->GetEntry(i) ) ;

		if( Top1bcharge*Top2bcharge > 0 ){  //the signs of charges are same.
			samesignnum++ ;
		}else if( Top1bcharge==0 && Top2bcharge==0 ){  //both charges are 0.
			both0num++ ;
		}else{  // used for analysis

			if( Top1bcharge < 0 ){  //Top1 = top, Top2 = antitop, Top1b = antibottom, Top2b = bottom
				htop->Fill( Top1costheta ) ;
				htop->Fill( -Top2costheta ) ;
				hbottom->Fill( -Top1bcostheta ) ;
				hbottom->Fill( Top2bcostheta ) ;
			}else if( Top1bcharge > 0 ){  //Top1 = antitop, Top2 = top, Top1b = bottom, Top2b = antibottom
				htop->Fill( -Top1costheta ) ;
				htop->Fill( Top2costheta ) ;
				hbottom->Fill( Top1bcostheta ) ;
				hbottom->Fill( -Top2bcostheta ) ;
			}else if( Top2bcharge > 0 ){  //Top1 = top, Top2 = antitop, Top1b = antibottom, Top2b = bottom
				htop->Fill( Top1costheta ) ;
				htop->Fill( -Top2costheta ) ;
				hbottom->Fill( -Top1bcostheta ) ;
				hbottom->Fill( Top2bcostheta ) ;
			}else{  //Top1 = antitop, Top2 = top, Top1b = bottom, Top2b = antibottom
				htop->Fill( -Top1costheta ) ;
				htop->Fill( Top2costheta ) ;
				hbottom->Fill( Top1bcostheta ) ;
				hbottom->Fill( -Top2bcostheta ) ;
			}

		}

	}

	usednum = afterchi2 - samesignnum - both0num ;
	cout << "used number = " << usednum << endl;
	cout << "same charge sign number = " << samesignnum << endl;
	cout << "both charge 0 number = " << both0num << endl;

	htop->Draw() ;
	hbottom->Draw( "SAME" ) ;

	c1->WaitPrimitive() ;

	return 0;

}




