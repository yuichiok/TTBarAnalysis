#include <iostream>
#include <string>
#include <TFile.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1.h>
#include <TCut.h>
#include <TEventList.h>
#include <TStyle.h>

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

	gStyle->SetOptStat(0) ;

	int Top1bcharge, Top2bcharge ;
	float Top1costheta, Top1bcostheta, Top2costheta, Top2bcostheta ;
	float MCTopcostheta, MCTopBarcostheta ;
	float qMCBcostheta[2] ;

	TCut btag = " ( Top1btag > 0.80 ) && ( Top2btag > 0.30 ) " ;
	TCut chi2_1 = " chiTopMass1 + chiTopE1 + chiPbstar1 < 30 " ;
	TCut chi2_2 = " chiTopMass2 + chiTopE2 + chiPbstar2 < 30 " ;
	TCut chi2 = chi2_1 + chi2_2 ;
	TCut kinematic = " ( Top1mass > 140 ) && ( Top1mass < 210 ) && ( Top2mass > 140 ) && ( Top2mass < 210 ) " ;
	TCut samecharge = " Top1bcharge * Top2bcharge > 0 " ;
	TCut both0charge = " ( Top1bcharge == 0 ) && ( Top2bcharge == 0 ) " ;

	TTree* tree = (TTree*)file.Get( treename.c_str() ) ;
	int eventnum = tree->GetEntries() ;
	cout << "eventnum            = " << eventnum << " (100%)" << endl ;
	int afterbtag = tree->GetEntries( btag ) ;
	cout << "after b-tag cut     = " << afterbtag << " (" << (float)100*afterbtag/eventnum << "%)" << endl ;
	int afterkinematic = tree->GetEntries( btag && kinematic ) ;
	cout << "atfer kinematic cut = " << afterkinematic << " (" << (float)100*afterkinematic/eventnum << "%)" << endl;
	int afterchi2 = tree->GetEntries( btag && kinematic && chi2 ) ;
	cout << "after chi2 cut      = " << afterchi2 << " (" << (float)100*afterchi2/eventnum << "%)" << endl;

	int samesignnum = tree->GetEntries( btag && kinematic && chi2 && samecharge ) ;
	int both0num = tree->GetEntries( btag && kinematic && chi2 && both0charge ) ;
	int usednum = afterchi2 - samesignnum - both0num ;
	cout << endl ;
	cout << "used number = " << usednum << endl ;
	cout << "same charge sign number = " << samesignnum << endl ;
	cout << "both charge 0 number = " << both0num << endl ;

	TEventList *elist = new TEventList( "elist", "Reconstructed Event List" ) ;
	tree->Draw( ">>elist", btag && kinematic && chi2 && !samecharge && !both0charge ) ;
	//usednum = elist->GetN() ;
	//cout << "used number = " << usednum << endl ;

	tree->SetBranchAddress( "Top1bcharge", &Top1bcharge ) ;
	tree->SetBranchAddress( "Top1costheta", &Top1costheta ) ;
	tree->SetBranchAddress( "Top1bcostheta", &Top1bcostheta ) ;
	tree->SetBranchAddress( "Top2bcharge", &Top2bcharge ) ;
	tree->SetBranchAddress( "Top2costheta", &Top2costheta ) ;
	tree->SetBranchAddress( "Top2bcostheta", &Top2bcostheta ) ;

	TCanvas* c1 = new TCanvas( "c1", "c1", 1280, 480 ) ;
	TH1F* htop = new TH1F( "htop", "top polar angle (reconstructed)", 40, -1, 1 ) ;
	htop->SetLineColor(2) ;
	TH1F* hbottom = new TH1F( "hbottom", "bottom polar angle (reconstructed)", 40, -1, 1 ) ;
	hbottom->SetLineColor(2) ;
	TH1F* htopgen = new TH1F( "htopgen", "top polar angle (generated)", 40, -1, 1 ) ;
	htopgen->SetLineColor(4) ;
	TH1F* hbottomgen = new TH1F( "hbottomgen", "bottom polar angle (generated)", 40, -1, 1 ) ;
	hbottomgen->SetLineColor(4) ;

	double fillweight = (double)1/usednum ;

	for( int i=0 ; i<usednum ; i++ ){

		tree->GetEntry( elist->GetEntry(i) ) ;

		if( Top1bcharge < 0 ){  //Top1 = top, Top2 = antitop, Top1b = bottom, Top2b = antibottom
			htop->Fill( Top1costheta, fillweight ) ;
			htop->Fill( -Top2costheta, fillweight ) ;
			hbottom->Fill( Top1bcostheta, fillweight ) ;
			hbottom->Fill( -Top2bcostheta, fillweight ) ;
		}else if( Top1bcharge > 0 ){  //Top1 = antitop, Top2 = top, Top1b = antibottom, Top2b = bottom
			htop->Fill( -Top1costheta, fillweight ) ;
			htop->Fill( Top2costheta, fillweight ) ;
			hbottom->Fill( -Top1bcostheta, fillweight ) ;
			hbottom->Fill( Top2bcostheta, fillweight ) ;
		}else if( Top2bcharge > 0 ){  //Top1 = top, Top2 = antitop, Top1b = bottom, Top2b = antibottom
			htop->Fill( Top1costheta, fillweight ) ;
			htop->Fill( -Top2costheta, fillweight ) ;
			hbottom->Fill( Top1bcostheta, fillweight ) ;
			hbottom->Fill( -Top2bcostheta, fillweight ) ;
		}else{  //Top1 = antitop, Top2 = top, Top1b = antibottom, Top2b = bottom
			htop->Fill( -Top1costheta, fillweight ) ;
			htop->Fill( Top2costheta, fillweight ) ;
			hbottom->Fill( -Top1bcostheta, fillweight ) ;
			hbottom->Fill( Top2bcostheta, fillweight ) ;
		}

	}

	tree->SetBranchStatus( "Top*", 0 ) ;
	tree->SetBranchAddress( "MCTopcostheta", &MCTopcostheta ) ;
	tree->SetBranchAddress( "MCTopBarcostheta", &MCTopBarcostheta ) ;
	tree->SetBranchAddress( "qMCBcostheta", &qMCBcostheta ) ;

	TEventList* elistgen = new TEventList( "elistgen", "Generated Event List" ) ;
	tree->Draw( ">>elistgen", " ( MCTopcostheta != -2 ) && ( MCTopBarcostheta != -2 ) && ( qMCBcostheta[0] != -2 ) && ( qMCBcostheta[1] != -2 ) " ) ;
	int usedgennum = elistgen->GetN() ;
	cout << endl ;
	cout << "used genarated number = " << usedgennum << endl ;

	fillweight = (double)1/usedgennum ;

	for( int i=0 ; i<usedgennum ; i++ ){

		tree->GetEntry( elistgen->GetEntry(i) ) ;

		htopgen->Fill( MCTopcostheta, fillweight ) ;
		htopgen->Fill( -MCTopBarcostheta, fillweight ) ;
		hbottomgen->Fill( qMCBcostheta[0], fillweight ) ;
		hbottomgen->Fill( -qMCBcostheta[1], fillweight ) ;

	}

	c1->Divide(2,1) ;
	c1->cd(1) ;
	htopgen->SetTitle( "top polar angle (Rec:Red Gen:Blue)" ) ;
	htopgen->Draw() ;
	htop->Draw( "SAME" ) ;

	c1->cd(2) ;
	hbottomgen->SetTitle( "bottom polar angle (Rec:Red Gen:Blue)" ) ;
	hbottomgen->Draw() ;
	hbottom->Draw( "SAME" ) ;

	c1->WaitPrimitive() ;

//	c1->Print( "picture/newpicture.png" ) ;


	return 0;

}




