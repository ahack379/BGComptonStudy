
void cosmics(){

	TFile * _file0 = new TFile("sixCosmicFiles.root","READ") ;
	TTree * t = ( TTree *)_file0->Get("ana_tree");

	//Photon Conversions (no fiducial cut): All Comptons and pair productions in the detector, parent gamma out of detector
	TH2D * h = new TH2D("h","Vertices of Comptons/PP in Detector, Ancestor Out",100,0,256.35,100,-116.5,116.5) ;
	TCanvas * c = new TCanvas("c","c");
	t->Draw("_Y:_X>>h","(_PDG==3 || _PDG ==4) && _inActiveVolume && !_ancestorInActiveVolume","COLZ");
	h->GetXaxis()->SetTitle("X (cm)");
	h->GetYaxis()->SetTitle("Y (cm)");

	TH2D * h0 = new TH2D("h0","Vertices of Comptons/PP in Detector, Ancestor In",100,0,256.35,100,-116.5,116.5) ;
	TCanvas * c0 = new TCanvas("c0","c0");
	t->Draw("_Y:_X>>h","(_PDG==3 || _PDG ==4) && _inActiveVolume && _ancestorInActiveVolume","COLZ");
	h0->GetXaxis()->SetTitle("X (cm)");
	h0->GetYaxis()->SetTitle("Y (cm)");

	TH2D * h5 = new TH2D("h5","Vertices of Comptons in Detector with Cyl cut(10cm), BDtW cut(50cm)",100,0,256.35,100,-116.5,116.5) ;
	TCanvas * c5 = new TCanvas("c5","c5");
	t->Draw("_Y:_X>>h","_PDG==3  && _inActiveVolume && _minMuDist > 10 && _distBackAlongTraj > 50","COLZ");
	h5->GetXaxis()->SetTitle("X (cm)");
	h5->GetYaxis()->SetTitle("Y (cm)");

	TH2D * h6 = new TH2D("h6","Vertices of Comptons in Detector with Cyl cut(10cm)",100,0,256.35,100,-116.5,116.5) ;
	TCanvas * c6 = new TCanvas("c6","c6");
	t->Draw("_Y:_X>>h","_PDG==3  && _inActiveVolume && _minMuDist > 10","COLZ");
	h6->GetXaxis()->SetTitle("X (cm)");
	h6->GetYaxis()->SetTitle("Y (cm)");




////******ADD the YX verticies plots for entries with all cuts included.

	//Energy things, Compton and PP, parent out
	TH1D * h1 = new TH1D("h1","Energy of Comptons in Detector, Parent Out",100,0,0.2) ;
	TCanvas * c1 = new TCanvas("c1","c1");
    t->Draw("_E>>h1","_PDG==3 && _inActiveVolume && !_ancestorInActiveVolume") ;
	c1->SetLogy() ;
	h1->GetXaxis()->SetTitle("E (GeV)");

	TH1D * h2 = new TH1D("h2","Energy of PP/Compton in Detector, Ancestor In",100,0,1) ;
	TCanvas * c2 = new TCanvas("c2","c2");
    t->Draw("_E>>h2","_PDG ==4 && _inActiveVolume && _ancestorInActiveVolume") ;
    t->Draw("_E","_PDG==3 && _inActiveVolume && _ancestorInActiveVolume","same") ;
	c2->SetLogy() ;
	h2->GetXaxis()->SetTitle("E (GeV)");
    h2->SetLineColor(kRed);

	TH1D * h2b = new TH1D("h2b","Energy of PP/Compton in Detector, Ancestor Out",100,0,1) ;
	TCanvas * c2b = new TCanvas("c2b","c2b");
    t->Draw("_E>>h2","_PDG==4 && _inActiveVolume && !_ancestorInActiveVolume") ;
    t->Draw("_E","_PDG==3 && _inActiveVolume && !_ancestorInActiveVolume","same") ;
	c2b->SetLogy() ;
	h2b->GetXaxis()->SetTitle("E (GeV)");
    h2b->SetLineColor(kRed);


	TH1D * h3 = new TH1D("h3","Distribution of e- Distance to Closest Track",500,0,107);
	TCanvas * c3 = new TCanvas("c3","c3");
	t->Draw("_minMuDist");
	c3->SetLogy();
	h3->GetXaxis()->SetTitle("Distance [cm]");
	h3->GetYaxis()->SetTitle("Count");

	TH1D * h4 = new TH1D("h4","Impact Parameter to Closest Track in Event - e-",500,0,100);
	TCanvas * c4 = new TCanvas("c4","c4");
	t->Draw("_minMuIP");
	c4->SetLogy();
	h4->GetXaxis()->SetTitle("Minimum Impact Parameter [cm]");
	h4->GetYaxis()->SetTitle("Count");

	//TLegend *l = new TLegend(0.6,0.6,0.7,0.75) ;
	//l->AddEntry(h10,"Z, No Energy Cut","fl");   // h1 and h2 are histogram pointers
	//l->Draw();

} //end testScript.C
