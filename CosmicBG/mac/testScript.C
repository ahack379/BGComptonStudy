
void testScript(){

	TFile * _file0 = new TFile("sixCosmicFiles.root","READ") ;
	TTree * t = ( TTree *)_file0->Get("ana_tree");

	//Photon Conversions (no fiducial cut): All Comptons and pair productions in the detector, parent gamma out of detector
	TH2D * h = new TH2D("h","Vertices of Comptons in Detector, Parent Out",50,0,256.35,50,-116.5,116.5) ;
	TCanvas * c = new TCanvas("c","c");
	t->Draw("_Y:_X>>h","_PDG==3 && _distBackAlongTraj >0 && (_parentX <0 || _parentX >256.35 || _parentY <-116.5 || _parentY > 116.5 ||  _parentZ <0 || _parentZ > 1036.8)","COLZ");
	h->GetXaxis()->SetTitle("X (cm)");
	h->GetYaxis()->SetTitle("Y (cm)");

	TH2D * h0 = new TH2D("h0","Vertices of PP in Detector, Parent Out",20,0,256.35,20,-116.5,116.5) ;
	TCanvas * c0 = new TCanvas("c0","c0");
	t->Draw("_Y:_X>>h0"," _PDG ==4 && _distBackAlongTraj >0 && (_parentX <0 || _parentX >256.35 || _parentY <-116.5 || _parentY > 116.5 ||  _parentZ <0 || _parentZ > 1036.8)","COLZ");
	h0->GetXaxis()->SetTitle("X (cm)");
	h0->GetYaxis()->SetTitle("Y (cm)");



	//Energy things, Compton and PP, parent out
	TH1D * h1 = new TH1D("h1","Energy of Comptons in Detector, Parent Out",100,0,0.2) ;
	TCanvas * c1 = new TCanvas("c1","c1");
    t->Draw("_E>>h1","_PDG==3 && _distBackAlongTraj >0 && (_parentX <0 || _parentX >256 || _parentY <-116 || _parentY >116 ||  _parentZ <0 || _parentZ > 1037)");
	c1->SetLogy() ;
	h1->GetXaxis()->SetTitle("E (GeV)");
	

	TH1D * h2 = new TH1D("h2","Energy of PP in Detector, Parent Out",100,0,1) ;
	TCanvas * c2 = new TCanvas("c2","c2");
    t->Draw("_E>>h2","_PDG==4 && _distBackAlongTraj >0 && (_parentX <0 || _parentX >256 || _parentY <-116 || _parentY >116 ||  _parentZ <0 || _parentZ > 1037)");
	c2->SetLogy() ;
	h2->GetXaxis()->SetTitle("E (GeV)");



	//Distance back to wall, Compton, parent out 
	TH1D * h3 = new TH1D("h3","Distance Back to Wall for Compton, Parent Out",100,0,1036.35) ;
	TCanvas * c3 = new TCanvas("c3","c3") ;	
    t->Draw("_distBackAlongTraj>>h3","_PDG ==3 && _distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)");
    t->Draw("_distBackAlongTraj","_PDG ==3 && E > 0.05 &&_distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)","same");
	c3->SetLogy() ;
	h3->GetXaxis()->SetTitle("Distance (cm)");
	h3->SetLineColor(kRed);

	TH1D * h4 = new TH1D("h4","Distance Back to Wall for PP, Parent Out",50,0,600) ;
	TCanvas * c4 = new TCanvas("c4","c4") ;	
    t->Draw("_distBackAlongTraj>>h4","_PDG ==4 && _distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)");
    t->Draw("_distBackAlongTraj","_PDG ==4 && E > 0.05 &&_distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)","same");
	c4->SetLogy() ;
	h4->GetXaxis()->SetTitle("Distance (cm)");
	h4->SetLineColor(kRed);



	//Vertices of Compton, Parent Out
	TH1D * h5 = new TH1D("h5","Compton Vertex (X), Parent Out",100,0,256.35) ;
	TCanvas * c5 = new TCanvas("c5","c5") ;	
    t->Draw("_X>>h5","_PDG == 3 && _distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)");
    t->Draw("_X","_PDG == 3 && _E > 0.005 &&_distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)","same");
	h5->GetXaxis()->SetTitle("Distance (cm)");
	h5->SetLineColor(kRed);

	TH1D * h6 = new TH1D("h6","Compton Vertex (Y), Parent Out",150,-116.5,116.5) ;
	TCanvas * c6 = new TCanvas("c6","c6") ;	
    t->Draw("_Y>>h6","_PDG == 3 && _distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)");
    t->Draw("_Y","_PDG == 3 && _E > 0.005 &&_distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)","same");
	h6->GetXaxis()->SetTitle("Distance (cm)");
	h6->SetLineColor(kRed);

	TH1D * h7 = new TH1D("h7","Compton Vertex (Z), Parent Out",150,0,1036.8) ;
	TCanvas * c7 = new TCanvas("c7","c7") ;	
  t->Draw("_Z>>h7","_PDG == 3 && _distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)");
	t->Draw("_Z","_PDG == 3 && _E > 0.005 &&_distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)","same");
	h7->GetXaxis()->SetTitle("Distance (cm)");
	h7->SetLineColor(kRed);



	//Vertices of PP, Parent Out
	TH1D * h8 = new TH1D("h8","PP Vertex (X), Parent Out",30,0,256.35) ;
	TCanvas * c8 = new TCanvas("c8","c8") ;	
    t->Draw("_X>>h8","_PDG == 4 && _distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)");
    t->Draw("_X","_PDG == 4 && _E > 0.05 &&_distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)","same");
	h8->GetXaxis()->SetTitle("Distance (cm)");
	h8->SetLineColor(kRed);

	TH1D * h9 = new TH1D("h9","PP Vertex (Y), Parent Out",30,-116.5,116.5) ;
	TCanvas * c9 = new TCanvas("c9","c9") ;	
    t->Draw("_Y>>h9","_PDG == 4 && _distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)");
    t->Draw("_Y","_PDG == 4 && _E > 0.05 &&_distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)","same");
	h9->GetXaxis()->SetTitle("Distance (cm)");
	h9->SetLineColor(kRed);

	TH1D * h10 = new TH1D("h10","PP Vertex (Z), Parent Out",30,0,1036.8) ;
	TCanvas * c10 = new TCanvas("c10","c10") ;	
  t->Draw("_Z>>h10","_PDG == 4 && _distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)");
	t->Draw("_Z","_PDG == 4 && _E > 0.05 &&_distBackAlongTraj > 0 && (_parentX < 0 || _parentX > 256.35 || _parentY < -116.5 || _parentY > 116.5 || _parentZ < 0 || _parentZ > 1036.8)","same");
	h10->GetXaxis()->SetTitle("Distance (cm)");
	h10->SetLineColor(kRed);

	TLegend *l = new TLegend(0.6,0.6,0.7,0.75) ;
	l->AddEntry(h10,"Z, No Energy Cut","fl");   // h1 and h2 are histogram pointers
	l->Draw();

} //end testScript.C
