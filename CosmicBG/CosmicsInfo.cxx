#ifndef COSMICSINFO_CXX
#define COSMICSINFO_CXX

#include "CosmicsInfo.h"
#include <set>

namespace larlite {

  bool CosmicsInfo::initialize() {

    _ana_tree=0;

    PrepareTTree();

   return true;
  }
  
  bool CosmicsInfo::analyze(storage_manager* storage) {

    auto my_mcpart = storage->get_data<event_mcpart>("largeant") ;

	if(!my_mcpart) {
		  std::cout<<"No mcpart made by largeant"<<std::endl;
		  return true;
		}

	_inVol.SetVolume(0,256.35,-116.5,116.5,0,1036.8) ;

	for(auto const & mcp : * my_mcpart){

		Reset();

		//Create mcshower-like things
		if( mcp.PdgCode() == 11 && mcp.Process() == "compt")
			_PDG = 3 ;
		
		if( mcp.PdgCode() == -11 && mcp.Process() == "conv"){
			_PDG = 4 ;
			int posMother = mcp.Mother();
			for(auto const & mcp2 : * my_mcpart){
				
				if(posMother == mcp2.TrackId()){
					_E = mcp2.Trajectory().at(mcp2.Trajectory().size()-2).E() ;
					_Px = mcp2.Trajectory().at(mcp2.Trajectory().size()-2).Px() ;
					_Py = mcp2.Trajectory().at(mcp2.Trajectory().size()-2).Py() ;
					_Pz = mcp2.Trajectory().at(mcp2.Trajectory().size()-2).Pz() ;
				  }
			   }
			}//if positron, conv

		//Only save particles which are gammas, pi0s, compton scatters and pair productions
//		if( mcp.PdgCode() == 111 || mcp.PdgCode() == 22 || _PDG == 3 || _PDG == 4 ){
			
			_run = my_mcpart->run() ;
			_subrun = my_mcpart->subrun();
			_event = my_mcpart->event_id(); 
	
			_PDG = mcp.PdgCode() ;
			_trackID = mcp.TrackId() ;
			_process= mcp.Process() ;
	
			_X = mcp.Trajectory().at(0).X() ;
			_Y = mcp.Trajectory().at(0).Y() ;
			_Z = mcp.Trajectory().at(0).Z() ;
			_T = mcp.Trajectory().at(0).T() ;
	
			_Px = mcp.Trajectory().at(0).Px() ;
			_Py = mcp.Trajectory().at(0).Py() ;
			_Pz = mcp.Trajectory().at(0).Pz() ;
			_E = mcp.Trajectory().at(0).E() ;
	
			std::vector<double> vtx = { _X, _Y, _Z } ;
			std::vector<double> mom = { _Px, _Py, _Pz } ;

			if(_inVol.PointInVolume(vtx ))
				_inActiveVolume = 1 ; 
	
			_distAlongTraj     = _showerObject.DistanceToWall(vtx,mom,1);
			_distBackAlongTraj = _showerObject.DistanceToWall(vtx,mom,0);
	
			//Get Parent info as well
			int motherID = mcp.Mother(); 
			for(auto const & mcp3 : * my_mcpart){
			
				if( motherID == mcp3.TrackId()){
					_parentPDG = mcp3.PdgCode() ;
	
					_parentX = mcp3.Trajectory().at(0).X();
					_parentY = mcp3.Trajectory().at(0).Y();
					_parentZ = mcp3.Trajectory().at(0).Z();
					_parentT = mcp3.Trajectory().at(0).T();
	
					_parentPx = mcp3.Trajectory().at(0).Px();
					_parentPy = mcp3.Trajectory().at(0).Py();
					_parentPz = mcp3.Trajectory().at(0).Pz();
					_parentE = mcp3.Trajectory().at(0).E();
	
					std::vector<double> pVtx = { _parentX, _parentY, _parentZ } ; 
					std::vector<double> pMom = { _parentPx, _parentPy, _parentPz } ;
				
					if(_inVol.PointInVolume(pVtx ))
						_inActiveVolume = 0;  

				 	}
				}

			//Adding muons trajectory information for DavidC's SaLSA cuts,etc
			if(mcp.PdgCode() == 13 || mcp.PdgCode() == -13) {
				std::cout<<"MUON! "<<std::endl;
			    std::vector< std::vector< std::vector<double> > > muonTracks;
				std::vector< std::vector<double> > muonTraj ;
			    muonTracks.clear();
				muonTraj.resize(0);

			    for (size_t h=0; h < mcp.Trajectory().size(); h++){

					std::vector<double> muonVtx = { mcp.Trajectory().at(h).X(),
													mcp.Trajectory().at(h).Y(), 
													mcp.Trajectory().at(h).Z() };

					if(_inVol.PointInVolume(muonVtx)){   
					   muonTraj.push_back(muonVtx) ;

						std::cout<<"\n1Muon things : "<<muonVtx[0]<<std::endl;
						std::cout<<"2Muon things : "<<muonTraj.at(0).at(0) <<std::endl; 
					 }
				 }
				MuonTraj = muonTraj ;

				if(muonTraj.size() >=1 )
				    muonTracks.push_back(muonTraj);

				}//if muon/ antimuon
	
				_ana_tree->Fill();
		//	} //End if statement conditions for which pdgs to save
		}				
	

    return true;
  }

void CosmicsInfo::PrepareTTree() {

    if(!_ana_tree) {
      _ana_tree = new TTree("ana_tree","");

	  _ana_tree->Branch("_run",&_run,"run/I");
	  _ana_tree->Branch("_subrun",&_subrun,"subrun/I");
	  _ana_tree->Branch("_event",&_event,"event/I");

	  _ana_tree->Branch("_process","std::string",&_process);
	  _ana_tree->Branch("_PDG",&_PDG,"PDG/I");
	  _ana_tree->Branch("_trackID",&_trackID,"trackID/I");

	  _ana_tree->Branch("_X",&_X,"X/D");
	  _ana_tree->Branch("_Y",&_Y,"Y/D");
	  _ana_tree->Branch("_Z",&_Z,"Z/D");
	  _ana_tree->Branch("_T",&_T,"T/D");

	  _ana_tree->Branch("_Px",&_Px,"Px/D");
	  _ana_tree->Branch("_Py",&_Py,"Py/D");
	  _ana_tree->Branch("_Pz",&_Pz,"Pz/D");
	  _ana_tree->Branch("_E",&_E,"E/D");

	  _ana_tree->Branch("_inActiveVolume",&_inActiveVolume,"inActiveVolume/D");
	  
	  _ana_tree->Branch("_distAlongTraj",&_distAlongTraj,"distAlongTraj/D") ;
	  _ana_tree->Branch("_distBackAlongTraj",&_distBackAlongTraj,"distBackAlongTraj/D") ;

	  ////PARENT INFO
	  //
	  _ana_tree->Branch("_parentPDG",&_parentPDG,"parentPDG/I");
	  _ana_tree->Branch("_parentX",&_parentX,"parentX/D");
	  _ana_tree->Branch("_parentY",&_parentY,"parentY/D");
	  _ana_tree->Branch("_parentZ",&_parentZ,"parentZ/D");
	  _ana_tree->Branch("_parentT",&_parentT,"parentT/D");

	  _ana_tree->Branch("_parentPx",&_parentPx,"parentPx/D");
	  _ana_tree->Branch("_parentPy",&_parentPy,"parentPy/D");
	  _ana_tree->Branch("_parentPz",&_parentPz,"parentPz/D");
	  _ana_tree->Branch("_parentE",&_parentE,"parentE/D");

	  _ana_tree->Branch("_parentInActiveVolume",&_parentInActiveVolume,"parentInActiveVolume/D");

	  _ana_tree->Branch("MuonTraj",&MuonTraj) ;
	
	}

 }

void CosmicsInfo::Clear() {

  }

void CosmicsInfo::Reset(){


   _run     = -1;
   _subrun  = -1;
   _event 	= -1;

   _process = "NONE";
   _PDG  	= -1;
   _trackID = -1;

   _X  		= -9999999;
   _Y 		= -9999999;
   _Z 		= -9999999;
   _T 		= -9999999;

   _Px 		= -9999999;
   _Py 		= -9999999;
   _Pz 		= -9999999;
   _E 		= -9999999;

   _distAlongTraj     = -9999999;
   _distBackAlongTraj = -9999999;

   //parent
   _parentPDG = -1;
   _parentX   = -9999999;
   _parentY   = -9999999;
   _parentZ   = -9999999;
   _parentT   = -9999999;

   _parentPx  = -9999999;
   _parentPy  = -9999999;
   _parentPz  = -9999999;
   _parentE   = -9999999;

   _inActiveVolume = -99 ;
   _parentInActiveVolume = -99 ;

 }

bool CosmicsInfo::finalize() {

    if(_fout) {

      _fout->cd();
      if(_ana_tree) 
        _ana_tree->Write();
		}
      
     else
       print(larlite::msg::kERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");

      delete _ana_tree;
	
    return true;
  }

}
#endif
