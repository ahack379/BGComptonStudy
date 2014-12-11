#ifndef BGSHOWERINFO_CXX
#define BGSHOWERINFO_CXX

#include "BGShowerInfo.h"
#include <set>

namespace larlite {

  bool BGShowerInfo::initialize() {

    _ana_tree=0;
	_gamma_tree = 0;
	_gamma_tree2 = 0;
	_pp_tree = 0;
	_muon_tree = 0;
    PrepareTTree();

	_count0 = 0;
	_count1 = 0;

energy = {1.000E-03, 1.500E-03, 2.000E-03, 3.000E-03, 3.203E-03, 3.203E-03, 4.000E-03, 5.000E-03, 6.000E-03, 8.000E-03, 1.000E-02, 1.500E-02, 2.000E-02,
	 3.000E-02, 4.000E-02, 5.000E-02, 6.000E-02, 8.000E-02, 1.000E-01, 1.500E-01, 2.000E-01, 3.000E-01, 4.000E-01, 5.000E-01, 6.000E-01, 8.000E-01,
	1.00E+00, 1.02E+00, 1.25E+00, 1.50E+00, 2.00E+00, 2.04E+00, 3.00E+00, 4.00E+00, 5.00E+00, 6.00E+00, 7.00E+00, 8.00E+00, 9.00E+00,
	1.00E+01, 1.10E+01, 1.20E+01, 1.30E+01, 1.40E+01, 1.50E+01, 1.60E+01, 1.80E+01, 2.00E+01, 2.20E+01, 2.40E+01, 2.60E+01, 2.80E+01,
	3.00E+01, 4.00E+01, 5.00E+01, 6.00E+01, 8.00E+01, 1.00E+02, 1.50E+02, 2.00E+02, 3.00E+02, 4.00E+02, 5.00E+02, 6.00E+02, 8.00E+02,
	1.00E+03, 1.50E+03, 2.00E+03, 3.00E+03, 4.00E+03, 5.00E+03, 6.00E+03, 8.00E+03, 1.00E+04 }; 
	
pp = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000827184, 0.00483919, 0.0217575, 0.0236998, 0.0741492, 0.134803, 0.195148,
	0.252186, 0.304672, 0.352269, 0.395028, 0.433578, 0.467831, 0.498749, 0.526898, 0.552102, 0.57489, 0.59584, 0.632771, 0.663934, 0.690429, 0.712962, 0.732882,
	0.750193, 0.765527, 0.820963, 0.855809, 0.87924, 0.909057, 0.926784, 0.950663, 0.962507, 0.974383, 0.98046, 0.98407, 0.986292, 0.989478, 0.991496, 0.994241,
	0.995459, 0.997023, 0.997629, 0.998227, 0.998205, 0.998823, 0.999004 } ;



   return true;
  }
  
  bool BGShowerInfo::analyze(storage_manager* storage) {

//    auto my_mctruth = storage->get_data<event_mctruth>("generator") ;
    auto my_mcpart = storage->get_data<event_mcpart>("largeant") ;

	if(!my_mcpart) {
		  std::cout<<"No mcpart made by largeant"<<std::endl;
		  return true;
		}

	double detHalfWidth = ::larutil::Geometry::GetME()->DetHalfWidth();
	double detHalfHeight = ::larutil::Geometry::GetME()->DetHalfHeight();
	double detLength = ::larutil::Geometry::GetME()->DetLength();
	
	/// Set TPC Boundaries
	_TpcBox = geoalgo::AABox(0, -detHalfHeight, 0,
				 2*detHalfWidth, detHalfHeight, detLength);

	//   Loop over all mcparticles and look at various processes 
	//1) Compton electron with vtx in TPC
	//2) Pair produced gamma in TPC
	//3) Check that Pair production rates should be as low as they
	//	 are given the vectors defined above.


		//Check if muon track makes it through detector
		//First add traj points to vector
	/*	if(mcp.PdgCode() ==13  || mcp.PdgCode() == -13){
			_count1++ ;
			std::vector<double> _muonTrajPt ; 

			int j=0;
			for(; j<mcp.Trajectory().size(); j++ ){
				_muonTrajPt	= {mcp.Trajectory().at(j).X(), mcp.Trajectory().at(j).Y(), mcp.Trajectory().at(j).Z() };

					if(inVol.PointInVolume(_muonTrajPt) ){
		//				std::cout<<"Tr crosses : "<<_muonTrajPt[0]<<" " <<_muonTrajPt[1]<<" "<<_muonTrajPt[2]<<std::endl;
						_timeMuonEnter = mcp.Trajectory().at(j).T() ;
						
						if(_muon_tree)
							_muon_tree->Fill();
						break ;

					  }
			 }	

		}*/
	for(auto const & mcp : * my_mcpart){

	   if (mcp.PdgCode() == 11 && mcp.Process() == "compt"){

			_elecVtx = { mcp.Trajectory().at(0).X(), mcp.Trajectory().at(0).Y(), mcp.Trajectory().at(0).Z() };

			if ( _TpcBox.Contain(geoalgo::Vector(_elecVtx)) ){
				int	MotherOfComp =  mcp.Mother() ;

				for(auto const& mcp2 : * my_mcpart){

				    if(mcp2.TrackId() == MotherOfComp ){ 

						std::vector<double> temp_gamma_vtx(3,0) ; 
						temp_gamma_vtx = { mcp2.Trajectory().at(0).X(), mcp2.Trajectory().at(0).Y(),  mcp2.Trajectory().at(0).Z() };

						if(!_TpcBox.Contain(geoalgo::Vector(temp_gamma_vtx))) {
							_gammaCompX = temp_gamma_vtx[0]; 
							_gammaCompY = temp_gamma_vtx[1]; 
							_gammaCompZ = temp_gamma_vtx[2]; 
	
	
							_energyElec = mcp.Trajectory().at(0).E() ;
							_elecMom = { mcp.Trajectory().at(0).Px(), mcp.Trajectory().at(0).Py(), mcp.Trajectory().at(0).Pz()};
							//_timeElecEnter = mcp.Trajectory().at(0).T() ;
							
							_pxElec = _elecMom[0]; 
							_pyElec = _elecMom[1]; 
							_pzElec = _elecMom[2]; 
			
							_inVolElecX = _elecVtx[0];
							_inVolElecY = _elecVtx[1];
							_inVolElecZ = _elecVtx[2];
		
							geoalgo::Point_t eVtx(_elecVtx); 
							geoalgo::HalfLine_t sDir(_elecVtx,_elecMom);
							_dist_ToWall = sqrt(_geoAlgo.SqDist(eVtx,_TpcBox));
							_dist_AlongTraj = eVtx.Dist(_geoAlgo.Intersection(_TpcBox,sDir)[0]);
							_dist_BackAlongTraj = eVtx.Dist(_geoAlgo.Intersection(_TpcBox,sDir,true)[0]);

							if(_ana_tree)	
								_ana_tree->Fill();
							break;
			
						}
	

					//	std::cout<<"trackID and mother: "<<mcp2.PdgCode()<<" " <<mcp2.Mother() <<" " <<mcp2.TrackId() <<"  "<<MotherOfComp<<std::endl ; }
					
		/*				for(auto const & mcp3 : *my_mcpart){

						//This block was a quick check on the mothers of gamma--leave it commented out
						//unless totally necessary--it contains misleading variables.
							if(mcp3.TrackId() == MotherOfGamma ){
								_parents = mcp3.PdgCode() ;
								_process = mcp3.Process() ;
								_energyGammaTotal = mcp2.Trajectory().at(0).E() ;
								_energyMother = mcp3.Trajectory().at(0).E() ;
								if(_gamma_tree)
									_gamma_tree->Fill() ;
								if(mcp3.Trajectory().at(0).E() > 0.05){
									_parentsCut = mcp3.PdgCode();

								if(_muon_tree)
									_muon_tree->Fill();
								  }
								break; */
							
						//This block checks the case where mother of the mother gamma is a muon.  Too specific.
						/*	if(mcp3.TrackId() == MotherOfGamma &&( mcp3.PdgCode() == 13 || mcp3.PdgCode()==-13 ) && mcp3.Trajectory().at(0).E() > 0.05){
								std::cout<<"Muon? "<<std::endl;

								std::vector<double> _muonTrajPt ; 
	
								int j=0;
								for(; j<mcp3.Trajectory().size(); j++ ){
									_muonTrajPt	= {mcp3.Trajectory().at(j).X(), mcp3.Trajectory().at(j).Y(), mcp3.Trajectory().at(j).Z() };

									if(inVol.PointInVolume(_muonTrajPt) ){
									_timeMuonEnter = mcp3.Trajectory().at(j).T() ;
				//					std::cout<<"muon and electron enter at : "<< _timeMuonEnter<<" "<<_timeElecEnter<<std::endl ;
									std::cout<<"muon IDs : "<<MotherOfGamma<<std::endl ; 
									std::cout<<"Origin : "<<mcp2.Trajectory().at(0).X() << " "<<mcp2.Trajectory().at(0).Y() << " "<<mcp2.Trajectory().at(0).Z()<<std::endl;
									//std::cout<<"gamma process and ID: "<<mcp2.Process()<<" " <<mcp2.TrackId() <<std::endl;
									_timeDiff = - _timeMuonEnter + _timeElecEnter ;
									std::cout<<"time difference: "<<_timeDiff<<std::endl;
						
										if(_muon_tree)
											_muon_tree->Fill();
									break ;
									  }
								  }
							
							*/ //	}	
					
						}

					} 


				}
		 	}   

	else if(mcp.PdgCode() == -11 && mcp.Process()=="conv" ) {

			  _positronVtx = { mcp.Trajectory().at(0).X(), mcp.Trajectory().at(0).Y(), mcp.Trajectory().at(0).Z() };

			  if ( _TpcBox.Contain(geoalgo::Vector(_positronVtx)) ){
				int	MotherID =  mcp.Mother() ;

				for(auto const& mcp2 : * my_mcpart){
				    if(mcp2.TrackId() ==MotherID ){

						std::vector<double> temp_gamma_vtx(3,0) ; 
						temp_gamma_vtx = { mcp2.Trajectory().at(0).X(), mcp2.Trajectory().at(0).Y(),  mcp2.Trajectory().at(0).Z() };

						if(!_TpcBox.Contain(geoalgo::Vector(temp_gamma_vtx))) {
							_energyGammaBegin = mcp2.Trajectory().at(0).E() ;
							_energyGammaEnd = mcp2.Trajectory().at(mcp2.Trajectory().size()-2).E() ;
							_inVolPosX = _positronVtx[0];
							_inVolPosY = _positronVtx[1];
							_inVolPosZ = _positronVtx[2];

							_pxPP = mcp2.Trajectory().at(mcp2.Trajectory().size()-2).Px() ;
							_pyPP = mcp2.Trajectory().at(mcp2.Trajectory().size()-2).Py() ;
							_pzPP = mcp2.Trajectory().at(mcp2.Trajectory().size()-2).Pz() ;
		
							_ppMom = {_pxPP, _pyPP, _pzPP };
							geoalgo::Point_t ppStart(_positronVtx);
							geoalgo::HalfLine_t ppDir(_positronVtx,_ppMom);
							_dist_BackAlongTrajPP = ppStart.Dist(_geoAlgo.Intersection(_TpcBox,ppDir,true)[0]);

							if(_pp_tree)
								_pp_tree->Fill();
							break;
					 }
				  }
				}

				  }   

			 	 }
			
	/*		else if ( mcp.PdgCode()==22){
			////////////////////////PairProd check
				_rand = drand48() ;

				std::vector<double> temp_gamma_vtx(3,0) ; 
				temp_gamma_vtx = { mcp.Trajectory().at(0).X(), mcp.Trajectory().at(0).Y(), temp_gamma_vtx[2] = mcp.Trajectory().at(0).Z() };

				if(inVol.PointInVolume(temp_gamma_vtx)) {
					_energyGammaTotal = mcp.Trajectory().at(0).E() ;
					_inVolGammaX = temp_gamma_vtx[0]; 
					_inVolGammaY = temp_gamma_vtx[1]; 
					_inVolGammaZ = temp_gamma_vtx[2]; 


					int j=0;
					while(j < energy.size()){
						if(_energyGammaTotal < energy[j]/1000)
							break;
						j++;
						}
					_prob = pp[j] ;

					if( _rand < _prob ){ 
						_checkGammaPP = _energyGammaTotal ;
						_gamma_tree2->Fill() ;
						}
			//////////////////////////PairProd check end


					if(_gamma_tree)
						_gamma_tree->Fill() ;
					}

				} */
		}				

    return true;
  }

void BGShowerInfo::PrepareTTree() {

    if(!_ana_tree) {
      _ana_tree = new TTree("ana_tree","");

	  _ana_tree->Branch("_timeElecEnter",&_timeElecEnter,"timeElecEnter/D");
	  _ana_tree->Branch("_timing",&_timing,"timing/D"); 
   
      _ana_tree->Branch("_energyElec",&_energyElec,"energyElec/D") ;
	  _ana_tree->Branch("_pxElec",&_pxElec,"pxElec/D");
	  _ana_tree->Branch("_pyElec",&_pyElec,"pyElec/D");
	  _ana_tree->Branch("_pzElec",&_pzElec,"pzElec/D");
	  _ana_tree->Branch("_inVolElecX",&_inVolElecX,"inVolElecX/D") ;
	  _ana_tree->Branch("_inVolElecY",&_inVolElecY,"inVolElecY/D") ;
	  _ana_tree->Branch("_inVolElecZ",&_inVolElecZ,"inVolElecZ/D") ;
	  //_ana_tree->Branch("_pdgCode",&_pdgCode,"pdgCode/D") ;


	  _ana_tree->Branch("_dist_ToWall",&_dist_ToWall,"dist_ToWall/D") ;
	  _ana_tree->Branch("_dist_AlongTraj",&_dist_AlongTraj,"dist_AlongTraj/D") ;
	  _ana_tree->Branch("_dist_BackAlongTraj",&_dist_BackAlongTraj,"dist_BackAlongTraj/D") ;

    }
	if(!_pp_tree){

      _pp_tree = new TTree("pp_tree","");	

      _pp_tree->Branch("_energyGammaBegin",&_energyGammaBegin,"energyGammaBegin/D") ;
      _pp_tree->Branch("_energyGammaEnd",&_energyGammaEnd,"energyGammaEnd/D") ;
      _pp_tree->Branch("_energyElec",&_energyElec,"energyElec/D") ;
	  _pp_tree->Branch("_pxPP",&_pxPP,"pxPP/D");
	  _pp_tree->Branch("_pyPP",&_pyPP,"pyPP/D");
	  _pp_tree->Branch("_pzPP",&_pzPP,"pzPP/D");
	  _pp_tree->Branch("_inVolPosX",&_inVolPosX,"inVolPosX/D") ;
	  _pp_tree->Branch("_inVolPosY",&_inVolPosY,"inVolPosY/D") ;
	  _pp_tree->Branch("_inVolPosZ",&_inVolPosZ,"inVolPosZ/D") ;
	}	

	if(!_gamma_tree){

	  _gamma_tree = new TTree("gamma_tree","");

	  _gamma_tree->Branch("_process","std::string",&_process);
	_gamma_tree->Branch("_energyMother",&_energyMother,"energyMother/D");
	_gamma_tree->Branch("_parents",&_parents,"parents/D");
	  _gamma_tree->Branch("_energyGammaTotal",&_energyGammaTotal,"energyGammaTotal/D") ;
	  _gamma_tree->Branch("_inVolGammaX",&_inVolGammaX,"inVolGammaX/D") ;
	  _gamma_tree->Branch("_inVolGammaY",&_inVolGammaY,"inVolGammaY/D") ;
	  _gamma_tree->Branch("_inVolGammaZ",&_inVolGammaZ,"inVolGammaZ/D") ;



		}
	
	if(!_gamma_tree2){
	
	  _gamma_tree2 = new TTree("gamma_tree2","") ;

	  _gamma_tree2->Branch("_rand",&_rand,"rand/D") ;
	  _gamma_tree2->Branch("_checkGammaPP",&_checkGammaPP,"checkGammaPP/D") ;
	
	}

	if(!_muon_tree){
	
	_muon_tree = new TTree("muon_tree","");

	_muon_tree->Branch("_parentsCut",&_parentsCut,"parentsCut/D");
	_muon_tree->Branch("_timeMuonEnter",&_timeMuonEnter,"timeMuonEnter/D") ;
	_muon_tree->Branch("_timeElecEnter",&_timeElecEnter,"timeElecEnter/D") ;
	_muon_tree->Branch("_timeDiff",&_timeDiff,"timeDiff/D");
	
	
	}

 }

void BGShowerInfo::Clear() {

  }

void BGShowerInfo::Reset(){

 }

bool BGShowerInfo::finalize() {

    if(_fout) {

      _fout->cd();
      if(_ana_tree && _gamma_tree && _pp_tree && _gamma_tree2 && _muon_tree){
        _ana_tree->Write();
		_gamma_tree->Write();
		_gamma_tree2->Write();
		_pp_tree->Write();
		_muon_tree->Write() ;
		}
      }
     else
       print(larlite::msg::kERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");

	 std::cout<<"Muons ! : "<<_count1<<std::endl;

      delete _ana_tree;
	  delete _gamma_tree;
	  delete _gamma_tree2;
	  delete _pp_tree;
	  delete _muon_tree ;
	
    return true;
  }

}
#endif
