#ifndef BGSHOWERINFO_CXX
#define BGSHOWERINFO_CXX

#include "BGShowerInfo.h"
#include <set>

namespace larlite {

  bool BGShowerInfo::initialize() {

    _ana_tree=0;
    PrepareTTree();

	_count0 = 0;
	_count1 = 0;

   return true;
  }
  
  bool BGShowerInfo::analyze(storage_manager* storage) {

//    auto my_mctruth = storage->get_data<event_mctruth>("generator") ;
    auto my_mcpart = storage->get_data<event_mcpart>("largeant") ;

/*    if(!my_mctruth) {
      std::cout<<"Hey I did not find mctruth made by generator!"<<std::endl;
      return true;
    }
*/
	if(!my_mcpart) {
	  std::cout<<"No mcpart made by largeant"<<std::endl;
	  return true;
	}

	    geoalgo::DistToBoxWall showerObject ;
		TrajectoryInVolume inVol ;
		inVol.SetVolume(0,256,-116,116,0,1037) ;

		for(auto const & mcp : * my_mcpart){

               if (mcp.PdgCode() == 11 && mcp.Process() == "compt"){

					_elecVtx = { mcp.Trajectory().at(0).X(), mcp.Trajectory().at(0).Y(), mcp.Trajectory().at(0).Z() };

				if ( inVol.PointInVolume(_elecVtx) && mcp.Process()=="compt" ){

					_energyElec = mcp.Trajectory().at(0).E() ;
                    _elecMom = { mcp.Trajectory().at(0).Px(), mcp.Trajectory().at(0).Py(), mcp.Trajectory().at(0).Pz()};

					_inVolElecX = _elecVtx[0];
					_inVolElecY = _elecVtx[1];
					_inVolElecZ = _elecVtx[2];

   	                _dist_ToWall        = showerObject.DistanceToWall(_elecVtx) ;
	                _dist_AlongTraj     = showerObject.DistanceToWall(_elecVtx,_elecMom,1);
    	            _dist_BackAlongTraj = showerObject.DistanceToWall(_elecVtx,_elecMom,0);


					}
                 }   

				else if(mcp.PdgCode() == -11 && mcp.Process()=="conv" ) {

				  _positronVtx = { mcp.Trajectory().at(0).X(), mcp.Trajectory().at(0).Y(), mcp.Trajectory().at(0).Z() };

				  if ( inVol.PointInVolume(_positronVtx) && mcp.Process()=="conv" ){
					int	MotherID =  mcp.Mother() ;
                    _positronMom = { mcp.Trajectory().at(0).Px(), mcp.Trajectory().at(0).Py(), mcp.Trajectory().at(0).Pz()};
					for(auto const& mcp2 : * my_mcpart){
					  if(mcp2.TrackId() ==MotherID ){
                        _energyGammaBegin = mcp2.Trajectory().at(0).E() ;
                        _energyGammaEnd = mcp2.Trajectory().at(mcp2.Trajectory().size()-2).E() ;
							}
						}
    				if(_ana_tree)
          				_ana_tree->Fill();
                	 }   

    		 	  }


			}
		

    return true;
  }

void BGShowerInfo::PrepareTTree() {

    if(!_ana_tree) {
      _ana_tree = new TTree("ana_tree","");

   
      _ana_tree->Branch("_energyGammaBegin",&_energyGammaBegin,"energyGammaBegin/D") ;
      _ana_tree->Branch("_energyGammaEnd",&_energyGammaEnd,"energyGammaEnd/D") ;
      _ana_tree->Branch("_energyElec",&_energyElec,"energyElec/D") ;

	  _ana_tree->Branch("_inVolElecX",&_inVolElecX,"inVolElecX/D") ;
	  _ana_tree->Branch("_inVolElecY",&_inVolElecY,"inVolElecY/D") ;
	  _ana_tree->Branch("_inVolElecZ",&_inVolElecZ,"inVolElecZ/D") ;

	  _ana_tree->Branch("_dist_ToWall",&_dist_ToWall,"dist_ToWall/D") ;
	  _ana_tree->Branch("_dist_AlongTraj",&_dist_AlongTraj,"dist_AlongTraj/D") ;
	  _ana_tree->Branch("_dist_BackAlongTraj",&_dist_BackAlongTraj,"dist_BackAlongTraj/D") ;


    }
 }

void BGShowerInfo::Clear() {


  }

void BGShowerInfo::Reset(){

 }

bool BGShowerInfo::finalize() {

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
