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
	_count2 = 0;

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


 //   for( auto const & mct : * my_mctruth){

//		std::cout<<"Size of GetPArticles: "<<mct.GetParticles().size() <<std::endl; 
        
   //     for(auto const& mcp : mct.GetParticles()) {
    
		for(auto const & mcp : * my_mcpart){

//			std::cout<<"Daughter thing: "<<mcp.Daughters().count(10)<<std::endl;

               if (mcp.PdgCode() == 11){


				    geoalgo::DistToBoxWall showerObject ;
					TrajectoryInVolume inVol ;
					
					inVol.SetVolume(0,256,-116,116,0,1037) ;
					_elecVtx = { mcp.Trajectory().at(0).X(), mcp.Trajectory().at(0).Y(), mcp.Trajectory().at(0).Z() };

				if ( inVol.PointInVolume(_elecVtx)&& mcp.Process()=="compt" ){

                    _elecOrGamma = 0;
                    _energyElec = mcp.Trajectory().at(0).E() ;

                    _elecMom = { mcp.Trajectory().at(0).Px(), mcp.Trajectory().at(0).Py(), mcp.Trajectory().at(0).Pz()};

					_inVolElecX = _elecVtx[0];
					_inVolElecY = _elecVtx[1];
					_inVolElecZ = _elecVtx[2];

					_count1 ++ ;
   	                _dist_ToWall        = showerObject.DistanceToWall(_elecVtx) ;
	                _dist_AlongTraj     = showerObject.DistanceToWall(_elecVtx,_gammaMom,1);
    	            _dist_BackAlongTraj = showerObject.DistanceToWall(_elecVtx,_gammaMom,0);
    			if(_ana_tree)
          			_ana_tree->Fill();
					}
    

                 }   

                else if(mcp.PdgCode() == 22){
					
                    _elecOrGamma = 1;   
                    _energyGamma = mcp.Trajectory().at(0).E() ;

					_gammaVtx = { mcp.Trajectory().at(0).X(), mcp.Trajectory().at(0).Y(), mcp.Trajectory().at(0).Z() };
                    _gammaMom = { mcp.Trajectory().at(0).Px(), mcp.Trajectory().at(0).Py(), mcp.Trajectory().at(0).Pz()};

                }   


 				else if(mcp.PdgCode() == 13){
                    _muonParent = 1 ; 
                    _energyMuon = mcp.Trajectory().at(0).E() ;

                }    

	
    }

    return true;
  }

void BGShowerInfo::PrepareTTree() {

    if(!_ana_tree) {
      _ana_tree = new TTree("ana_tree","");

      _ana_tree->Branch("_elecOrGamma",&_elecOrGamma,"elecOrGamma/D") ;
   
      _ana_tree->Branch("_energyGamma",&_energyGamma,"energyGamma/D") ;
      _ana_tree->Branch("_energyElec",&_energyElec,"energyElec/D") ;

	  _ana_tree->Branch("_inVolElecX",&_inVolElecX,"inVolElecX/D") ;
	  _ana_tree->Branch("_inVolElecY",&_inVolElecY,"inVolElecY/D") ;
	  _ana_tree->Branch("_inVolElecZ",&_inVolElecZ,"inVolElecZ/D") ;

      _ana_tree->Branch("_energyMuon",&_energyMuon,"energyMuon/D") ;

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
	
	std::cout<<" gammas in detector, gammas from comptons: "<<_count0<<", "<<_count1<<std::endl;

    return true;
  }

}
#endif
