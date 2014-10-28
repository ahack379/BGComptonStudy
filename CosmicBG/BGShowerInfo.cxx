#ifndef BGSHOWERINFO_CXX
#define BGSHOWERINFO_CXX

#include "BGShowerInfo.h"

namespace larlite {

  bool BGShowerInfo::initialize() {

    _ana_tree=0;
    PrepareTTree();

   return true;
  }
  
  bool BGShowerInfo::analyze(storage_manager* storage) {

    auto my_mctruth = storage->get_data<event_mctruth>("generator") ;

    if(!my_mctruth) {
      std::cout<<"Hey I did not find mctruth made by generator!"<<std::endl;
      return true;
    }

    for( auto const & mct : * my_mctruth){
      std::cout<<"This is a test "<<std::endl;
    }    
/*        for(auto const& mcp : mct.GetParticles()) {
    

//			_energyGamma = mcp.Trajectory().at(0).E() ;

               if (mcp.PdgCode() ==11){
                    _elecOrGamma = 0;
                    _energyGamma = mcp.Trajectory().at(0).E() ;
    
                 }   
                else if(mcp.PdgCode() ==22){
                    _elecOrGamma = 1;   
                    _energyElec = mcp.Trajectory().at(0).E() ;
                }   

      //Also look at case where gamma has muon parent
 				else if(mcp.PdgCode() == 13){
                    _muonParent = 1 ; 
                    _energyMuon = mcp.Trajectory().at(0).E() ;

                }    

            }*/



    //if(_ana_tree)
    //      _ana_tree->Fill();

    //}

    return true;
  }

void BGShowerInfo::PrepareTTree() {

    if(!_ana_tree) {
      _ana_tree = new TTree("ana_tree","");

      _ana_tree->Branch("_elecOrGamma",&_elecOrGamma,"elecOrGamma/D") ;
   
      _ana_tree->Branch("_energyGamma",&_energyGamma,"energyGamma/D") ;
      _ana_tree->Branch("_energyElec",&_energyElec,"energyElec/D") ;

      _ana_tree->Branch("_energyMuon",&_energyMuon,"energyMuon/D") ;


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
