#ifndef COSMICSINFO_CXX
#define COSMICSINFO_CXX

#include "CosmicsInfo.h"
#include <set>

namespace larlite {

  bool CosmicsInfo::initialize() {

    _ana_tree=0;

    PrepareTTree();

    detHalfWidth = ::larutil::Geometry::GetME()->DetHalfWidth();
    detHalfHeight = ::larutil::Geometry::GetME()->DetHalfHeight();
    detLength = ::larutil::Geometry::GetME()->DetLength();

   return true;
  }
  
  bool CosmicsInfo::analyze(storage_manager* storage) {

    auto my_mcpart = storage->get_data<event_mcpart>("largeant") ;

    if(!my_mcpart) {
      std::cout<<"No mcpart made by largeant"<<std::endl;
      return true;
    }

    //get particles with mcgetter
    _MCgetter.Reset(my_mcpart);
    std::vector<TreeNode> result = _MCgetter.getTreeNodelist().at(0);
    std::vector<TreeNode> muon = _MCgetter.getTreeNodelist().at(1);
    std::vector<TreeNode> antimuon = _MCgetter.getTreeNodelist().at(2);


    // first get all mu+/mu-
    _allMuonTracksInTPC.clear();
    for (size_t y=0; y < muon.size(); y++){
      mcpart mu = my_mcpart->at(_MCgetter.searchParticleMap(muon.at(y).getNodeIndex()));
      addMuonTrack(&(my_mcpart->at(y)));
    }//for all mu-
    for (size_t k=0; k < antimuon.size(); k++){
      mcpart mu = my_mcpart->at(_MCgetter.searchParticleMap(antimuon.at(k).getNodeIndex()));
      addMuonTrack(&(my_mcpart->at(k)));
    }//for all mu+

    
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
      
      // do muon-track cuts
      getNearestMuonCutParams( &vtx, &mom);
      
      
      _ana_tree->Fill();
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

	  _ana_tree->Branch("_minMuDist",&_minMuDist,"minMuDist/D");
	  _ana_tree->Branch("_minMuIP",&_minMuIP,"minMuIP/D");
	  _ana_tree->Branch("_distToIP",&_distToIP,"distToIP/D");

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

//	  _ana_tree->Branch("MuonTraj",&MuonTraj) ;
	
	}

 }

  void CosmicsInfo::addMuonTrack(mcpart *part){

    mctrajectory traj = part->Trajectory();
    std::vector<std::vector<double> > thistrack;

    for (size_t i=0; i < traj.size(); i++){

      // check if point in TPC
      if ( (traj.at(i).X() > 0) and (traj.at(i).X() < 2*detHalfWidth) and
	   (traj.at(i).Y() > -detHalfHeight) and (traj.at(i).Y() < detHalfHeight) and
	   (traj.at(i).Z() > 0) and (traj.at(i).Z() < detLength) )
	   thistrack.push_back( {traj.at(i).X(), traj.at(i).Y(), traj.at(i).Z()} );
				     
    }//for all points
    
    if (thistrack.size() > 0)
      _allMuonTracksInTPC.push_back(thistrack);

    return;
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

   _minMuDist = -1;
   _minMuIP = -1;
   _distToIP = -1;

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



  void CosmicsInfo::getNearestMuonCutParams( std::vector<double> *shrStart,
					     std::vector<double> *shrDir){

    // use shower's start point + direction and list of muon trajectory to find:
    // 1) PoCA and PoCA distance to shower start point w/ PoCA GeoAlgo class
    // 2) MuonDist for cylinder cut
    
    // initialize values for PoCA cut
    double minIP = 10000;
    std::vector<double> IP_onShr = {-1000, -1000, -1000};
    std::vector<double> c1 = {-1000,-1000,-1000};
    std::vector<double> c2 = {-1000,-1000,-1000};
    std::vector<double> PoCAPointMU = {-1000,-1000,-1000};
    std::vector<double> PoCAPointE = {-1000,-1000,-1000};
    // initialize values for muon proximity
    double minDist = 10000;
    
    // use shower start & momentum to define a segment
    // which starts 3 meters before start point
    // and ends 10 cm after start point
    // aligned with shr momentum. Use for Impact Parameter
    std::vector<double> shrOrigin = { shrStart->at(0)-300*shrDir->at(0),
				      shrStart->at(1)-300*shrDir->at(1),
				      shrStart->at(2)-300*shrDir->at(2) };
    
    std::vector<double> shrEnd = { shrStart->at(0)+10*shrDir->at(0),
				   shrStart->at(1)+10*shrDir->at(1),
				   shrStart->at(2)+10*shrDir->at(2) };
    
    // loop over all muon tracks and calculate value per-muon
  for (size_t u=0; u < _allMuonTracksInTPC.size(); u++){
    
    if (_allMuonTracksInTPC.at(u).size() > 1){
      // distance to muon track
      double tmpDist = _pointDist.DistanceToTrack(shrStart, &(_allMuonTracksInTPC.at(u)));
      // Impact parameter
      double tmpIP = _PoCA.ClosestApproachToTrajectory(&(_allMuonTracksInTPC.at(u)), &shrOrigin, &shrEnd, c1, c2);
      
      if (tmpDist < minDist) { minDist = tmpDist; }
      if (tmpIP < minIP) { 
	minIP = tmpIP; 
	IP_onShr = c2;
      }// if best muon so far
      
    }// if the muon's track is > 1
    
  }// for all muon tracks
  
  _distToIP = pow( (IP_onShr.at(0)-shrStart->at(0))*(IP_onShr.at(0)-shrStart->at(0)) +
		   (IP_onShr.at(1)-shrStart->at(1))*(IP_onShr.at(1)-shrStart->at(1)) +
		   (IP_onShr.at(2)-shrStart->at(2))*(IP_onShr.at(2)-shrStart->at(2)), 0.5 );
  
  
  _minMuDist = sqrt(minDist);
  _minMuIP = sqrt(minIP);
  
  return;
  }
  
  
}
#endif
