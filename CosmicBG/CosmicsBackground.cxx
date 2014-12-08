#ifndef COSMICSBACKGROUND_CXX
#define COSMICSBACKGROUND_CXX

#include "CosmicsBackground.h"
#include <set>

namespace larlite {

  bool CosmicsBackground::initialize() {

    _ana_tree=0;

    PrepareTTree();

    detHalfWidth = ::larutil::Geometry::GetME()->DetHalfWidth();
    detHalfHeight = ::larutil::Geometry::GetME()->DetHalfHeight();
    detLength = ::larutil::Geometry::GetME()->DetLength();

    /// Set TPC Boundaries
    _TpcBox = geoalgo::AABox(0, -detHalfHeight, 0,
			     2*detHalfWidth, detHalfHeight, detLength);

    _hAllMuonTrackLen = new TH1D("hAllMuonTrackLen","Sum Length of All Muon Tracks in TPC",100,0,40);

   return true;
  }
  
  bool CosmicsBackground::analyze(storage_manager* storage) {

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

    std::cout << "Muons: " << muon.size() << std::endl;
    // first get all mu+/mu-
    _allMuonTracksInTPC.clear();
    _allMuonTracksIDs.clear();
    for (size_t y=0; y < muon.size(); y++){
      mcpart mu = my_mcpart->at(_MCgetter.searchParticleMap(muon.at(y).getNodeIndex()));
      addMuonTrack(mu);
    }//for all mu-
    for (size_t k=0; k < antimuon.size(); k++){
      mcpart antimu = my_mcpart->at(_MCgetter.searchParticleMap(antimuon.at(k).getNodeIndex()));
      addMuonTrack(antimu);
    }//for all mu+

    // get sum length of all muon tracks
    double totLen = 0;
    for (size_t t=0; t < _allMuonTracksInTPC.size(); t++)
      totLen += getLen(_allMuonTracksInTPC.at(t));

    _hAllMuonTrackLen->Fill(totLen/100.);
    std::cout << "Sum Muon Track Length in meters: " << totLen/100. << std::endl;
    
    for (size_t j=0; j < result.size(); j++){
    
      mcpart part = my_mcpart->at(_MCgetter.searchParticleMap( result.at(j).getNodeIndex() ));

      // make sure above energy cut
      if (part.Trajectory().at(0).E() > _Ecut){
	
	Reset();
	
	_run = my_mcpart->run() ;
	_subrun = my_mcpart->subrun();
	_event = my_mcpart->event_id(); 
	
	_PDG = part.PdgCode() ;
	_trackID = part.TrackId() ;
	_process= part.Process() ;
	
	_X = part.Trajectory().at(0).X() ;
	_Y = part.Trajectory().at(0).Y() ;
	_Z = part.Trajectory().at(0).Z() ;
	_T = part.Trajectory().at(0).T() ;
	
	_Px = part.Trajectory().at(0).Px() ;
	_Py = part.Trajectory().at(0).Py() ;
	_Pz = part.Trajectory().at(0).Pz() ;
	_E  = part.Trajectory().at(0).E() ;

	//Create mcshower-like things
	if( part.PdgCode() == -11 && part.Process() == "compt")
	  _PDG = 3 ;
	
	if( part.PdgCode() == -11 && part.Process() == "conv"){
	  _PDG = 4 ;
	  
	  if (_MCgetter.searchParticleMap(result.at(j).getParentId()) >= 0){
	    
	    mcpart mother = my_mcpart->at(_MCgetter.searchParticleMap( result.at(j).getParentId() ));
	    
	    _E = mother.Trajectory().at(mother.Trajectory().size()-2).E() ;
	    _Px = mother.Trajectory().at(mother.Trajectory().size()-2).Px() ;
	    _Py = mother.Trajectory().at(mother.Trajectory().size()-2).Py() ;
	    _Pz = mother.Trajectory().at(mother.Trajectory().size()-2).Pz() ;
	    
	  }
	}//if positron, conv
	
	geoalgo::Point_t vtx(_X, _Y, _Z);
	geoalgo::Vector_t mom(_Px, _Py, _Pz);

	if(_TpcBox.Contain(vtx))
	  _inActiveVolume = 1 ; 
	else
	  _inActiveVolume = 0 ;

	if (_inActiveVolume){
	  double muonDist, muonIP, distToIP;
	  getNearestMuonCutParams(vtx, mom, muonDist, muonIP, distToIP, 0);
	  _minMuDist = muonDist;
	  _minMuIP = muonIP;
	  _distToIP = distToIP;

	  geoalgo::HalfLine_t sDir(vtx,mom);
	  
	  _distAlongTraj = vtx.Dist(_iAlgo.Intersection(_TpcBox,sDir));
	  _distBackAlongTraj = vtx.SqDist(_iAlgo.Intersection(_TpcBox,sDir,true));
	}	
	
	// get parent information
	if (_MCgetter.searchParticleMap(result.at(j).getParentId()) >= 0){
	  mcpart mother = my_mcpart->at(_MCgetter.searchParticleMap( result.at(j).getParentId() ));
	  
	  _parentPDG = mother.PdgCode() ;
	  
	  _parentX = mother.Trajectory().at(0).X();
	  _parentY = mother.Trajectory().at(0).Y();
	  _parentZ = mother.Trajectory().at(0).Z();
	  _parentT = mother.Trajectory().at(0).T();
	  
	  _parentPx = mother.Trajectory().at(0).Px();
	  _parentPy = mother.Trajectory().at(0).Py();
	  _parentPz = mother.Trajectory().at(0).Pz();
	  _parentE  = mother.Trajectory().at(0).E();
	  
	  geoalgo::Point_t pVtx(_parentX, _parentY, _parentZ); 
	  geoalgo::Vector_t pMom(_parentPx, _parentPy, _parentPz);
	  
	  if(_TpcBox.Contain(pVtx))
	    _parentInActiveVolume = 1; 
	  else
	    _parentInActiveVolume = 0;
	  
	}// if mother exists
	
	
	// get ancestor information
	if (_MCgetter.searchParticleMap(result.at(j).getAncestorId()) >= 0){
	  mcpart ancestor = my_mcpart->at(_MCgetter.searchParticleMap( result.at(j).getAncestorId() ));

	  if (_inActiveVolume){
	    double muonDistExceptAncestor, muonIPExceptAncestor, distToIPExceptAncestor;
	    getNearestMuonCutParams(vtx, mom, muonDistExceptAncestor, muonIPExceptAncestor,
				     distToIPExceptAncestor, result.at(j).getAncestorId());
	    _minMuDistExceptAncestor = muonDistExceptAncestor;
	    _minMuIPExceptAncestor = muonIPExceptAncestor;
	    _distToIPExceptAncestor = distToIPExceptAncestor;
	  }
	  
	  _ancestorPDG = ancestor.PdgCode() ;
	  
	  _ancestorX = ancestor.Trajectory().at(0).X();
	  _ancestorY = ancestor.Trajectory().at(0).Y();
	  _ancestorZ = ancestor.Trajectory().at(0).Z();
	  _ancestorT = ancestor.Trajectory().at(0).T();
	  
	  _ancestorPx = ancestor.Trajectory().at(0).Px();
	  _ancestorPy = ancestor.Trajectory().at(0).Py();
	  _ancestorPz = ancestor.Trajectory().at(0).Pz();
	  _ancestorE  = ancestor.Trajectory().at(0).E();

	  geoalgo::Point_t aVtx(_ancestorX, _ancestorY, _ancestorZ);
	  geoalgo::Vector_t aMom(_ancestorPx, _parentPy, _ancestorPz);
	  
	  if(_TpcBox.Contain(aVtx))
	    _ancestorInActiveVolume = 1;  
	  else
	    _ancestorInActiveVolume = 0;

	  // get full ancestor trajectory
	  geoalgo::Trajectory_t ancTraj = getMuonTrack(ancestor);

	  if (ancTraj.size() > 1){
	    
	    //if ancestor is pi+/pi-/mu+/mu-/proton/e+/e- then find distance to track
	    if ( (abs(_ancestorPDG) == 11) or (abs(_ancestorPDG) == 13) or
		 (abs(_ancestorPDG == 211)) or (_ancestorPDG == 2212) ){
	      
	      double shrP = sqrt( mom.at(0)*mom.at(0) + mom.at(1)*mom.at(1) + mom.at(2)*mom.at(2) );
	      std::vector<double> shrDir = { mom.at(0)/shrP, mom.at(1)/shrP, mom.at(2)/shrP };
	      
	      // use shower start & momentum to define a segment
	      // which starts 3 meters before start point
	      // and ends 10 cm after start point
	      // aligned with shr momentum. Use for Impact Parameter
	      geoalgo::Point_t shrOrigin(vtx.at(0)-300*shrDir.at(0),
					 vtx.at(1)-300*shrDir.at(1),
					 vtx.at(2)-300*shrDir.at(2));
	      
	      geoalgo::Point_t shrEnd(vtx.at(0)+10*shrDir.at(0),
				      vtx.at(1)+10*shrDir.at(1),
				      vtx.at(2)+10*shrDir.at(2));
	      
	      geoalgo::LineSegment_t shrSeg(shrOrigin,shrEnd);
	      
	      //used for PoCA
	      geoalgo::Point_t c1(3);
	      geoalgo::Point_t c2(3);
	      
	      _ancDist = _dAlgo.SqDist(vtx,ancTraj);
	      _ancIP = _dAlgo.SqDist(ancTraj,shrSeg,c1,c2);
	      _ancToIP = sqrt ( (c2.at(0)-vtx.at(0))*(c2.at(0)-vtx.at(0)) +
				(c2.at(1)-vtx.at(1))*(c2.at(1)-vtx.at(1)) +
				(c2.at(2)-vtx.at(2))*(c2.at(2)-vtx.at(2)) );
	      
	    }//if PDG is charged
	  }//if trajectory of ancestor > 1 in length
	}// if ancestor exists
	
	
	
	
	_ana_tree->Fill();

      }	// if above energy cut

    } // for all particles
    
    
    return true;
  }
  
  void CosmicsBackground::PrepareTTree() {
    
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

	  _ana_tree->Branch("_inActiveVolume",&_inActiveVolume,"inActiveVolume/I");
	  
	  _ana_tree->Branch("_distAlongTraj",&_distAlongTraj,"distAlongTraj/D") ;
	  _ana_tree->Branch("_distBackAlongTraj",&_distBackAlongTraj,"distBackAlongTraj/D") ;

	  _ana_tree->Branch("_minMuDist",&_minMuDist,"minMuDist/D");
	  _ana_tree->Branch("_minMuIP",&_minMuIP,"minMuIP/D");
	  _ana_tree->Branch("_distToIP",&_distToIP,"distToIP/D");

	  _ana_tree->Branch("_minMuDistExceptAncestor",&_minMuDistExceptAncestor,"minMuDistExceptAncestor/D");
	  _ana_tree->Branch("_minMuIPExceptAncestor",&_minMuIPExceptAncestor,"minMuIPExceptAncestor/D");
	  _ana_tree->Branch("_distToIPExceptAncestor",&_distToIPExceptAncestor,"distToIPExceptAncestor/D");


	  ////ANCESTOR INFO
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

	  _ana_tree->Branch("_parentInActiveVolume",&_parentInActiveVolume,"parentInActiveVolume/I");

	  _ana_tree->Branch("_ancDist",&_ancDist,"ancDist/D");
	  _ana_tree->Branch("_ancIP",&_ancIP,"ancIP/D");
	  _ana_tree->Branch("_ancToIP",&_ancToIP,"ancToIP/D");



	  ////PARENT INFO
	  //
	  _ana_tree->Branch("_ancestorPDG",&_ancestorPDG,"ancestorPDG/I");
	  _ana_tree->Branch("_ancestorX",&_ancestorX,"ancestorX/D");
	  _ana_tree->Branch("_ancestorY",&_ancestorY,"ancestorY/D");
	  _ana_tree->Branch("_ancestorZ",&_ancestorZ,"ancestorZ/D");
	  _ana_tree->Branch("_ancestorT",&_ancestorT,"ancestorT/D");

	  _ana_tree->Branch("_ancestorPx",&_ancestorPx,"ancestorPx/D");
	  _ana_tree->Branch("_ancestorPy",&_ancestorPy,"ancestorPy/D");
	  _ana_tree->Branch("_ancestorPz",&_ancestorPz,"ancestorPz/D");
	  _ana_tree->Branch("_ancestorE",&_ancestorE,"ancestorE/D");

	  _ana_tree->Branch("_ancestorInActiveVolume",&_ancestorInActiveVolume,"ancestorInActiveVolume/I");


//	  _ana_tree->Branch("MuonTraj",&MuonTraj) ;
	
	}

 }

  void CosmicsBackground::addMuonTrack(const mcpart& part){

    mctrajectory traj = part.Trajectory();
    geoalgo::Trajectory_t thistrack(0,3);

    for (size_t i=0; i < traj.size(); i++){

      // check if point in TPC
      if ( (traj.at(i).X() > 0) and (traj.at(i).X() < 2*detHalfWidth) and
	   (traj.at(i).Y() > -detHalfHeight) and (traj.at(i).Y() < detHalfHeight) and
	   (traj.at(i).Z() > 0) and (traj.at(i).Z() < detLength) )
	thistrack.push_back( geoalgo::Point_t(traj.at(i).X(), traj.at(i).Y(), traj.at(i).Z()) );
				     
    }//for all points

    if (thistrack.size() > 0){
      _allMuonTracksInTPC.push_back(thistrack);
      _allMuonTracksIDs.push_back(part.TrackId());
    }

    return;
  }



  geoalgo::Trajectory_t CosmicsBackground::getMuonTrack(const mcpart& part){

    mctrajectory traj = part.Trajectory();
    geoalgo::Trajectory_t thistrack(0,3);

    for (size_t i=0; i < traj.size(); i++){

      // check if point in TPC
      if ( (traj.at(i).X() > 0) and (traj.at(i).X() < 2*detHalfWidth) and
	   (traj.at(i).Y() > -detHalfHeight) and (traj.at(i).Y() < detHalfHeight) and
	   (traj.at(i).Z() > 0) and (traj.at(i).Z() < detLength) )
	thistrack.push_back( geoalgo::Point_t(traj.at(i).X(), traj.at(i).Y(), traj.at(i).Z()) );
				     
    }//for all points
    
    return thistrack;
  }


void CosmicsBackground::Clear() {

  }

void CosmicsBackground::Reset(){


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

   _minMuDistExceptAncestor = -1;
   _minMuIPExceptAncestor = -1;
   _distToIPExceptAncestor = -1;

   _ancDist = -1;
   _ancIP = -1;
   _ancToIP = -1;

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


   //ancestor
   _ancestorPDG = -1;
   _ancestorX   = -9999999;
   _ancestorY   = -9999999;
   _ancestorZ   = -9999999;
   _ancestorT   = -9999999;

   _ancestorPx  = -9999999;
   _ancestorPy  = -9999999;
   _ancestorPz  = -9999999;
   _ancestorE   = -9999999;


   _inActiveVolume = -99 ;
   _parentInActiveVolume = -99 ;
   _ancestorInActiveVolume = -99 ;

 }

bool CosmicsBackground::finalize() {

  if(_fout) {

    _fout->cd();
    _hAllMuonTrackLen->Write();
    if(_ana_tree) 
      _ana_tree->Write();
  }
  
  else
    print(larlite::msg::kERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
  
  delete _ana_tree;
  
  return true;
  }


  double CosmicsBackground::getLen(geoalgo::Trajectory_t& muonTrack){
    
    return muonTrack.Length();
  }


  void CosmicsBackground::getNearestMuonCutParams( const geoalgo::Point_t& shrStart,
						   const geoalgo::Vector_t& shrMom,
						   double &muonDist,
						   double &muonIP,
						   double &distToIP, int ancestorTrackID){

    // use shower's start point + direction and list of muon trajectory to find:
    // 1) PoCA and PoCA distance to shower start point w/ PoCA GeoAlgo class
    // 2) MuonDist for cylinder cut

    // initialize values for PoCA cut
    double minIP = 10000;
    geoalgo::Point_t IP_onShr(3);
    geoalgo::Point_t c1(3);
    geoalgo::Point_t c2(3);
    // initialize values for muon proximity
    double minDist = 10000;

    geoalgo::Vector_t shrP = shrMom/shrMom.Length();
    
    // use shower start & momentum to define a segment
    // which starts 3 meters before start point
    // and ends 10 cm after start point
    // aligned with shr momentum. Use for Impact Parameter
    geoalgo::LineSegment_t shrSeg(shrStart-shrP*300,shrStart+shrP*10);
    
    // loop over all muon tracks and calculate value per-muon
    for (size_t u=0; u < _allMuonTracksInTPC.size(); u++){

      if ( (_allMuonTracksInTPC.at(u).size() > 1) and (_allMuonTracksIDs.at(u) != ancestorTrackID) ){
	
	// distance to muon track
	double tmpDist = _dAlgo.SqDist(shrStart,_allMuonTracksInTPC.at(u));
	// Impact parameter
	double tmpIP = _dAlgo.SqDist(_allMuonTracksInTPC.at(u),shrSeg,c1,c2);
	
	if (tmpDist < minDist) { minDist = tmpDist; }
	if (tmpIP < minIP) { 
	  minIP = tmpIP; 
	  IP_onShr = c2;
	}// if best muon so far
	
      }// if the muon's track is > 1
      
    }// for all muon tracks

    if (minDist != 10000){
      distToIP = IP_onShr.Dist(shrStart);
      muonDist = sqrt(minDist);
      muonIP = sqrt(minIP);
      //need to figure out if IP point is "before" or "after" start point w.r.t. momentum direction
      if (distToIP > 0.001 ){
	double dot = ( (IP_onShr-shrStart)/((IP_onShr-shrStart).Length() ))*shrP;
	if (dot == 1 ) { distToIP *= -1; }
      }
    }

    return;
  }
  
  
}
#endif
