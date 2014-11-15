/**
 * \file CosmicsBackground.h
 *
 * \ingroup CosmicBG
 * 
 * \This class builds several trees to investigate pp and compton scattering from cosmic gammas. 
 *
 * @author arianaHackenburg
 */

/** \addtogroup CosmicBG

    @{*/

#ifndef LARLITE_COSMICSBACKGROUND_H
#define LARLITE_COSMICSBACKGROUND_H

#include "Analysis/ana_base.h"
#include "DistToBoxWall.h"
#include "TrajectoryInVolume.h"
#include "LArUtil/Geometry.h"
#include "SegmentPoCA.h"
#include "PointToLineDist.h"
#include "MCgetter.h"
#include <vector>
#include <string>

namespace larlite {

  class CosmicsBackground : public ana_base{
  
  public:

    CosmicsBackground(){ _name="CosmicsBackground"; _fout=0;};

    virtual ~CosmicsBackground(){};

    void PrepareTTree() ;

    void Clear() ;

    void Reset() ;   

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void addMuonTrack(mcpart *part);

    std::vector<std::vector<double> > getMuonTrack(mcpart *part);

    void setECut(double e) { _Ecut = e; }

    double getLen(std::vector<std::vector<double> > *muonTrack);

    /// Set MCgetter object
    void SetMCgetter(MCgetter mcgetter) { _MCgetter = mcgetter; }

    void getNearestMuonCutParams( std::vector<double> *shrStart,
				  std::vector<double> *shrMom,
				  double &muonDist,
				  double &muonIP,
				  double &distToIP,
				  int ancestorTrackID);
      

    protected:
	//ana_tree currently deals with electrons 
	//from compton scatters
	TTree * _ana_tree ;
		
	int _run ;
	int _subrun ;
	int _event ;

	std::string _process ;
	int _PDG ;
	int _trackID ;

	double _X ;
	double _Y ;
	double _Z ;
	double _T ;

	double _Px ;
	double _Py ;
	double _Pz ;
	double _E ;

	int _inActiveVolume ;

	double _distAlongTraj ;
	double _distBackAlongTraj ;

	double _minMuDist;
	double _minMuIP;
	double _distToIP;

	double _minMuDistExceptAncestor;
	double _minMuIPExceptAncestor;
	double _distToIPExceptAncestor;

	double _ancDist;
	double _ancIP;
	double _ancToIP;

	//Save info about parent as well
	int _parentPDG ;
	double _parentX ;
	double _parentY ;
	double _parentZ ;
	double _parentT ;

    double _parentPx;
	double _parentPy;
	double _parentPz ; 
	double _parentE ; 

	int _parentInActiveVolume ;


	//Save info about ancestor as well
	int _ancestorPDG ;
	double _ancestorX ;
	double _ancestorY ;
	double _ancestorZ ;
	double _ancestorT ;

    double _ancestorPx;
	double _ancestorPy;
	double _ancestorPz ; 
	double _ancestorE ; 

	int _ancestorInActiveVolume ;


	// MCgetter stuff
	MCgetter _MCgetter;

	//Not sure where this goes
//	int kCompton = 3 ;
//	int kPairProduction = 4; 

	geoalgo::DistToBoxWall 		 _showerObject;  
	geoalgo::TrajectoryInVolume	 _inVol ;  
	/// GeoAlg for PoCA cut
	geoalgo::SegmentPoCA _PoCA;
	/// GeoAlg for point to line dist
	geoalgo::PointToLineDist _pointDist;

	//Muon things
	std::vector<std::vector<std::vector<double> > > _allMuonTracksInTPC ; 
	// all tracks except ancestor
	std::vector<int> _allMuonTracksIDs ;

	// detector size
	double detHalfHeight;
	double detHalfWidth;
	double detLength;

	int _count0 ;


	double _Ecut;

	TH1D *_hAllMuonTrackLen;

  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
