/**
 * \file CosmicsInfo.h
 *
 * \ingroup CosmicBG
 * 
 * \This class builds several trees to investigate pp and compton scattering from cosmic gammas. 
 *
 * @author arianaHackenburg
 */

/** \addtogroup CosmicBG

    @{*/

#ifndef LARLITE_COSMICSINFO_H
#define LARLITE_COSMICSINFO_H

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

  class CosmicsInfo : public ana_base{
  
  public:

    CosmicsInfo(){ _name="CosmicsInfo"; _fout=0;};

    virtual ~CosmicsInfo(){};

    void PrepareTTree() ;

    void Clear() ;

    void Reset() ;   

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void addMuonTrack(mcpart *part);

    void getNearestMuonCutParams( std::vector<double> *shrStart,
				  std::vector<double> *shrDir);
      

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

	bool _inActiveVolume ;

	double _distAlongTraj ;
	double _distBackAlongTraj ;

	double _minMuDist;
	double _minMuIP;
	double _distToIP;

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

	bool _parentInActiveVolume ;


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

	// detector size
	double detHalfHeight;
	double detHalfWidth;
	double detLength;

	int _count0 ;

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
