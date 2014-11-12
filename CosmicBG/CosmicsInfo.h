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

	//Not sure where this goes
	int kCompton = 3 ;
	int kPairProduction = 4; 


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
