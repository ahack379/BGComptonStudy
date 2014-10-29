/**
 * \file BGShowerInfo.h
 *
 * \ingroup CosmicBG
 * 
 * \brief Class def header for a class BGShowerInfo
 *
 * @author arianaHackenburg
 */

/** \addtogroup CosmicBG

    @{*/

#ifndef LARLITE_BGSHOWERINFO_H
#define LARLITE_BGSHOWERINFO_H

#include "Analysis/ana_base.h"
#include "DistToBoxWall.h"
#include "TrajectoryInVolume.h"
#include <vector>

namespace larlite {

  class BGShowerInfo : public ana_base{
  
  public:

    BGShowerInfo(){ _name="BGShowerInfo"; _fout=0;};

    virtual ~BGShowerInfo(){};

    void PrepareTTree() ;

    void Clear() ;

    void Reset() ;   

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    protected:

	TTree * _ana_tree ;

	//0 for electron 1 for gamma
    bool _elecOrGamma ;
    bool _muonParent ;

    double _energyGamma;
    double _energyElec ;
    double _energyMuon ;    

	double _dist_ToWall ;
	double _dist_AlongTraj ;
	double _dist_BackAlongTraj ;

	double _inVolElecX ;
	double _inVolElecY ;
	double _inVolElecZ ;

	std::vector<double>  _elecVtx ;
	std::vector<double>  _gammaVtx ;
	std::vector<double>  _muonVtx ;
	
	std::vector<double> _elecMom ;
	std::vector<double> _gammaMom ;
	std::vector<double> _muonMom ;

	int _count0 ;
	int _count1 ;
	int _count2 ;

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
