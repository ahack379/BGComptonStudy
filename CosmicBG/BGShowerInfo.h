/**
 * \file BGShowerInfo.h
 *
 * \ingroup CosmicBG
 * 
 * \This class builds several trees to investigate pp and compton scattering from cosmic gammas. 
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
	TTree * _pp_tree ;
	TTree * _gamma_tree ;
	TTree * _gamma_tree2 ;

    double _energyGammaBegin;
    double _energyGammaEnd;
    double _energyElec ;

	double _energyGammaTotal ;

	double _dist_ToWall ;
	double _dist_AlongTraj ;
	double _dist_BackAlongTraj ;
	
	double _inVolElecX ;
	double _inVolElecY ;
	double _inVolElecZ ;

	double _inVolGammaX ;
	double _inVolGammaY ;
	double _inVolGammaZ ;

	std::vector<double>  _elecVtx ;
	std::vector<double>  _gammaVtx ;
	std::vector<double>  _positronVtx ;
	
	std::vector<double> _elecMom ;
	std::vector<double> _positronMom ;

	int _count0 ;
	int _count1 ;

	std::vector<double> energy ;
	std::vector<double> pp ;

	double _rand ;
	double _prob ;

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
