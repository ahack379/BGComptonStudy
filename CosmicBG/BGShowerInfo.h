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
#include "DistanceAlgo.h"
#include "LArUtil/Geometry.h"
#include "IntersectAlgo.h"
#include <vector>
#include <string>

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

        // GeoAlgo for Distance Algos
        geoalgo::DistanceAlgo _dAlgo;
	// geoalgo for Intersection Algos
	geoalgo::IntersectAlgo _iAlgo;
	
	// TPC AABox object
	geoalgo::AABox _TpcBox;

	//ana_tree currently deals with electrons 
	//from compton scatters
	TTree * _ana_tree ;
	TTree * _pp_tree ;
	TTree * _gamma_tree ;
	TTree * _gamma_tree2 ;
	TTree * _muon_tree ;

	//Info for PP tree
    double _energyGammaBegin;
    double _energyGammaEnd;

	double _inVolPosX ;
	double _inVolPosY ;
	double _inVolPosZ ;

	double _pxPP ;
	double _pyPP ;
	double _pzPP ;
	double _energyPP ;
	
	double _gammaPPX ;
	double _gammaPPY ;
	double _gammaPPZ ;

	double _dist_BackAlongTrajPP ;
	std::vector<double> _positronVtx ;
	std::vector<double> _ppMom ;


	//Info for _ana_tree
	double _timeElecEnter ;  //Time the compton scattered electron enters the TPC
	
	double _inVolElecX ;
	double _inVolElecY ;
	double _inVolElecZ ;

	double _pxElec ;
	double _pyElec ;
	double _pzElec ;
    double _energyElec ;

	double _gammaCompX ;
	double _gammaCompY ;
	double _gammaCompZ ;

	double _dist_ToWall ;
	double _dist_AlongTraj ;
	double _dist_BackAlongTraj ;
	
	std::vector<double>  _elecVtx ;
	std::vector<double> _elecMom ;

	//Info for _gamma_tree 
	std::string _process ;
	double _energyMother ;
	double _parents ;
	double _parentsCut ;
	double _pdgCode ;
	double _energyGammaTotal ;
	double _inVolGammaX ;
	double _inVolGammaY ;
	double _inVolGammaZ ;
	double _timing ;

	std::vector<double>  _gammaVtx ;

	//Info for _gamma_tree2
	double _rand ;
	double _prob ;
	double _checkGammaPP ; 

	//Muon Tree
	double _timeMuonEnter ; 
	double _timeDiff ;

	std::vector<double> energy ;
	std::vector<double> pp ;


	int _count0 ;
	int _count1 ;


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
