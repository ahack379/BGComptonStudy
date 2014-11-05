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
	std::vector<double> _positronVtx ;

	//Info for _gamma_tree2
	double _rand ;
	double _prob ;
	double _checkGammaPP ; 


	std::vector<double> energy ;
	std::vector<double> pp ;

	//Muon Tree
	double _timeMuonEnter ; 
	double _timeDiff ;


	//Info for _ana_tree
	double _timeElecEnter ;  //Time the compton scattered electron enters the TPC

	
	double _inVolElecX ;
	double _inVolElecY ;
	double _inVolElecZ ;
    double _energyElec ;

	double _dist_ToWall ;
	double _dist_AlongTraj ;
	double _dist_BackAlongTraj ;
	
	std::vector<double>  _elecVtx ;
	std::vector<double> _elecMom ;

	//Info for _gamma_tree 
	double _pdgCode ;
	double _energyGammaTotal ;
	double _inVolGammaX ;
	double _inVolGammaY ;
	double _inVolGammaZ ;
	double _timing ;

	std::vector<double>  _gammaVtx ;

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
