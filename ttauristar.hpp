/**
 * \file ttauristar.hpp
 *
 * \authors Cynthia Yan
 *
 * \brief Declares the TTauriStar class.
 */

#ifndef TTAURISTAR_HPP_INCLUDED
#define TTAURISTAR_HPP_INCLUDED

#include <cstddef>
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;
/**
 * \class TTauriStar
 *
 * \brief Represents a T Tauri Star through out its evolution.
 */
class TTauriStar {
public:
	TTauriStar(vector<vector<double>> cmktable, double mass, double age, double massdotfactor);
	
	double update();
	vector<double> getvector(int n);
	string getname(int n);
	string getunit(int n);
	void plot(int m, int n);
private:
	double calculatemassdot();
	double calculateradius();
	double calculatebfield();
	double calculaterm();
	double calculatediskdensity();
	void calculatemasses();
	void calculateperiods();
	
	// TTauriStar data members
	vector<vector<double>> cmktable_;
	double mass_;
	double mass0_;              /// mass that we are trying to converging to
	double mass2_;              /// mass calculated going forward 
	double age_;
	double massdotfactor_;
	double massdot_;            /// M_sun/yr
	double period_;
	double radius_;
	double bfield_;
	double rm_;
	double diskdensity_;
	double propendtime_;
	double acceff_;             /// fraction of deposition by accretion

	bool valid_;
	
	vector<double> ages_;
	vector<double> masses_;         /// mass of the protostar
    vector<double> massdots_;         /// mass of the protostar
	vector<double> periods_;       /// period of the protostar
	vector<double> acceffs_;
	vector<double> radii_;
	vector<double> rms_;
	vector<double> diskdensities_;
};

#endif
