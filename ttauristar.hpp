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
 *
 * \details The first parameter is a pointer matrix to the mass-age-radius table, 
 *          mass and age are the final values for a particular star 
 *          massdotfactor is the rendomization Mdot parameter. 
 */
class TTauriStar {
public:

	
	TTauriStar(vector<vector<double>> cmktable, double mass, double age, double massdotfactor);
	/**
	 * \brief update runs the simulation for each star until convergence.
	 *
	 * \param
	 * \returns period
	 */
	double update();

	/**
	 * \brief getvector fetches various evolutionary vectors for each star, 
	 *
	 * \param integer code of the vector you want
	 *
	 * \returns the vector you want
	 */ 
	vector<double> getvector(int n);

	/**
	 * \brief fetch names of these vectors for use in plotting
	 *
	 * \param integer code of the vector you want
	 *
	 * \returns a string telling the name
	 */
	string getname(int n);

	/**
	 * \brief fetch units of these vectors for use in plotting
	 *
	 * \param integer code of the vector you want
	 *
	 * \returns a string telling the unit
	 */
	string getunit(int n);

	/**
	 * \brief makes plots through gnuplot vector m vs. vector n is plotted.
	 *
	 * \param integer codes of the vectors you want
	 *
	 * \returns 
	 */
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
	double mass0_;              /// mass that we are trying to converge to, equal to the original input mass
	double mass2_;              /// mass calculated going forward to test for convergence against mass0
	double age_;
	double massdotfactor_;      /// equal to the original input parameter
	double massdot_;            /// M_sun/yr
	double period_;
	double radius_;
	double bfield_;
	double rm_;
	double diskdensity_;        /// surface density ar R_M
	double propendtime_;        /// end time for the propeller effect
	double acceff_;             /// fraction of deposition by accretion

	bool valid_;  /// if the mass is > 3, there are not values in the table, so we drop stars with those masses. 
	//This should be changed in the future with better data tables.
 	
	vector<double> ages_;
	vector<double> masses_;         /// mass of the protostar changed backward
	vector<double> massdots_;       /// mass accretion rate of the protostar
	vector<double> periods_;        /// period of the protostar
	vector<double> acceffs_;
	vector<double> radii_;
	vector<double> rms_;
	vector<double> diskdensities_;
};

#endif
