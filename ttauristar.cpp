/**
 * \file ttauristar.cpp
 *
 * \authors Cynthia Yan
 *
 * \brief Implements the TTauriStar class.
 */

#include "ttauristar.hpp"

using namespace std;

double const TIMESTEP = 0.004;       /// timestep
double const ALPHA = 0.01;           /// viscocity parameter
double const BETA = 1.35;            /// R_M/R_A
double const GAMMA = 1.0;            /// R_C/R_M
//double const SHAPEFACTOR = 0.17;     /// f in the rotational inertia I=fMR^2
double const PROPEFF = 0.3;          /// propeller coefficient
double const CRITICALDENSITY = 6e-4; /// accretion disk density cutoff 6e-4
double const PROPSTARTTIME = 0.05;   /// simulation starts at this time
double const TURNONTIME = 0.05;
// double const PROPTIMESPREAD = 0.002; /// introduce randomness for PROPSTARTTIME
double const BFIELD = 1.67;           /// fieldstrength is set to be a constant
double const DELTAM = 0.01;            /// deltam to determine when to stop

TTauriStar::TTauriStar(vector<vector<double>> cmktable, 
	double mass, double age, double massdotfactor)
    :cmktable_(cmktable), mass_(mass), age_(age), massdotfactor_(massdotfactor)
{
	// set propeller endtime to be the same as starttime
	propendtime_ = PROPSTARTTIME;
	acceff_ = 1;
	massi_ = 0;
}

double TTauriStar::calculatemassdot()
{
	// massdot has unit M_sun/yr
	return 7.0e-8*massdotfactor_*pow(age_,-2.1)*pow(mass_,2.43);
}

double TTauriStar::calculateradius()
{
	size_t index1 = 0;
	size_t index3 = 0;

	// find masslower and massupper such that masslower <= mass < massupper
	double masslower = cmktable_[0][0];
	double massupper = cmktable_[0][0];
	while (cmktable_[0][index3] <= mass_) {
        if (cmktable_[0][index3] > masslower) {
        	index1 = index3;
        	masslower = cmktable_[0][index3];
        }
		++index3;
	}
	massupper = cmktable_[0][index3];
	// find agelower1 and ageupper1 such that agelower1 <= age < ageupper1
	double agelower1 = cmktable_[1][index1];
	double ageupper1 = cmktable_[1][index1];
	size_t index2 = index1;

	while (cmktable_[1][index2] <= age_) {
		++index2;
	}
	index1 = index2 - 1;
	agelower1 = cmktable_[1][index1];
	ageupper1 = cmktable_[1][index2];

	// find agelower2 and ageupper2 such that agelower2 <= age < ageupper2
	double agelower2 = cmktable_[1][index3];
	double ageupper2 = cmktable_[1][index3];
	size_t index4 = index3;
	while (cmktable_[1][index4] <= age_) {
		++index4;
	}
	index3 = index4 - 1;
	agelower2 = cmktable_[1][index3];
	ageupper2 = cmktable_[1][index4];

	// four coefficients
	double coeff1 = (massupper-mass_)*(ageupper1-age_)/(massupper-masslower)/(ageupper1-agelower1);
	double coeff2 = (massupper-mass_)*(age_-agelower1)/(massupper-masslower)/(ageupper1-agelower1);
	double coeff3 = (mass_-masslower)*(ageupper2-age_)/(massupper-masslower)/(ageupper2-agelower2);
	double coeff4 = (mass_-masslower)*(age_-agelower2)/(massupper-masslower)/(ageupper2-agelower2);

	// return the interpolated radius
	return coeff1*cmktable_[2][index1]+coeff2*cmktable_[2][index2]+coeff3*cmktable_[2][index3]+coeff4*cmktable_[2][index4];
}

double TTauriStar::calculatebfield()
{
	// turn on a constant dipole magnetic field at some input time
	
	if (age_ > TURNONTIME) {
		return BFIELD;
	} else {
		return 0;
	}
}

double TTauriStar::calculaterm()
{
	// radius at which rom pressure equals magnetic pressure
	return radius_*7.1883*BETA*pow(mass_/0.5,-1./7.)*pow(bfield_,4./7.)*pow(radius_/2.,5./7.)*pow(massdot_/1.e-8,-2./7.);
}

double TTauriStar::calculatediskdensity()
{
	// ratio of moment of inertia and mass*radius^2
	double f = pow(1.0-pow(radius_/rm_,0.5),1.0/4.0);
	// density of accretion disk at rm
	return 5.2*pow(ALPHA,-4.0/5.0)*pow(massdot_*6.307e9,0.7)*pow(mass_,0.25)*pow(rm_*6.996,-0.75)*pow(f,14.0/5.0);
}

void TTauriStar::calculatemasses()
{
	// clear vectors involved
	masses_.clear();
	ages_.clear();
	acceffs_.clear();
	// initial values
	masses_.push_back(mass_);
	ages_.push_back(age_);
	acceffs_.push_back(acceff_);
	// go backwards in time
	while (age_ > PROPSTARTTIME) {
		// calculate mass accretion rate
		massdot_ = calculatemassdot();
		// determine the acceff
		if (age_ > propendtime_) {
			// no accretion happens at propeller phase
			acceff_ = 1;
		} else {
			acceff_ = 0;
		}
		// calculate new mass
		mass_ -= 1.0e6*massdot_*acceff_*TIMESTEP;	
		// push_back the mass into the masses vector
		// will reversed after all the push_backs are done for efficiency
		masses_.push_back(mass_);
		// calculate new age
		age_ -= TIMESTEP;
		// push_back the age into the ages vector
		// will also be reversed later
		ages_.push_back(age_);
		acceffs_.push_back(acceff_);
	}
	// reverse the two vectors to be in forward time order
	reverse(masses_.begin(),masses_.end());
	reverse(ages_.begin(),ages_.end());
	reverse(acceffs_.begin(),acceffs_.end());
}

void TTauriStar::calculateperiods()
{
	// clear vectors involved
	periods_.clear();
	// initialze massi_
	massi_ = masses_[0];
	// go forward in time
	for (size_t i = 0; i < ages_.size(); ++i) {
		// retrieve mass and age stored in vectors
		mass_ = masses_[i];
		age_ = ages_[i];
		acceff_ = acceffs_[i];
		// update data members
		massdot_ = calculatemassdot();
		radius_ = calculateradius();
		bfield_ = calculatebfield();
		rm_ = calculaterm();
		diskdensity_ = calculatediskdensity();
		// calculate the period at rm
	    double periodrm = 0.1159*pow(rm_,3./2.)*pow(mass_,-1./2.);
	    // 1. spin at break-up period
	    if (i == 0) {
	    	period_ = 0.1159*pow(radius_,3./2.)*pow(mass_,-1./2.);	
	    	cout << 1 << endl;
		// 2. spin down due to propeller effect		
		} else if (period_ < periodrm) {
			period_ += TIMESTEP*PROPEFF*0.972*pow(BETA,-3.)*pow(period_,2.)*pow(mass_,-4./7.)*pow(bfield_,2./7.)*pow(radius_,-8./7.)*pow(massdot_/1.e-8,6./7.);
			// keep track of the propeller endtime
			propendtime_ = age_;
			cout << 2 << endl;
		// 3. disk-locked
		} else if (diskdensity_ > CRITICALDENSITY) {
			period_ = 8.*pow(GAMMA*BETA/0.9288,3./2.)*pow(massdot_/1.0e-8,-3./7.)*pow(mass_/0.5,-5./7.)*pow(radius_/2.,18./7.)*pow(bfield_,6./7.);
		    cout << 3 << endl;
		// 4. unlocked
		} else {
			period_ += TIMESTEP*acceff_*period_*massdot_/mass_*(1-pow(rm_,0.5)*pow(radius_,-2.)*period_*8.08e6);
		    cout << 4 << endl;
		}
		// store the period
		periods_.push_back(period_);
		//radii_.push_back(radius_);
		//rms_.push_back(rm_);
	}
	cout << acceff_ << endl;
}

void TTauriStar::update()
{
	// keep track of the number of iterations
	int i = 1;
	calculatemasses();
	while (abs((masses_[0]-massi_)/masses_[0]) > DELTAM) {	
	    calculateperiods();
	    calculatemasses();
	    ++i;
	}
	calculateperiods();
	cout << "iterate " << i << " times" << endl;	
}

vector<vector<double>> TTauriStar::getvectors(string vectorname1, string vectorname2)
{
	vector<double> vector1;
	vector<double> vector2;
	if (vectorname1 == "Age") {
		vector1 = ages_;
	}
	if (vectorname1 == "Mass") {
		vector1 = masses_;
	}
	if (vectorname1 == "Period") {
		vector1 = periods_;
	}
	if (vectorname1 == "Acceff") {
		vector1 = acceffs_;
	}
	if (vectorname1 == "Radius") {
		vector1 = radii_;
	}
	if (vectorname1 == "RM") {
		vector1 = rms_;
	}
	if (vectorname2 == "Age") {
		vector2 = ages_;
	}
	if (vectorname2 == "Mass") {
		vector2 = masses_;
	}
	if (vectorname2 == "Period") {
		vector2 = periods_;
	}
	if (vectorname2 == "Acceff") {
		vector2 = acceffs_;
	}
	if (vectorname2 == "Radius") {
		vector2 = radii_;
	}
	if (vectorname2 == "RM") {
		vector2 = rms_;
	}
	vector<vector<double>> output;
	output.push_back(vector1);
	output.push_back(vector2);
	return output;
}
