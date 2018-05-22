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
double const CRITICALDENSITY = 250; /// accretion disk density cutoff 6e-4
double const PROPSTARTTIME = 0.05;   /// simulation starts at this time
double const TURNONTIME = 0.05;
// double const PROPTIMESPREAD = 0.002; /// introduce randomness for PROPSTARTTIME
double const BFIELD = 1.67;           /// fieldstrength is set to be a constant
double const DELTAM = 0.0002;            /// deltam to determine when to stop

TTauriStar::TTauriStar(vector<vector<double>> cmktable, 
	double mass, double age, double massdotfactor)
    :cmktable_(cmktable), mass_(mass), age_(age), massdotfactor_(massdotfactor)
{
	// iniliatize validity
	if (mass > 3) {
		valid_ = false;
	} else {
		valid_ = true;
	}
	// initialize mass0_
	mass0_ = mass;
	// initialize mass2_
	mass2_ = 0;
	// set propeller endtime to be the same as starttime
	propendtime_ = PROPSTARTTIME;
	acceff_ = 1.0;
	// starting values of age and acceff
	ages_.push_back(age_);
	acceffs_.push_back(acceff_);
	// go backwards in time
	while (age_ > PROPSTARTTIME) {
		// calculate new age
		age_ -= TIMESTEP;
		// push_back the age into the ages vector
		// will also be reversed later
		ages_.push_back(age_);
		acceffs_.push_back(1.0);
	}
	// reverse the two vectors to be in forward time order
	reverse(ages_.begin(),ages_.end());
	reverse(acceffs_.begin(),acceffs_.end());
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
	return 8.79e6*pow(ALPHA,-4.0/5.0)*pow(massdot_,0.7)*pow(mass_,0.25)*pow(rm_,-0.75)*pow(f,14.0/5.0);
}

void TTauriStar::calculatemasses()
{
	// clear vectors involved
	masses_.clear();
	// initial values
	masses_.push_back(mass0_);
	// go backwards in time
	for (size_t i = ages_.size() - 1; i >= 1; --i) {
		// retrieve age and acceff
		age_ = ages_[i];
		acceff_ = acceffs_[i];
		// calculate mass accretion rate
		massdot_ = calculatemassdot();
		// calculate new mass
		mass_ -= 1.0e6*massdot_*acceff_*TIMESTEP;	
		// push_back the mass into the masses vector
		// will reversed after all the push_backs are done for efficiency
		masses_.push_back(mass_);
	}
	// reverse the two vectors to be in forward time order
	reverse(masses_.begin(),masses_.end());
}

void TTauriStar::calculateperiods()
{
	// clear vectors involved
	periods_.clear();
	massdots_.clear();
	radii_.clear();
	rms_.clear();
	diskdensities_.clear();
	acceffs_.clear();
	// initialze mass2_
	mass2_ = masses_[0];
	// initial phase
	size_t phase = 1;
	// go forward in time
	for (size_t i = 0; i < ages_.size(); ++i) {
		// retrieve mass and age stored in vectors
		mass_ = masses_[i];
		// cannot handle mass > 3 solar mass
		if (mass_ > 3) {
			valid_ = false;
			break;
		}
		age_ = ages_[i];
		acceff_ = 1;
		// update data members
		massdot_ = calculatemassdot();
		double radius = radius_;
		radius_ = calculateradius();
		bfield_ = calculatebfield();
		rm_ = calculaterm();
		diskdensity_ = calculatediskdensity();
		// calculate the period at rm
	    double periodrm = 0.1159*pow(rm_,3./2.)*pow(mass_,-1./2.);
	    // 1. spin at break-up period
	    if (i == 0 && phase <= 1) {
	    	period_ = 0.1159*pow(radius_,3./2.)*pow(mass_,-1./2.);	
	    	// doesn't accrete
	    	acceff_ = 0;
		// 2. spin down due to propeller effect		
		} else if (period_ < periodrm && phase <= 2) {
			period_ += TIMESTEP*PROPEFF*0.972*pow(BETA,-3.)*pow(period_,2.)*pow(mass_,-4./7.)*pow(bfield_,2./7.)*pow(radius_,-8./7.)*pow(massdot_/1.e-8,6./7.);
			// keep track of the propeller endtime
			propendtime_ = age_;
			// doesn't accrete
			acceff_ = 0;
			phase = 2;
		// 3. disk-locked
		} else if (diskdensity_ > CRITICALDENSITY && phase <= 3) {
			period_ = 8.*pow(GAMMA*BETA/0.9288,3./2.)*pow(massdot_/1.0e-8,-3./7.)*pow(mass_/0.5,-5./7.)*pow(radius_/2.,18./7.)*pow(bfield_,6./7.);
		    phase = 3;
		// 4. unlocked
		} else {
			// G in units of (solar radius^3)/(day^2 solarmass) G = 2937.5
			period_ += period_*2*(radius_-radius)/radius_
			    +TIMESTEP*acceff_*period_*massdot_/mass_
				-50.74*TIMESTEP*acceff_*pow(period_,2)*massdot_/pow(mass_,0.5)*pow(rm_,0.5)*pow(radius_,-2.);
		    phase = 4;
		}
		// calculate mass moving forward
		if (i < ages_.size() - 1) {
			mass2_ += 1.0e6*massdot_*acceff_*TIMESTEP;
		}	
		// store the period
		periods_.push_back(period_);
		massdots_.push_back(massdot_);
		radii_.push_back(radius_);
		rms_.push_back(rm_);
		diskdensities_.push_back(diskdensity_);
		acceffs_.push_back(acceff_);
	}
}

double TTauriStar::update()
{
	if (valid_) {
		// keep track of the number of iterations
		int i = 0;
		while (abs((mass2_-mass0_)/mass0_) > DELTAM && i < 20) {	
			calculatemasses();
		    calculateperiods();	
		    // cout << "mass2" << mass2_ << endl;
		    // cout << "mass" << mass_ << endl;    
		    ++i;
		}
		cout << "iterate " << i << " times" << endl;
		return period_;	
	} else {
		return 0;
	}
	
}

vector<double> TTauriStar::getvector(int n)
{
	vector<double> output;
	if (n == 1) {
		output = ages_;
	}
	if (n == 2) {
		output = masses_;
	}
	if (n == 3) {
		output = periods_;
	}
	if (n == 4) {
		output = massdots_;
	}
	if (n == 5) {
		output = radii_;
	}
	if (n == 6) {
		output = rms_;
	}
	if (n == 7) {
		output = diskdensities_;
	}
	if (n == 8) {
		output = acceffs_;
	}

	return output;
}

string TTauriStar::getname(int n)
{
	string output;
	if (n == 1) {
		output = "Age";
	}
	if (n == 2) {
		output = "Mass";
	}
	if (n == 3) {
		output = "Period";
	}
	if (n == 4) {
		output = "Accretion rate";
	}
	if (n == 5) {
		output = "Radius";
	}
	if (n == 6) {
		output = "R_M";
	}
	if (n == 7) {
		output = "Disk density";
	}
	if (n == 8) {
		output = "Acceff";
	}

	return output;
}

string TTauriStar::getunit(int n)
{
	string output;
	if (n == 1) {
		output = "Myr";
	}
	if (n == 2) {
		output = "solar mass";
	}
	if (n == 3) {
		output = "day";
	}
	if (n == 4) {
		output = "solar mass/year";		
	}
	if (n == 5) {
		output = "solar radius";
	}
	if (n == 6) {
		output = "solar radius";
	}
	if (n == 7) {
		output = "g/cm^2";
	}
	if (n == 8) {
		output = "";
	}

	return output;
}

void TTauriStar::plot(int m, int n)
{
	vector<double> vector1 = getvector(m);
	vector<double> vector2 = getvector(n);
    FILE * temp1 = fopen("data.temp", "w");
    FILE* gp1=popen("gnuplot -persistent","w");
    for(size_t k=0;k<vector1.size();k++) {
        fprintf(temp1,"%f %f \n",vector1[k],vector2[k]);
    }
    fprintf(gp1, "%s%s %s %s%s\n", "set title \"",getname(n).data(),"vs",getname(m).data(),"\"");

    fprintf(gp1, "%s%s %s%s%s\n", "set xlabel \"",getname(m).data(),"(",getunit(m).data(),")\"");
    fprintf(gp1, "%s%s %s%s%s\n", "set ylabel \"",getname(n).data(),"(",getunit(n).data(),")\"");
    fprintf(gp1, "%s \n", "plot 'data.temp'");
}
