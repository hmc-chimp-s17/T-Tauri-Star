/**
 * \file onestar.cpp
 *
 * \authors Cynthia Yan
 *
 * \brief Provides the main() function for simulating one ttauri star
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <random>
#include <stdio.h>
#include "ttauristar.hpp"

using namespace std;
/**
 * \brief read cmk data file
 * 
 * \param fname      string representing the cmk data filename
 * \returns          a vector of vectors of doubles     
 */
vector<vector<double>> readcmk(string fname)
{
	ifstream inputFile(fname);

	if (!inputFile.good()) {
        throw invalid_argument( "Couldn't open cmk file for reading" );
    }

    // throw the first line
    string line;
    getline(inputFile, line);

    // allocate space on the heap to store the cmk data as a vector of vectors of doubles
    vector<vector<double>> cmktable;

    vector<double> cmkmasses;
    vector<double> cmkages;
    vector<double> cmkradii;

    const double PI = 3.141592653589793;
    const double SIGMAB = 5.67e-5;   // stefan-boltzmann constatn in erg/K^4 cm^2 s
    const double SOLARRADIUS = 7e10; // Solar radius in cm
    const double SOLARLUMINOSITY = 3.86e33;  // Solar luminosity in ergs per second
    // Solar tempertaure in Kelvins
    const double SOLARTEMPERATURE = pow(SOLARLUMINOSITY/(4*PI*pow(SOLARRADIUS,2)*SIGMAB),0.25);


    double mass;
    double logage;
    double age;
    double logL; // log(Luminosity)
    double logT; // log(Temperature)
    double radius;
    
    // read one line at a time
    while (getline(inputFile, line)){
    	stringstream lineStream(line);
        // read the four doubles separated by spaces 
        // and put them in the corresponding vector
        lineStream >> mass;
        lineStream >> logage;
        lineStream >> logL;
        lineStream >> logT;
        // calculate age
        age = pow(10.0,logage - 6.0);
        // calculate radius
        radius = pow(10,logL/2) * pow(pow(10,logT)/SOLARTEMPERATURE,-2);
        cmkmasses.push_back(mass);
        cmkages.push_back(age); 
        cmkradii.push_back(radius);  
    }
    // store mass, age, radius vectors in cmktable
    cmktable.push_back(cmkmasses);
    cmktable.push_back(cmkages);
    cmktable.push_back(cmkradii);

    // close file
    inputFile.close();
    return cmktable;
}

/**
 * \brief read cluster mass and age file
 * 
 * \param fname      string representing the cluster filename
 * \returns          a vector of vectors of doubles     
 */
vector<vector<double>> readcluster(string fname)
{
    ifstream inputFile(fname);

    if (!inputFile.good()) {
        throw invalid_argument( "Couldn't open cmk file for reading" );
    }

    // throw the first line
    string line;
    getline(inputFile, line);

    // allocate space on the heap to store the cmk data as a vector of vectors of doubles
    vector<vector<double>> startable;

    vector<double> logmasses;
    vector<double> logages;

    double mass;
    double age;
    
    // read one line at a time
    while (getline(inputFile, line)){
        stringstream lineStream(line);
        // read the four doubles separated by spaces 
        // and put them in the corresponding vector
        lineStream >> mass;
        lineStream >> age;
        // push back log10 of masses and ages
        logmasses.push_back(log10(mass));
        logages.push_back(log10(age)); 
    }
    // store mass, age, radius vectors in cmktable
    startable.push_back(logmasses);
    startable.push_back(logages);

    // close file
    inputFile.close();
    return startable;
}

/**
 * \brief generate the distribution corresponding to a cluster
 * 
 * \param fname      vectors of logmasses and logages
 * \returns          a distribution     
 */
vector<vector<double>> generatedistribution(vector<vector<double>> startable)
{
    // create mass bins
    vector<double> logmassbins{-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1};
    // number of mass bins
    size_t nummassbin = logmassbins.size() - 1;
    // initialize weights with zeros
    vector<double> logmassweights(nummassbin,0.0);
    // create age bins 
    vector<double> logagebins{-1.2,-1,-0.8,-0.4,-0.2,0,0.2,0.4};
    // number of age bins
    size_t numagebin = logagebins.size() - 1;
    // initilize ages with zeros
    vector<vector<double>> logageweights;
    for (size_t i = 0; i < nummassbin; ++i) {
        logageweights.push_back(vector<double>(numagebin,0.0));
    }
    // iterate through startable to get the weight
    // get logmass and logage from the table
    vector<double> logmasses = startable[0];
    vector<double> logages = startable[1];
    // the number of stars in the table
    size_t numstar = logmasses.size();
    // iteratre over the stars in the cluster
    for (size_t i = 0; i < numstar; ++i) {
        // find the interval that the logmass belongs to
        for (size_t j = 0; j < nummassbin; ++j) {
            if (logmassbins[j] <= logmasses[i] && logmasses[i] < logmassbins[j+1]) {
                // increase the mass count in the appropriate interval
                logmassweights[j] += 1.0;
                // look at its age
                for (size_t k = 0; k < numagebin; ++k) {
                    if (logagebins[k] <= logages[i] && logages[i] < logagebins[k+1]) {
                        // increase the age weight
                        logageweights[j][k] += 1.0;
                    }
                }
            }
        }
    }
    /*
    for (size_t i = 0; i < nummassbin; ++i) {
        cout << "m" << logmassweights[i] << endl;
        for (size_t j = 0; j < numagebin; ++j) {
            cout << logageweights[i][j] << endl;
        }
    }
    */
    // create distributions
    piecewise_constant_distribution<> logmassdist(logmassbins.begin(),
        logmassbins.end(),logmassweights.begin());
    vector<piecewise_constant_distribution<>> logagedists;
    for (size_t i = 0; i < nummassbin; ++i) {
        logagedists.push_back(piecewise_constant_distribution<>(logagebins.begin(),
        logagebins.end(),logageweights[i].begin()));
    }
    // vectors to store simulated mass
    vector<double> simlogmasses;
    vector<double> simlogages;
    // random number generator
    random_device rd;
    mt19937 gen(rd());
    for (size_t i = 0; i < 1000; ++i){
        double logmass = logmassdist(gen);
        double logage = 0.0;
        for (size_t i = 0; i < nummassbin; ++i) {
            if (logmassbins[i] <= logmass && logmass < logmassbins[i+1]) {
                logage = logagedists[i](gen);
            }
        }
        simlogmasses.push_back(logmass);
        simlogages.push_back(logage);
    }
    // pack the two vectors as a table
    vector<vector<double>> simstartable;
    simstartable.push_back(simlogmasses);
    simstartable.push_back(simlogages);
    return simstartable;
}

int main()
{     
    /*
    vector<vector<double>> cmktable = readcmk("cmkdata.txt");
    TTauriStar star = TTauriStar(cmktable, 0.68, 1, 1);
    star.update();
    string vectorname1 = "Age";
    string vectorname2 = "Period";
    vector<vector<double>> table = star.getvectors(vectorname1,vectorname2);
    
    FILE * temp = fopen("data.temp", "w");
    FILE* gp=popen("gnuplot -persistent","w");
    for(size_t k=0;k<table[0].size();k++) {
        fprintf(temp,"%f %f \n",table[0][k],table[1][k]);
    }
    fprintf(gp, "%s%s %s %s%s\n", "set title \"",vectorname2.data(),"vs",vectorname1.data(),"\"");    
    fprintf(gp, "%s \n", "set xlabel \"Age (Myr)\"");
    fprintf(gp, "%s \n", "set ylabel \"Period (days)\"");
    fprintf(gp, "%s \n", "plot 'data.temp'");
    */
    
    vector<vector<double>> startable = readcluster("dahm.txt");
    /*
    FILE * temp = fopen("data.temp", "w");
    FILE* gp=popen("gnuplot -persistent","w");
    for(size_t k=0;k<startable[0].size();k++) {
        fprintf(temp,"%f %f \n",startable[0][k],startable[1][k]);
    }
    fprintf(gp, "%s%s %s \n", "set title \"","ONC","\"");    
    fprintf(gp, "%s \n", "set xlabel \"LogMass\"");
    fprintf(gp, "%s \n", "set ylabel \"LogAge\"");
    fprintf(gp, "%s \n", "plot 'data.temp'");
    */
    vector<vector<double>> simstartable = generatedistribution(startable);
    FILE * temp = fopen("data.temp", "w");
    FILE* gp=popen("gnuplot -persistent","w");
    for(size_t k=0;k<simstartable[0].size();k++) {
        fprintf(temp,"%f %f \n",simstartable[0][k],simstartable[1][k]);
    }
    fprintf(gp, "%s%s %s \n", "set title \"","ONC","\"");    
    fprintf(gp, "%s \n", "set xlabel \"LogMass\"");
    fprintf(gp, "%s \n", "set ylabel \"LogAge\"");
    fprintf(gp, "%s \n", "plot 'data.temp'");
	return 0;
}
