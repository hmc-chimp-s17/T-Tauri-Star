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
//#include "ttauristar.hpp"

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
    const double SOLARTEMPERATURE = SOLARLUMINOSITY/pow(4*PI*pow(SOLARRADIUS,2)*SIGMAB,0.25);


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


int main()
{
    vector<vector<double>> cmktable = readcmk("cmkdata.txt");
    for (vector<double> vec : cmktable) {
        for (double value : vec) {
            cout << value << endl;
        }
    }
    for (vector<double> vec : cmktable) {
        cout << vec.size() << endl;
    }
	// simulateOneStar();
	return 0;
}