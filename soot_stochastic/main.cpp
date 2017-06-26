//
//  main.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 22/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "mixing.hpp"
#include "tools.hpp"
#include "streams.hpp"
#include "nucleation.hpp"

using namespace std;

int main()
{
    cout << "Soot Stochastic" << endl << endl;
    
    //user inputs
    string pathProject("/Users/bouaniche/Xcode_projects/soot_stochastic");
    int Np0 = 100;                // initial stochastic particles at inlet0
    int Np1 = 100;                // initial stochastic particles at inlet1
    double c0 = 0.0;              // initial progress variable at inlet0
    double c1 = 1.0;              // initial progress variable at inlet1
    double l0 = 1.0;              // initial soot size at inlet0
    double l1 = 50.0;             // initial soot size at inlet1
    double lp0 = 1.0;             // nascent particles size
    int it = 10;                 // number of iteration
    double cSensitivity(0.05);    // distance between two c bins for graphic representation of P(c)
    double h = 3.0e5;               // constant used for source term of nucleation
    double nT0 = 1e5;             // initial total soot number density
    double deltaL = 0.5;           // size of I_\ell^* bins for Lpdf
    
    // time and mixing parameters
    double deltaT(1);             // iteration step time
    double tau(5);                // characteristic mixing time
    
    // initiate particles
    vector<vector<double> > allParticles;
    allParticles = initParticles(Np0, Np1, c0, c1, l0, l1);
     
    // Initial state. write in output files
    double t(0);
    double nT = nT0;
    writeCpdf(pathProject, "/outputs/Cpdf.dat", allParticles, cSensitivity);
    writeCpdft(pathProject, "/outputs/Cpdf_t/Cpdf", t, allParticles, cSensitivity);
    updateCpdf(pathProject, "/outputs/Cpdf.dat", t, allParticles, cSensitivity);
    
    // advancing t, mixing (Cpdf), source terms, advancing nT and Lpdf
    int j;
    for(j=0; j<it; j++ )
    {
        t = t+deltaT;                                   // advancing t
        mix(allParticles, deltaT, tau, t);              // advancing Cpdf = mixing
        double dotH = nuclSource(allParticles, h);      // calculation of nucleation source term
        nT = nT + dotH;                                 // advancing nT
                                                        // advancing Lpdf = reallocating soot particles
        
        cout << "nT = " << nT << "   " << "dotH = " << dotH << endl;
        
        LpdfAlphaH(allParticles, nT, dotH, deltaL, lp0, t);
        
        printParticles(allParticles, t);
        updateCpdf(pathProject, "/outputs/Cpdf.dat", t, allParticles, cSensitivity);
        writeCpdft(pathProject, "/outputs/Cpdf_t/Cpdf", t, allParticles, cSensitivity);
        
        // advancing nT, source terms...
    }
    
    
    return 0;
}
